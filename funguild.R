#' Assign FungalTraits/FUNGuild fields to a phyloseq object + compute guild proportions
assign_funguild_to_phyloseq <- function(
    ps,
    traits_to_coalesce = c(
      "guild_fg",
      "trophic_mode_fg",
      "culture_media",
      "notes_fg",
      "source_funguild_fg",
      "speciesMatched",
      "confidence_fg"
    ),
    trait_cols = c(
      "guild_fg",
      "guild_source",
      "trophic_mode_fg",
      "notes_fg",
      "source_funguild_fg",
      "culture_media",
      "speciesMatched",
      "confidence_fg",
      "n_trait_rows_genus"
    ),
    mutualist_pattern  = "Ectomycorrhizal",
    saprotroph_pattern = "[S,s]aprotroph",
    pathogen_pattern   = "Plant [P,p]athogen|Fungal Parasite",
    exclude_pattern    = "Ectomycorrhizal",
    return_debug = FALSE
) {
  
  # ---- packages ----
  req_pkgs <- c("phyloseq","dplyr","stringr","tibble","fungaltraits","microbiome")
  missing <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing packages: ", paste(missing, collapse = ", "),
         "\nInstall them and try again.")
  }
  
  `%>%` <- dplyr::`%>%`
  coalesce <- dplyr::coalesce
  
  fung <- ps
  
  # pull genus/species from origional tax table (assumes Genus=col 6, Species=col 7)
  genera  <- fung@tax_table[, 6] %>% stringr::str_remove("^g__")
  species <- fung@tax_table[, 7] %>% stringr::str_remove("^s__")
  
  # download FungalTraits database
  traits_db <- fungaltraits::fungal_traits()
  
  # normalize db for joining
  traits_db_norm <- traits_db %>%
    dplyr::mutate(
      species = species %>%
        stringr::str_trim() %>%
        stringr::str_to_lower() %>%
        stringr::str_replace_all("\\s+", "_"),
      Genus = Genus %>% stringr::str_to_lower()
    )
  
  # ---- species-level assignment ----
  fungal_traits_sp <-
    data.frame(genus = genera) %>%
    dplyr::mutate(
      species = paste(genus, species, sep = "_") %>% stringr::str_to_lower(),
      genus   = stringr::str_to_lower(genus)
    ) %>%
    dplyr::filter(species != "na_na") %>%
    dplyr::distinct(species, .keep_all = TRUE) %>%
    dplyr::left_join(traits_db_norm, by = "species", multiple = "all")
  
  # ---- genus-level fallback ----
  fungal_traits_genus <-
    data.frame(genus = genera) %>%
    dplyr::mutate(
      Genus_key = stringr::str_trim(genus) %>% stringr::str_to_lower()
    ) %>%
    dplyr::filter(!stringr::str_detect(Genus_key, "gen_incertae_sedis")) %>%
    dplyr::distinct(Genus_key, .keep_all = TRUE) %>%
    dplyr::left_join(traits_db_norm, by = c("Genus_key" = "Genus"), multiple = "all")
  
  # collapse function for character columns with multiple assignments
  collapse_chars <- function(x) {
    x <- unique(stats::na.omit(x))
    if (length(x) == 0) NA_character_ else paste(sort(x), collapse = "|")
  }
  
  # ---- combine species- and genus-level traits ----
  traits_by_species <- fungal_traits_sp %>%
    dplyr::group_by(species, genus) %>%
    dplyr::summarise(
      n_trait_rows_species = dplyr::n(),
      dplyr::across(dplyr::where(is.numeric), ~ mean(.x, na.rm = TRUE)),
      dplyr::across(dplyr::where(is.character), collapse_chars),
      .groups = "drop"
    )
  
  traits_just_genus <- fungal_traits_genus %>%
    dplyr::group_by(Genus_key) %>%
    dplyr::summarise(
      n_trait_rows_genus = dplyr::n(),
      dplyr::across(dplyr::where(is.numeric), ~ mean(.x, na.rm = TRUE)),
      dplyr::across(dplyr::where(is.character), collapse_chars),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(speciesMatched))
  
  cols <- traits_to_coalesce
  cols_genus <- paste0(cols, ".genus")
  
  traits_one <- traits_by_species %>%
    dplyr::left_join(
      traits_just_genus %>% dplyr::select(Genus_key, n_trait_rows_genus, dplyr::all_of(cols)),
      by = c("genus" = "Genus_key"),
      suffix = c("", ".genus")
    ) %>%
    dplyr::mutate(
      guild_source = dplyr::case_when(
        !is.na(guild_fg) ~ "species",
        is.na(guild_fg) & !is.na(guild_fg.genus) ~ "genus",
        TRUE ~ "none"
      )
    ) %>%
    dplyr::mutate(
      dplyr::across(dplyr::all_of(cols), ~ dplyr::na_if(.x, ""))
    )
  
  # coalesce cols species vs genus
  traits_one[cols] <- Map(
    coalesce,
    traits_one[cols],
    traits_one[cols_genus]
  )
  
  # ---- finalize traits df ----
  traits_df <- traits_one %>%
    dplyr::select(-dplyr::ends_with(".genus")) %>%
    dplyr::select(species, genus, guild_source, n_trait_rows_species, n_trait_rows_genus,
                  dplyr::all_of(traits_to_coalesce))
  
  # ---- append traits to tax_table ----
  tt <- phyloseq::tax_table(fung) %>%
    as("matrix") %>%
    as.data.frame(stringsAsFactors = FALSE)
  
  taxon_id <- rownames(tt)
  
  tax_key <- tibble::tibble(
    taxon_id    = taxon_id,
    genus_raw   = stringr::str_remove(tt$Genus, "^g__") %>% stringr::str_trim(),
    species_raw = stringr::str_remove(tt$Species, "^s__") %>% stringr::str_trim()
  ) %>%
    dplyr::mutate(
      species_clean = dplyr::case_when(
        is.na(species_raw) | species_raw == "" ~ NA_character_,
        TRUE ~ species_raw
      ),
      species_key = paste(genus_raw, coalesce(species_clean, "na"), sep = "_") %>%
        stringr::str_to_lower()
    )
  
  tax_traits <- tax_key %>%
    dplyr::left_join(
      traits_df %>% dplyr::select(species, dplyr::any_of(trait_cols)),
      by = c("species_key" = "species")
    )
  
  tt2 <- tt
  for (nm in trait_cols) tt2[[nm]] <- tax_traits[[nm]]
  rownames(tt2) <- tax_traits$taxon_id
  
  fung_traits <- fung
  phyloseq::tax_table(fung_traits) <- phyloseq::tax_table(as.matrix(tt2))
  
  return(fung_traits)
}

a <- readRDS("test_data/fung_clean_physeq.RDS")

b <- assign_funguild_to_phyloseq(a)

c <- data.frame(tax_table(b))

                