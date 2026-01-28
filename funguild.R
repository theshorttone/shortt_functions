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
      "guild_fg", "guild_source", "trophic_mode_fg",
      "notes_fg", "source_funguild_fg", "culture_media"
    ),
    mutualist_pattern  = "Ectomycorrhizal",
    saprotroph_pattern = "[S,s]aprotroph",
    pathogen_pattern   = "Plant [P,p]athogen|Fungal Parasite",
    exclude_pattern    = "Ectomycorrhizal",
    return_debug = FALSE
) {
  
  # ---- packages (fail fast with clear messages) ----
  req_pkgs <- c("phyloseq","dplyr","stringr","tibble","fungaltraits","microbiome")
  missing <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing packages: ", paste(missing, collapse = ", "),
         "\nInstall them and try again.")
  }
  
  `%>%` <- dplyr::`%>%`
  coalesce <- dplyr::coalesce
  
  fung <- ps
  
  # pull genus/species from ORIGINAL tax table (assumes Genus=col 6, Species=col 7)
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
      n_trait_rows = dplyr::n(),
      dplyr::across(dplyr::where(is.numeric), ~ mean(.x, na.rm = TRUE)),
      dplyr::across(dplyr::where(is.character), collapse_chars),
      .groups = "drop"
    )
  
  traits_just_genus <- fungal_traits_genus %>%
    dplyr::group_by(Genus_key) %>%
    dplyr::summarise(
      n_trait_rows = dplyr::n(),
      dplyr::across(dplyr::where(is.numeric), ~ mean(.x, na.rm = TRUE)),
      dplyr::across(dplyr::where(is.character), collapse_chars),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(speciesMatched))
  
  cols <- traits_to_coalesce
  cols_genus <- paste0(cols, ".genus")
  
  traits_one <- traits_by_species %>%
    dplyr::left_join(
      traits_just_genus %>% dplyr::select(Genus_key, n_trait_rows, dplyr::all_of(cols)),
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
    dplyr::select(species, genus, guild_source, n_trait_rows, dplyr::all_of(traits_to_coalesce))
  
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
  
  # ---- guild vector from tax_table (named by taxa) ----
  tt_mat <- as(phyloseq::tax_table(fung_traits), "matrix")
  if ("guild_fg" %in% colnames(tt_mat)) {
    guild_col <- tt_mat[, "guild_fg"]
  } else {
    guild_col <- tt_mat[, 8]  # fallback (your original assumption)
  }
  guild_col <- as.character(guild_col)
  
  guild_by_taxon <- stats::setNames(guild_col, phyloseq::taxa_names(fung_traits))
  
  # ---- define which guild strings count ----
  mutualist_guilds <- grep(mutualist_pattern, guild_col, value = TRUE) %>% unique()
  
  saprotroph_guilds <-
    grep(saprotroph_pattern, guild_col, value = TRUE) %>%
    grep(pattern = exclude_pattern, x = ., value = TRUE, invert = TRUE) %>%
    unique()
  
  pathogen_guilds <-
    grep(pathogen_pattern, guild_col, value = TRUE) %>%
    grep(pattern = exclude_pattern, x = ., value = TRUE, invert = TRUE) %>%
    unique()
  
  # ---- counts -> proportions (NO relative abundance transform needed) ----
  lib_sizes <- phyloseq::sample_sums(fung)
  
  guild_props <- function(ps_counts, guild_by_taxon_named, guild_string_set) {
    tn <- phyloseq::taxa_names(ps_counts)
    this_guild <- guild_by_taxon_named[tn]  # align to this phyloseq's taxa order
    keep_taxa <- tn[!is.na(this_guild) & (this_guild %in% guild_string_set)]
    
    # if no taxa match, proportions are 0 for every sample
    if (length(keep_taxa) == 0) {
      out <- rep(0, phyloseq::nsamples(ps_counts))
      names(out) <- phyloseq::sample_names(ps_counts)
      return(out)
    }
    
    counts <- phyloseq::sample_sums(phyloseq::prune_taxa(keep_taxa, ps_counts))
    props  <- counts / phyloseq::sample_sums(ps_counts)
    
    # ensure alignment by sample name
    props <- props[phyloseq::sample_names(ps_counts)]
    props[is.na(props)] <- 0
    props
  }
  
  mutualist_proportions  <- guild_props(fung, guild_by_taxon, mutualist_guilds)
  saprotroph_proportions <- guild_props(fung, guild_by_taxon, saprotroph_guilds)
  pathogen_proportions   <- guild_props(fung, guild_by_taxon, pathogen_guilds)
  
  # ---- build guild_df ----
  mutualism_df <-
    microbiome::meta(fung_traits) %>%
    dplyr::mutate(proportion_mutualist = mutualist_proportions)
  
  saprotroph_df <-
    microbiome::meta(fung_traits) %>%
    dplyr::mutate(proportion_saprotroph = saprotroph_proportions)
  
  pathogen_df <-
    microbiome::meta(fung_traits) %>%
    dplyr::mutate(proportion_pathogen = pathogen_proportions)
  
  # keep your join style; make it robust by joining on shared columns
  guild_df <-
    mutualism_df %>%
    dplyr::full_join(saprotroph_df, by = intersect(names(mutualism_df), names(saprotroph_df))) %>%
    dplyr::full_join(pathogen_df,  by = intersect(names(mutualism_df), names(pathogen_df)))
  
  out <- list(
    ps_traits = fung_traits,
    guild_df  = guild_df,
    guild_sets = list(
      mutualist_guilds  = mutualist_guilds,
      saprotroph_guilds = saprotroph_guilds,
      pathogen_guilds   = pathogen_guilds
    )
  )
  
  if (return_debug) {
    out$traits_df  <- traits_df
    out$tax_traits <- tax_traits
  }
  
  return(out)
}
