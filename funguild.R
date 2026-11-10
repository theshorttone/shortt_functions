#' Assign FungalTraits (Polme et al. 2020) guilds to a phyloseq object
#'
#' Matches each taxon's genus against the FungalTraits genus-level database
#' (Polme et al. 2020, Fungal Diversity 105:1-16, doi:10.1007/s13225-020-00466-2)
#' and appends the trait columns to the tax_table.
#'
#' By default the database comes from the microeco package
#' (`microeco::fungi_func_FungalTraits`), which is identical to the published
#' FungalTraits spreadsheet. You can pass your own copy (e.g. the xlsx read
#' with readxl) via `traits_db` to pin a version.
#'
#' FungalTraits is genus-level only, so there is no species-level matching.
#' A few genera appear on more than one row with different values; for those,
#' the distinct values are pasted together with "|" and `ft_multiple_rows`
#' is TRUE so you can inspect them.
#'
#' Replaces the older version of this function that used
#' `fungaltraits::fungal_traits()` (Zanne et al. 2020 compilation; its FUNGuild
#' fields are a frozen ~2018 snapshot and matched far fewer genera).
#'
#' @param ps          phyloseq object with a genus rank in its tax_table.
#' @param genus_rank  name of the tax_table column holding genus (default "Genus").
#'                    Prefixes like "g__" are stripped.
#' @param trait_cols  FungalTraits columns to append (original database names).
#' @param traits_db   optional data.frame of FungalTraits; NULL loads it from microeco.
#' @return the phyloseq object with trait columns added to its tax_table. Added
#'   columns are lower-cased and prefixed "ft_", plus `ft_guild_simple`, a
#'   coarse grouping of primary_lifestyle (see `simplify_fungaltraits_guild()`).
assign_fungaltraits_to_phyloseq <- function(
    ps,
    genus_rank = "Genus",
    trait_cols = c(
      "primary_lifestyle",
      "Secondary_lifestyle",
      "Endophytic_interaction_capability_template",
      "Plant_pathogenic_capacity_template",
      "Decay_substrate_template",
      "Decay_type_template",
      "Growth_form_template",
      "Ectomycorrhiza_exploration_type_template",
      "Ectomycorrhiza_lineage_template"
    ),
    traits_db = NULL
) {

  # ---- packages ----
  req_pkgs <- c("phyloseq", "dplyr", "stringr", "tibble")
  if (is.null(traits_db)) req_pkgs <- c(req_pkgs, "microeco")
  missing <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing packages: ", paste(missing, collapse = ", "),
         "\nInstall them and try again.")
  }
  `%>%` <- dplyr::`%>%`

  # ---- load FungalTraits ----
  if (is.null(traits_db)) {
    env <- new.env()
    utils::data("fungi_func_FungalTraits", package = "microeco", envir = env)
    traits_db <- env$fungi_func_FungalTraits
    db_source <- paste0("FungalTraits (Polme et al. 2020) via microeco ",
                        utils::packageVersion("microeco"))
  } else {
    db_source <- "FungalTraits (Polme et al. 2020), user-supplied traits_db"
  }

  bad_cols <- setdiff(c("GENUS", trait_cols), names(traits_db))
  if (length(bad_cols)) {
    stop("Not found in the FungalTraits table: ", paste(bad_cols, collapse = ", "))
  }

  # collapse genera that appear on more than one row
  collapse_chars <- function(x) {
    x <- unique(stats::na.omit(x))
    if (length(x) == 0) NA_character_ else paste(sort(x), collapse = "|")
  }

  traits_genus <- traits_db %>%
    dplyr::select(GENUS, dplyr::all_of(trait_cols)) %>%
    dplyr::mutate(
      dplyr::across(dplyr::everything(), ~ dplyr::na_if(stringr::str_trim(as.character(.x)), "")),
      genus_key = stringr::str_to_lower(GENUS)
    ) %>%
    dplyr::filter(!is.na(genus_key)) %>%
    dplyr::group_by(genus_key) %>%
    dplyr::summarise(
      ft_multiple_rows = dplyr::n() > 1,
      dplyr::across(dplyr::all_of(trait_cols), collapse_chars),
      .groups = "drop"
    )

  new_names <- paste0("ft_", stringr::str_to_lower(stringr::str_remove(trait_cols, "_template$")))
  names(traits_genus)[match(trait_cols, names(traits_genus))] <- new_names

  # ---- match to tax_table ----
  tt <- as.data.frame(as(phyloseq::tax_table(ps), "matrix"), stringsAsFactors = FALSE)
  if (!genus_rank %in% names(tt)) {
    stop("tax_table has no column named '", genus_rank, "'. Ranks: ",
         paste(names(tt), collapse = ", "))
  }

  tax_key <- tibble::tibble(
    taxon_id  = rownames(tt),
    genus_key = tt[[genus_rank]] %>%
      stringr::str_remove("^g__") %>%
      stringr::str_trim() %>%
      stringr::str_to_lower()
  )

  tax_traits <- tax_key %>%
    dplyr::left_join(traits_genus, by = "genus_key") %>%
    dplyr::mutate(
      ft_guild_simple = simplify_fungaltraits_guild(ft_primary_lifestyle),
      ft_multiple_rows = dplyr::coalesce(ft_multiple_rows, FALSE)
    )
  stopifnot(identical(tax_traits$taxon_id, rownames(tt)))

  added <- c(new_names, "ft_guild_simple", "ft_multiple_rows")
  for (nm in added) tt[[nm]] <- as.character(tax_traits[[nm]])

  out <- ps
  phyloseq::tax_table(out) <- phyloseq::tax_table(as.matrix(tt))

  n_match <- sum(!is.na(tax_traits$ft_primary_lifestyle))
  message(db_source, ": primary_lifestyle assigned to ", n_match, " of ",
          nrow(tax_traits), " taxa (", round(100 * n_match / nrow(tax_traits), 1), "%).")
  out
}

#' Coarse guild groups from FungalTraits primary_lifestyle
#'
#' ectomycorrhizal, arbuscular_mycorrhizal, saprotroph (any *_saprotroph),
#' plant_pathogen, endophyte (foliar/root), lichenized, other (any other
#' assigned lifestyle), NA (no match). Genera with multiple primary lifestyles
#' ("a|b") are classed by their first listed value.
simplify_fungaltraits_guild <- function(primary_lifestyle) {
  x <- sub("\\|.*$", "", primary_lifestyle)
  dplyr::case_when(
    is.na(x)                                ~ NA_character_,
    x == "ectomycorrhizal"                  ~ "ectomycorrhizal",
    x == "arbuscular_mycorrhizal"           ~ "arbuscular_mycorrhizal",
    grepl("saprotroph$", x)                 ~ "saprotroph",
    x == "plant_pathogen"                   ~ "plant_pathogen",
    x %in% c("foliar_endophyte", "root_endophyte") ~ "endophyte",
    x == "lichenized"                       ~ "lichenized",
    TRUE                                    ~ "other"
  )
}

# ---- example (runs only with `Rscript funguild.R`, not when source()d) ----
if (sys.nframe() == 0) {
  a <- readRDS("test_data/fung_clean_physeq.RDS")
  b <- assign_fungaltraits_to_phyloseq(a)
  c <- data.frame(phyloseq::tax_table(b))
  print(table(c$ft_guild_simple, useNA = "ifany"))
}
