library(ladybird)
library(here)
min_occurrences <- 1000
min_species <- 3
output_dir <- here("..", "Analysis", "models", "smoother")
maps <- here("..", "Analysis", "shape")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

for (country in c("BE", "NL")) {
  ds <- load_relevant(
    min_occurrences = min_occurrences, min_species = min_species,
    country = country
  )
  available_species <- head(tail(colnames(ds), -2), -1)
  rm(ds)
  gc()
  for (species in available_species) {
    output <- sprintf(
      "%s/%s_%s_%i_%i.rds", output_dir, tolower(species), tolower(country),
      min_occurrences, min_species
    )
    if (file.exists(output)) {
      next
    }
    message("smoother: ", country, " ", species, " ", Sys.time())
    bm <- try(smooth_model(
      species = species, min_occurrences = min_occurrences, country = country,
      min_species = min_species, path = maps, cellsize = 10e3
    ))
    if (!inherits(bm, "try-error")) {
      saveRDS(bm, output)
    }
    rm(bm)
    gc()
  }
}
