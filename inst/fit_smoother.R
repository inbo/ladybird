library(ladybird)
library(here)
min_occurrences <- 1000
min_species <- 3
buffer_distance <- 50e3
buffer_locations <- 50
output_dir <- here("..", "Analysis", "models", "smoother")
maps <- here("..", "Analysis", "shape")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

for (country in c("BE", "NL", "GB")) {
  ds <- load_relevant(
    min_occurrences = min_occurrences, min_species = min_species,
    country = country, buffer_distance = buffer_distance,
    buffer_locations = buffer_locations
  )
  available_species <- sample(head(tail(colnames(ds), -3), -1))
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
      min_species = min_species, path = maps, buffer_distance = buffer_distance,
      cellsize = ifelse(country == "GB", 15e3, 10e3),
      buffer_locations = buffer_locations
    ))
    if (!inherits(bm, "try-error")) {
      saveRDS(bm, output)
    }
    rm(bm)
    gc()
  }
}
