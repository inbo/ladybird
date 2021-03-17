library(ladybird)
library(here)
min_occurrences <- 1000
min_species <- 3
output_dir <- here("..", "Analysis", "models", "base")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
ds <- load_relevant(
  min_occurrences = min_occurrences, min_species = min_species
)
available_species <- head(tail(colnames(ds), -2), -1)

for (species in available_species) {
  output <- sprintf(
    "%s/%s_%i_%i.rds", output_dir, tolower(species), min_occurrences,
    min_species
  )
  if (file.exists(output)) {
    next
  }
  message("base: ", species, " ", Sys.time())
  bm <- try(base_model(
    species = species, min_occurrences = min_occurrences,
    min_species = min_species, knots = c(1990, 2000, 2010, 2020)
  ))
  if (!inherits(bm, "try-error")) {
    saveRDS(bm, output)
  }
}
bm <- readRDS(
  sprintf(
    "%s/%s_%i_%i.rds", output_dir, "harm_axyr", min_occurrences, min_species
  )
)

dir.create(file.path(output_dir, "..", "probability"), showWarnings = FALSE)
for (species in available_species) {
  if (species == "Harm_axyr") {
    next
  }
  output <- sprintf(
    "%s/%s_%i_%i.rds", file.path(output_dir, "..", "probability"),
    tolower(species), min_occurrences, min_species
  )
  if (file.exists(output)) {
    next
  }
  message("probability: ", species, " ", Sys.time())
  sm <- try(probability_model(
    species = species, min_occurrences = min_occurrences,
    min_species = min_species, secondary = bm, knots = c(1990, 2000, 2010, 2020)
  ))
  if (!inherits(sm, "try-error")) {
    saveRDS(sm, output)
  }
}

dir.create(file.path(output_dir, "..", "cumulative"), showWarnings = FALSE)
for (species in available_species) {
  if (species == "Harm_axyr") {
    next
  }
  output <- sprintf(
    "%s/%s_%i_%i.rds", file.path(output_dir, "..", "cumulative"),
    tolower(species), min_occurrences, min_species
  )
  if (file.exists(output)) {
    next
  }
  message("cumulative: ", species, " ", Sys.time())
  sm <- cumulative_model(
    species = species, min_occurrences = min_occurrences,
    min_species = min_species, secondary = bm, knots = c(1990, 2000, 2010, 2020)
  )
  saveRDS(sm, output)
}
