library(ladybird)
library(here)
min_occurrences <- 1000
min_species <- 3
first_order <- TRUE
output_dir <- here("..", "Analysis", "models", "base")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
ds <- load_relevant(
  min_occurrences = min_occurrences, min_species = min_species
)
available_species <- tail(colnames(ds), -2)

for (species in available_species) {
  output <- sprintf(
    "%s/%s_%i_%i_%i.rds", output_dir, tolower(species), min_occurrences,
    min_species, first_order
  )
  if (file.exists(output)) {
    next
  }
  message("base: ", species)
  bm <- base_model(
    species = species, min_occurrences = min_occurrences,
    min_species = min_species, first_order = first_order
  )
  saveRDS(bm, output)
}
dir.create(file.path(output_dir, "..", "secondary"), showWarnings = FALSE)
bm <- readRDS(
  sprintf(
    "%s/%s_%i_%i_%i.rds", output_dir, "harm_axyr", min_occurrences, min_species,
    first_order
  )
)
for (species in available_species) {
  if (species == "Harm_axyr") {
    next
  }
  output <- sprintf(
    "%s/%s_%i_%i_%i.rds", file.path(output_dir, "..", "secondary"),
    tolower(species), min_occurrences, min_species, first_order
  )
  if (file.exists(output)) {
    next
  }
  message("secondary: ", species)
  sm <- secondary_model(
    species = species, min_occurrences = min_occurrences,
    min_species = min_species, secondary = bm, first_order = first_order
  )
  saveRDS(sm, output)
}

first_order <- FALSE
for (species in available_species) {
  output <- sprintf(
    "%s/%s_%i_%i_%i.rds", output_dir, tolower(species), min_occurrences,
    min_species, first_order
  )
  if (file.exists(output)) {
    next
  }
  message("base: ", species)
  bm <- base_model(
    species = species, min_occurrences = min_occurrences,
    min_species = min_species, first_order = first_order
  )
  saveRDS(bm, output)
}
bm <- readRDS(
  sprintf(
    "%s/%s_%i_%i_%i.rds", output_dir, "harm_axyr", min_occurrences, min_species,
    first_order
  )
)
for (species in available_species) {
  if (species == "Harm_axyr") {
    next
  }
  output <- sprintf(
    "%s/%s_%i_%i_%i.rds", file.path(output_dir, "..", "secondary"),
    tolower(species), min_occurrences, min_species, first_order
  )
  if (file.exists(output)) {
    next
  }
  message("secondary: ", species)
  sm <- secondary_model(
    species = species, min_occurrences = min_occurrences,
    min_species = min_species, secondary = bm, first_order = first_order
  )
  saveRDS(sm, output)
}
