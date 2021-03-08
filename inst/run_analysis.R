library(ladybird)
library(here)
min_occurrences <- 1000
min_species <- 3
output_dir <- here("..", "Analysis", "models", "base")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
ds <- load_relevant(
  min_occurrences = min_occurrences, min_species = min_species
)
available_species <- tail(colnames(ds), -2)

for (species in available_species) {
  output <- sprintf(
    "%s/%s_%i_%i.rds", output_dir, tolower(species), min_occurrences,
    min_species
  )
  if (file.exists(output)) {
    next
  }
  message("base: ", species)
  bm <- base_model(
    species = species, min_occurrences = min_occurrences,
    min_species = min_species
  )
  saveRDS(bm, output)
}
dir.create(file.path(output_dir, "..", "secondary"), showWarnings = FALSE)
bm <- readRDS(
  sprintf(
    "%s/%s_%i_%i.rds", output_dir, "harm_axyr", min_occurrences, min_species
  )
)
for (species in available_species) {
  if (species == "Harm_axyr") {
    next
  }
  output <- sprintf(
    "%s/%s_%i_%i.rds", file.path(output_dir, "..", "secondary"),
    tolower(species), min_occurrences, min_species
  )
  if (file.exists(output)) {
    next
  }
  message("secondary: ", species)
  sm <- secondary_model(
    species = species, min_occurrences = min_occurrences,
    min_species = min_species, secondary = bm
  )
  saveRDS(sm, output)
}
