library(ladybird)
library(here)
min_occurrences <- 1000
min_species <- 3
first_order <- TRUE
center_year <- 2001
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
  message("base: ", species, " ", first_order, " ", Sys.time())
  bm <- base_model(
    species = species, min_occurrences = min_occurrences,
    min_species = min_species, first_order = first_order,
    center_year = center_year
  )
  saveRDS(bm, output)
}
bm <- readRDS(
  sprintf(
    "%s/%s_%i_%i_%i.rds", output_dir, "harm_axyr", min_occurrences, min_species,
    first_order
  )
)

dir.create(file.path(output_dir, "..", "probability"), showWarnings = FALSE)
for (species in available_species) {
  if (species == "Harm_axyr") {
    next
  }
  output <- sprintf(
    "%s/%s_%i_%i_%i.rds", file.path(output_dir, "..", "probability"),
    tolower(species), min_occurrences, min_species, first_order
  )
  if (file.exists(output)) {
    next
  }
  message("probability: ", species, " ", first_order, " ", Sys.time())
  sm <- probability_model(
    species = species, min_occurrences = min_occurrences,
    min_species = min_species, secondary = bm, first_order = first_order,
    center_year = center_year
  )
  saveRDS(sm, output)
}

dir.create(file.path(output_dir, "..", "cumulative"), showWarnings = FALSE)
for (species in available_species) {
  if (species == "Harm_axyr") {
    next
  }
  output <- sprintf(
    "%s/%s_%i_%i_%i.rds", file.path(output_dir, "..", "cumulative"),
    tolower(species), min_occurrences, min_species, first_order
  )
  if (file.exists(output)) {
    next
  }
  message("cumulative: ", species, " ", first_order, " ", Sys.time())
  sm <- cumulative_model(
    species = species, min_occurrences = min_occurrences,
    min_species = min_species, secondary = bm, first_order = first_order,
    center_year = center_year
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
  message("base: ", species, " ", first_order, " ", Sys.time())
  bm <- base_model(
    species = species, min_occurrences = min_occurrences,
    min_species = min_species, first_order = first_order,
    center_year = center_year
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
    "%s/%s_%i_%i_%i.rds", file.path(output_dir, "..", "probability"),
    tolower(species), min_occurrences, min_species, first_order
  )
  if (file.exists(output)) {
    next
  }
  message("probability: ", species, " ", first_order, " ", Sys.time())
  sm <- probability_model(
    species = species, min_occurrences = min_occurrences,
    min_species = min_species, secondary = bm, first_order = first_order,
    center_year = center_year
  )
  saveRDS(sm, output)
}

dir.create(file.path(output_dir, "..", "cumulative"), showWarnings = FALSE)
for (species in available_species) {
  if (species == "Harm_axyr") {
    next
  }
  output <- sprintf(
    "%s/%s_%i_%i_%i.rds", file.path(output_dir, "..", "cumulative"),
    tolower(species), min_occurrences, min_species, first_order
  )
  if (file.exists(output)) {
    next
  }
  message("cumulative: ", species, " ", first_order, " ", Sys.time())
  sm <- try(cumulative_model(
    species = species, min_occurrences = min_occurrences,
    min_species = min_species, secondary = bm, first_order = first_order,
    center_year = center_year
  ))
  if (!inherits(sm, "try-error")) {
    saveRDS(sm, output)
  }
}
