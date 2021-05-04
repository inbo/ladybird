library(ladybird)
library(here)
prepare_smoother_results(
  models = here("..", "Analysis", "models", "smoother"),
  maps = here("..", "Analysis", "shape"), output = here("..", "Analysis")
)
