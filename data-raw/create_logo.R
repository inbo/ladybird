checklist::create_hexsticker(
  package_name = "ladybird", scale = 0.45, x = 144, y = -130,
  filename = file.path("man", "figures", "logo.svg"),
  icon = file.path("data-raw", "ladybird.svg")
)
pkgdown::build_favicons(overwrite = TRUE)
