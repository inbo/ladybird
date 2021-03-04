#' Display a leaflet map with the occurrences for a given species
#' @export
#' @importFrom shiny shinyApp
occurrence_map <- function() {
  shinyApp(ui = occurrence_map_ui, server = occurrence_map_server)
}

#' @importFrom leaflet leafletOutput
#' @importFrom shiny absolutePanel bootstrapPage selectInput sliderInput tags
occurrence_map_ui <- function() {
  bootstrapPage(
    tags$style(type = "text/css", "html, body {width:100%;height:100%}"),
    leafletOutput("map", width = "100%", height = "100%"),
    title = "Occurrences",
    absolutePanel(
      bottom = 30, right = 10,
      sliderInput(
        "min_occurrences", label = "minimum occurrences", min = 100, max = 2500,
        step = 100, value = 1000
      ),
      sliderInput(
        "min_species", min = 1, max = 5, step = 1, value = 3,
        label = "minimum number of species per location-year combination"
      ),
      selectInput("select_species", label = "species", choices = ""),
      sliderInput(
        "year", min = 1991, max = 2020, step = 1, value = 1991,
        label = "year"
      )
    )
  )
}

#' @importFrom dplyr filter inner_join select %>%
#' @importFrom git2rdata read_vc
#' @importFrom leaflet addTiles addAwesomeMarkers awesomeIcons clearMarkers
#' fitBounds leaflet leafletProxy renderLeaflet
#' @importFrom rlang .data !!
#' @importFrom shiny observeEvent reactiveValues updateSelectInput
#' @importFrom utils tail
occurrence_map_server <- function(input, output, session) {
  data <- reactiveValues(
    icons = NULL,
    occurrences = NULL,
    selected_species = NULL,
    selected_year = NULL
  )

  observeEvent(
    input$min_occurrences, {
      data$occurrences <- load_relevant(
        min_occurrences = input$min_occurrences,
        min_species = input$min_species
      )
    }
  )

  observeEvent(
    input$min_species, {
      data$occurrences <- load_relevant(
        min_occurrences = input$min_occurrences,
        min_species = input$min_species
      )
    }
  )

  locations <- read_vc("location", system.file(package = "ladybird"))

  observeEvent(
    data$occurrences, {
      new_choices <- sort(tail(colnames(data$occurrences), -2))
      updateSelectInput(inputId = "select_species", choices = new_choices)
    }
  )

  observeEvent(
    input$select_species, {
      if (input$select_species == "") {
        return(NULL)
      }
      data$occurrences %>%
        inner_join(locations, by = "location") %>%
        select(
          .data$long, .data$lat, .data$year, occurrence = !!input$select_species
        ) -> data$selected_species
    }
  )

  observeEvent(
    data$selected_species, {
      data$selected_species %>%
        filter(.data$year == input$year) -> data$selected_year
      data$icons <- awesomeIcons(
        icon = "bug", iconColor = "black", library = "fa",
        markerColor = ifelse(
          data$selected_year$occurrence == 1, "blue", "lightgray"
        )
      )
    }
  )

  observeEvent(
    input$year, {
      if (is.null(data$selected_species)) {
        return(NULL)
      }
      data$selected_species %>%
        filter(.data$year == input$year) -> data$selected_year
      data$icons <- awesomeIcons(
        icon = "bug", iconColor = "black", library = "fa",
        markerColor = ifelse(
          data$selected_year$occurrence == 1, "blue", "lightgray"
        )
      )
    }
  )

  output$map <- renderLeaflet({
      leaflet() %>%
        addTiles() %>%
        fitBounds(
          min(locations$long), min(locations$lat),
          max(locations$long), max(locations$lat)
        )
    }
  )

  observeEvent(
    data$selected_year, {
      leafletProxy("map", data = data$selected_year) %>%
        clearMarkers() %>%
        addAwesomeMarkers(lng = ~long, lat = ~lat, icon = data$icons)
    }
  )
}
