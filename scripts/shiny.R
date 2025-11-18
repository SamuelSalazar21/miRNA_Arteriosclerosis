

#tutorial de shiny
# install.packages("shiny")
# library(shiny)
# library(bslib)
# # runExample("01_hello")
# # runApp("App-1")
# # dir()
# # 
# # runApp("App-1", display.mode = "showcase")
# 
# 
# 
# 
# ui <- page_sidebar(
#   title = "title panel",
#   sidebar=sidebar("sidebar"),
#   "main contents"
# )
# server <- function(input, output){
# 
# }
# 
# 
# # file.copy("C:/Users/PC/Documents/www/example.jpg", "example.jpg", overwrite = T)
# # 
# # ui <- page_sidebar(
# #   title="title panel",
# #   sidebar=sidebar("Sidebar"),
# #   card(
# #     card_header("Card header"),
# #     "cuerpo"),
# #   card_image(
# #     file = "example.jpg"
# #   
# #   )
# # )
# 
# 
# 
# ui <- page_sidebar(
#   title="titulo",
#   sidebar=sidebar("al lado"),
#   value_box(
#     tittle="caja de valores",
#     value=100,
#     showcase=bsicons::bs_icon("bar-chart"),
#     theme="teal"
#   ),
#   card("carta")
# )
# 
# server <- function(input, output){
# 
#   }
# 
# shinyApp(ui=ui, server=server)
# 
# 
# 
# 
# ui <- page_sidebar(
#   title = "title panel",
#   sidebar = sidebar("Sidebar"),
#   value_box(
#     title = "Value box",
#     value = 100,
#     showcase = bsicons::bs_icon("bar-chart"),
#     theme = "teal"
#   ),
#   card("Card"),
#   card("Another card")
# )
# getwd()
# 
# 
# 
# 
# 
# 
# ### ejercicio
# 
# 
# 
# file.copy("C:/Users/PC/Downloads/shiny-thumb.png", "shinyimg.png", overwrite=T)
# 
# ui <- page_sidebar(
#   title = "My shiny app",
#   sidebar = sidebar("cran", code("getwd()")), code("dir()"),
# card(
#   card_header("Introducing Shiny"),
#   card_image(file="shinyimg.png"),
#   card_footer("Shiny is a product of Posit"),
#   
# ))

# # 
# ui <- page_sidebar(
#   title = "censusVis",
#   sidebar = sidebar(
#     helpText(
#       "Create demographic maps with information from the 2010 US Census."
#     ),
#     selectInput(
#       "var",
#       label = "Choose a variable to display",
#       choices = 
#         c("Percent White",
#           "Percent Black",
#           "Percent Hispanic",
#           "Percent Asian"),
#       selected = "Percent White"
#     ),
#     sliderInput(
#       "range",
#       label = "Range of interest:",
#       min = 0, 
#       max = 100, 
#       value = c(0, 100)
#     )
#   ),
#   textOutput("selected_var"),
#   textOutput("selected_range")
# )
# 
# 
# server <- function(input, output){
#   output$selected_var <- renderText({paste(
#     "you have selected this", input$var)
#   })
#   output$selected_range <- renderText({paste("usted ha seleccionado", min( input$range), "y", max(input$range))})
#   
# }
# 
# shinyApp(ui=ui, server=server)
# getwd()
# setwd("C:/Users/PC/Documents")
# 
# 
# 
# condados <- readRDS("C:/Users/PC/Documents/App-1/data/counties.rds")
# head(condados)
# 
# 
# library(maps)
# library(mapproj)
# 
# 
# source("C:/Users/PC/Documents/App-1/helpers.R")
# percent_map(condados$white, "darkgreen", "% White")
# 
# 
# 
# 
# ui <-page_fluid(
#   textInput(
#     inputId = "custom_text",
#     label = "input some text here"
#   ),
#   strong("text is shown below:"),
#   
#   textOutput(outputId="texto")
# )
#   
#   server <- function(input, output, session){
#     output$texto <- renderText({input$custom_text})
#   }
#   
#   shinyApp(ui=ui, server=server)
  

# Load packages ----------------------------------------------------------------

library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(DT)

# Load data --------------------------------------------------------------------

load("movies.RData")

# Define UI --------------------------------------------------------------------

n_total=nrow(movies)
all_movies=sort(unique(movies$studio))


ui <- page_sidebar(
  
  sidebar = sidebar(
    
    HTML(paste("Enter a value between 1 and", n_total)),
    
    numericInput(inputId = "n",
                 label = "select a number",
                 value = 30,
                 step = 1,
                 min=1,
                 max=n_total),
    selectInput(
      inputId = "Studio",
      label="please select studio",
      choices = all_movies,
      selected = NULL,
      multiple = TRUE
      
    )
    
  ),
  
  card(
    DT::dataTableOutput(outputId = "moviestable")
  )
)

# Define server ----------------------------------------------------------------

server <- function(input, output, session) {
  
  
  output$moviestable <- DT::renderDataTable({
    req(input$n)
    movies_sample <- movies %>%
      sample_n(input$n) %>%
      select(title:studio)
    DT::datatable(data = movies_sample,
                  options = list(pageLength = 10),
                  rownames = FALSE)
  })
  
}

# Create the Shiny app object --------------------------------------------------

shinyApp(ui = ui, server = server)
  





