library(shiny)
shinyUI(fluidPage(
  titlePanel("sRNA miner"),
  sidebarLayout(
    sidebarPanel(
      wellPanel(
        helpText(h4("Input:")),
        uiOutput("uibam"),        
        uiOutput("uichr")
      ),
      wellPanel(
        helpText(h4("Options:")),   
        sliderInput(
          inputId = "bg", 
          label = "Minimal Reads per End", 
          min = 5, 
          max = 200, 
          value = 40, 
          step = 1),
        sliderInput(
          inputId = "se", 
          label = "Sharpness of End",
          min = 0.30, 
          max = 1.00, 
          value = 0.75,
          step = 0.01),
        selectInput(
          inputId = "end", 
          label = "select End", 
          choices = list("5' end" = 1, "3' end" = 2, "5' or 3' end" = 3, "5' and 3' ends" = 4), 
          selected = 3)
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Help", 
          includeMarkdown("docs/Help.md")
        ),
        
        tabPanel(
          title = "Report", 
          wellPanel(
            helpText("It may takes several minutes to get the coverage infomation."),
            downloadButton("report", "Download Report")
          ),
          dataTableOutput("bed")
        )
      )
    )
  )
))
