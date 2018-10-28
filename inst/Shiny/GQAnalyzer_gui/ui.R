
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(GQAnalyzer)
library(DT)
library(shinyalert)
#
## ui for GQAnalyzer
#
shinyUI(
  navbarPage(
    "GQAnalyzer",
    tabPanel("Input Data",
             sidebarLayout(
               sidebarPanel(width = 3,
                            useShinyalert(),
                            imageOutput("uptc.logo", inline=TRUE),
                            p(HTML("<h5>This is GQAnalyzer-GUI, the Shiny interface for analysis and
                                    plotting and analysis of Groundwater Quality and Hydrochemical
                                    characteristics of geothermal waters in <strong>R</strong>.</h5>
                                    This application can be used for the identification of the type of
                                    Groundwater, identification of hydrogeochemical groups, comparison of
                                    Groundwater Quality samples, and statistical analysis of Groundwater Quality
                                    data in space and time using the  <strong>R</strong>
                                    package  <a href='http://www.github.com/khaors/GQAnalyzer'>GQAnalyzer</a>.
                                    The analysis of a pumping/slug test is achieved in simple steps using the
                                    panels on the right."))),
               #
               mainPanel(
                 fluidRow(
                   column(width = 4,
                          h3("GQAnalyzer Import"),
                          br(),
                          br(),
                          #textInput(inputId = "project.name", label = "Project Name",
                          #          value = "test.GQp"),
                          #br(),
                          checkboxInput('header', ' Header?', TRUE),
                          checkboxInput('rownames', ' Row names?', FALSE),
                          selectInput('sep', 'Separator:',
                                      c("Comma","Semicolon","Tab","Space"), 'Comma'),
                          selectInput('quote', 'Quote:',
                                      c("None","Double Quote","Single Quote"),
                                      'Double Quote'),
                          selectInput('dec', 'Decimal mark', c("Period", "Comma"),
                                      'Period'),
                          numericInput('nrow.preview','Number of rows in the preview:',20),
                          numericInput('ncol.preview', 'Number of columns in the preview:',
                                       10),
                          fileInput('file1', 'Choose CSV/TXT File')
                          ),#end column
                   column(width = 8,
                          h3("Geochemical Dataset"),
                          br(),
                          helpText("Note: Even if the preview only shows a restricted
                          number of observations, the pumping_test object will be created
                          based on the full dataset."),
                          dataTableOutput(outputId = "view")
                          )#,
                   #column(width = 4,
                   #      h3("Statistical Summaries"),
                   #      br(),
                   #       verbatimTextOutput("summary"))
                 )#end fluid row

               )
             )
    ),
    #
    tabPanel(title = "Create Geochemical Dataset",
             sidebarLayout(
               sidebarPanel(width = 4,
                            p(HTML("<h5>The geochemical dataset is created in this tab</h5>"))),
               mainPanel(width = 8,
                         p(HTML("<h4>Fit Columns by names</h4>")),
                         br(),
                         actionButton(inputId = "col.fit", label = "Fit columns"),
                         br(),
                         p(HTML("<h4>Choose the columns in the data table corresponding to the
                                major ions.</h4><br><br>")),
                         h4("Physical Properties"),
                         br(),
                         fluidRow(
                           column(width = 3, uiOutput(outputId = "col.pH")),
                           column(width = 3, uiOutput(outputId = "col.TDS")),
                           column(width = 3, uiOutput(outputId = "col.EC")),
                           column(width = 3, uiOutput(outputId = "col.Temp"))),
                         br(),
                         h4("Cations"),
                         br(),
                         br(),
                         fluidRow(
                           column(width = 3, uiOutput(outputId = "col.Ca")),
                           column(width = 3, uiOutput(outputId = "col.Mg")),
                           column(width = 3, uiOutput(outputId = "col.Na")),
                           column(width = 3, uiOutput(outputId = "col.K"))),
                         br(),
                         h4("Anions"),
                         br(),
                         br(),
                         fluidRow(
                           column(width = 3, uiOutput(outputId = "col.HCO3")),
                           column(width = 3, uiOutput(outputId = "col.CO3")),
                           column(width = 3, uiOutput(outputId = "col.Cl")),
                           column(width = 3, uiOutput(outputId = "col.SO4"))),
                         br(),
                         actionButton("create.gdata", label = "Define Geochemical Dataset",
                                      icon = icon("bullseye"))
                         ))),
    #tabPanel(title = "Transformation"),
    #navbarMenu(title = "Transformation",
    #           tabPanel("Trans1"),
    #           tabPanel("Trans2"),
    #           tabPanel("Trans3")),
    #
    tabPanel(title = "EDA",
             sidebarLayout(
               sidebarPanel(
                 p(HTML("<h4>Create Statistical Plots</h4>")),
                 br(),
                 uiOutput(outputId = "eda.varselector"),
                 br(),
                 radioButtons(inputId = "eda.log", label = "Log",
                              choices = c("No", "Yes"),
                              selected = "No"),
                 br(),
                 uiOutput(outputId = "eda.nbins")
               ),
               mainPanel(
                 plotOutput(outputId = "eda.plot", height = "600px", width = "80%")
                 #br(),
                 #downloadButton('eda.downloadPlot', 'Download Plot')
               )
             )
    ),
    tabPanel(title = "CrossPlots",
             sidebarLayout(
               sidebarPanel(
                 p(HTML("<h4>Create Scatterplots</h4>")),
                 br(),
                 p(HTML("<h3>X Variable</h3>")),
                 br(),
                 uiOutput(outputId = "cross.varselector1"),
                 radioButtons(inputId = "cross.log1", label = "Log",
                              choices = c("No", "Yes"),
                              selected = "No"),
                 br(),
                 p(HTML("<h3>Y Variable</h3>")),
                 uiOutput(outputId = "cross.varselector2"),
                 radioButtons(inputId = "cross.log2", label = "Log",
                              choices = c("No", "Yes"),
                              selected = "No"),
                 br(),
                 uiOutput(outputId = "cross.varsize"),
                 uiOutput(outputId = "cross.varcolor")
               ),
               mainPanel(
                 plotOutput(outputId = "cross.plot", height = "600px", width = "80%")
               )
             )
    ),
    #navbarMenu(title = "EDA",
    #           tabPanel("EDA1"),
    #           tabPanel("EDA2"),
    #           tabPanel("EDA3")),
    #
    tabPanel(title = "Hydrogeochemical Plots",
             sidebarLayout(
               sidebarPanel(
                 h3("Choose Plot"),
                 br(),
                 uiOutput(outputId = "hplot.tselector"),
                 br(),
                 uiOutput(outputId = "hplot.tselector1"),
                 br(),
                 uiOutput(outputId = "hplot.option1"),
                 uiOutput(outputId = "hplot.option2"),
                 uiOutput(outputId = "hplot.option3"),
                 uiOutput(outputId = "hplot.option4"),
                 uiOutput(outputId = "hplot.option5")
               ),
               mainPanel(
                 plotOutput("hplot", height = "600px", width = "80%")
               )
             )),
    #tabPanel(title = "Maps"),
    #navbarMenu(title = "Maps",
    #           tabPanel("Piper"),
    #           tabPanel("Circular")),
    #
    #tabPanel(title = "Mixing Utilities"),
    #navbarMenu(title = "Mixing Utilities",
    #           tabPanel("M3 Mixing Model"),
    #           tabPanel("MIX Model"),
    #           tabPanel("Component Mixing Model")),
    #
    #tabPanel(title = "Anomaly Utilities"),
    #navbarMenu(title = "Anomaly Utilities",
    #           tabPanel("Statistical Methods"),
    #           tabPanel("Graphical Methods"),
    #           tabPanel("Concentration-Area Method"),
    #           tabPanel("Geostatistical Methods")),
    #
    tabPanel(title = "Report")
    #navbarMenu(title = "Report",
    #           tabPanel("Options"),
    #           tabPanel("Generation")),
    #navbarMenu(title = "Help",
    #           tabPanel("Contents"),
    #           tabPanel("About"))
  )
)

