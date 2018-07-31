#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(DT)
library(ggplot2)
library(gridExtra)
#
plot.types.all <-  c("None", "ternary", "piper", "modified_piper", "durov", "schoeller",
                     "multirectangular")
#
plot.types.single <- c("None", "stiff", "radial")
#
measures.type <- c("None", "conc", "meql")
#
shinyServer(function(input, output, session) {
  #
  server.env <- environment() # used to allocate in functions
  current.table <- NULL
  current.names <- NULL
  current.gdata <- NULL # this variable will contain the current geochemical dataset
  current.geonames <- NULL
  current.plot <- NULL
  current.names <- NULL
  first <- TRUE
  #
  output$uptc.logo <- renderImage(list(src = "uptc_jpg.jpg"),
                                  deleteFile = FALSE)
  #
  ## Panel 'Import data'
  dInput <- reactive({
    in.file <- input$file1
    #
    validate(
      need(input$file1, 'Check if file is loaded')
    )
    #
    if (is.null(in.file))
      return(NULL)
    #
    fname <- strsplit(input$file1$name, "\\.")
    the.sep <- switch(input$sep, "Comma"=",", "Semicolon"=";", "Tab"="\t",
                      "Space"="")
    the.quote <- switch(input$quote, "None"="","Double Quote"='"',
                        "Single Quote"="'")
    the.dec <- switch(input$dec, "Period"=".", "Comma"=",")
    if (input$rownames) {
      the.table <- read.table(in.file$datapath, header=input$header,
                              sep=the.sep, quote=the.quote, row.names=1,
                              dec=the.dec)
    } else {
      the.table <- read.table(in.file$datapath, header=input$header,
                              sep=the.sep, quote=the.quote, dec=the.dec)
    }
    if(!first){

    }
    if(first)
      first <- FALSE
    server.env$first <- first
    # return the table
    server.env$current.table <- the.table
    #print("Original Names")
    #print(names(the.table))
    server.env$current.names <- names(the.table)
    the.table
  })
  #
  # data preview table
  output$view <- renderDataTable({
    d.input <- dInput()
    if (is.null(d.input))
      return(NULL)
    if (ncol(d.input)>input$ncol.preview)
      d.input <- d.input[,1:input$ncol.preview]
    head(d.input, n=input$nrow.preview)
  },
  extensions = c('Buttons'),
  options = list(
    autoWidth = TRUE,
    pageLength = input$nrow.preview,
    dom = 'Bfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
    text = 'Download',
    scrollY = 200,
    scrollX = TRUE,
    scroller = TRUE))
  #
  output$summary <- renderPrint({
    d.input <- dInput()
    if (is.null(d.input))
      return(NULL)
    the.table <- server.env$current.table
    if(is.null(the.table))
      return(NULL)
    summary(the.table)
  })
  #######################################################################################
  #                                Create geochemical dataset tab
  #######################################################################################
  output$col.Ca <- renderUI({
    res <- NULL
    if(is.null(server.env$current.table))
      return(NULL)
    else{
      all.variables <- c("None", server.env$current.names)
      res <- selectInput(inputId = "col.Ca1", label = "Ca Column", choices = all.variables,
                         width = 100)
    }
    return(res)
  })
  #
  output$col.Mg <- renderUI({
    res <- NULL
    if(is.null(server.env$current.table)){
      return(NULL)
    }
    else{
      all.variables <- c("None", server.env$current.names)
      res <- selectInput(inputId = "col.Mg1", label = "Mg Column", choices = all.variables,
                         width = 100)
    }
    return(res)
  })
#
  output$col.Na <- renderUI({
    res <- NULL
    if(is.null(server.env$current.table)){
      return(NULL)
    }
    else{
      all.variables <- c("None", server.env$current.names)
      res <- selectInput(inputId = "col.Na1", label = "Na Column", choices = all.variables,
                         width = 100)
    }
    return(res)
  })
#
  output$col.K <- renderUI({
    res <- NULL
    if(is.null(server.env$current.table)){
      return(NULL)
    }
    else{
      all.variables <- c("None", server.env$current.names)
      res <- selectInput(inputId = "col.K1", label = "K Column", choices = all.variables,
                         width = 100)
    }
    return(res)
  })
  #
  output$col.HCO3 <- renderUI({
    res <- NULL
    if(is.null(server.env$current.table)){
      return(NULL)
    }
    else{
      all.variables <- c("None", server.env$current.names)
      res <- selectInput(inputId = "col.HCO31", label = "HCO3 Column", choices = all.variables,
                         width = 100)
    }
    return(res)
  })
  #
  output$col.CO3 <- renderUI({
    res <- NULL
    if(is.null(server.env$current.table)){
      return(NULL)
    }
    else{
      all.variables <- c("None", server.env$current.names)
      res <- selectInput(inputId = "col.CO31", label = "CO3 Column", choices = all.variables,
                         width = 100)
    }
    return(res)
  })
  #
  output$col.Cl <- renderUI({
    res <- NULL
    if(is.null(server.env$current.table)){
      return(NULL)
    }
    else{
      all.variables <- c("None", server.env$current.names)
      res <- selectInput(inputId = "col.Cl1", label = "Cl Column", choices = all.variables,
                         width = 100)
    }
    return(res)
  })
  #
  output$col.SO4 <- renderUI({
    res <- NULL
    if(is.null(server.env$current.table)){
      return(NULL)
    }
    else{
      all.variables <- c("None", server.env$current.names)
      res <- selectInput(inputId = "col.SO41", label = "SO4 Column", choices = all.variables,
                         width = 100)
    }
    return(res)
  })
#
  output$col.pH <- renderUI({
    res <- NULL
    if(is.null(server.env$current.table)){
      return(NULL)
    }
    else{
      all.variables <- c("None", server.env$current.names)
      res <- selectInput(inputId = "col.pH1", label = "pH Column", choices = all.variables,
                         width = 100)
    }
    return(res)
  })
  #
  output$col.TDS <- renderUI({
    res <- NULL
    if(is.null(server.env$current.table)){
      return(NULL)
    }
    else{
      all.variables <- c("None", server.env$current.names)
      res <- selectInput(inputId = "col.TDS1", label = "TDS Column", choices = all.variables,
                         width = 100)
    }
    return(res)
  })
  #
  output$col.EC <- renderUI({
    res <- NULL
    if(is.null(server.env$current.table)){
      return(NULL)
    }
    else{
      all.variables <- c("None", server.env$current.names)
      res <- selectInput(inputId = "col.EC1", label = "Elec.Cond. Column",
                         choices = all.variables, width = 100)
    }
    return(res)
  })
  #
  fit_columns <- function(){
    current.names <- server.env$current.names
    #print("FIT")
    #print(current.names)
    if(is.null(current.names))
      return(NULL)
    pos <- current.names == "Ca"
    #print(pos)
    if(sum(pos) == 1){
      updateSelectInput(session, inputId = "col.Ca1", selected = "Ca")
    }
    else{
      updateSelectInput(session, inputId = "col.Ca1", selected = "None")
    }
    #
    pos <- current.names == "Mg"
    if(sum(pos) == 1){
      updateSelectInput(session, inputId = "col.Mg1", selected = "Mg")
    }
    else{
      updateSelectInput(session, inputId = "col.Mg1", selected = "None")
    }
    #
    pos <- current.names == "Na"
    if(sum(pos) == 1){
      updateSelectInput(session, inputId = "col.Na1", selected = "Na")
    }
    else{
      updateSelectInput(session, inputId = "col.Na1", selected = "None")
    }
    #
    pos <- current.names == "K"
    if(sum(pos) == 1){
      updateSelectInput(session, inputId = "col.K1", selected = "K")
    }
    else{
      updateSelectInput(session, inputId = "col.K1", selected = "None")
    }
    #
    pos <- current.names == "HCO3"
    if(sum(pos) == 1){
      updateSelectInput(session, inputId = "col.HCO31", selected = "HCO3")
    }
    else{
      updateSelectInput(session, inputId = "col.HCO31", selected = "None")
    }
    #
    pos <- current.names == "CO3"
    if(sum(pos) == 1){
      updateSelectInput(session, inputId = "col.CO31", selected = "CO3")
    }
    else{
      updateSelectInput(session, inputId = "col.CO31", selected = "None")
    }
    #
    pos <- current.names == "Cl"
    if(sum(pos) == 1){
      updateSelectInput(session, inputId = "col.Cl1", selected = "Cl")
    }
    else{
      updateSelectInput(session, inputId = "col.Cl1", selected = "None")
    }
    #
    pos <- current.names == "SO4"
    if(sum(pos) == 1){
      updateSelectInput(session, inputId = "col.SO41", selected = "SO4")
    }
    else{
      updateSelectInput(session, inputId = "col.SO41", selected = "None")
    }
  }
  #
  create_gdata <- function(){
    the.table <- server.env$current.table
    #print("CO3")
    #print(input$col.CO31)
    if(input$col.CO31 == "None"){
      the.table$CO3 <- vector("numeric", length = nrow(the.table))
      chem.cols <- c(input$col.Ca1, input$col.Mg1, input$col.Na1, input$col.K1,
                     input$col.HCO31, "CO3", input$col.Cl1, input$col.SO41)
    }
    else{
      chem.cols <- c(input$col.Ca1, input$col.Mg1, input$col.Na1, input$col.K1,
                     input$col.HCO31, input$col.CO31, input$col.Cl1, input$col.SO41)
    }
    #print("chem")
    #print(chem.cols)
    phys.cols <- c(input$col.pH1, input$col.TDS1, input$col.EC1)
    all.cols <- c(chem.cols, phys.cols)
    #print("all")
    #print(all.cols)
    pos <- all.cols != "None"
    #print("POS")
    #print(pos)
    #print("names")
    #print(the.table[pos])
    all.cols <- all.cols[pos]
    #
    #print("current")
    #print( server.env$current.names)
    pos <- all.cols %in% names(the.table)
    #print("POS")
    #print(pos)
    the.table <- the.table[pos]
    #print("Names")
    #print(names(the.table))
    input$create.gdata
    server.env$current.gdata <- isolate(
      geochemical_dataset(name = "GeochemicalDataset", data = the.table)
    )
    server.env$current.geonames <- names(the.table)
  }
  #
  observeEvent(input$create.gdata, {
    create_gdata()
    shinyalert(title = "Geochemical Dataset Defined!!!", type = "success")
  })
  #
  observeEvent(input$col.fit,{
    fit_columns()
  })
  #######################################################################################
  #                               Transformation tab
  #######################################################################################

  #######################################################################################
  #                               EDA tab
  #######################################################################################
  output$eda.varselector <- renderUI({
    the.table <- server.env$current.table
    if(is.null(the.table)){
      return(NULL)
    }
    varnames <- c("None", server.env$current.geonames)
    res <- selectInput(inputId = "eda.varselector1", label = "Variable",
                       choices = varnames, selected = "None")
    return(res)
  })
  #
  output$eda.plot <- renderPlot({
    current.gdata <- server.env$current.gdata
    validate(
      need(!is.null(current.gdata), "The Geochemical Dataset is not defined")
    )
    if(is.null(current.gdata))
      return(NULL)
    #
    current.varname <- input$eda.varselector1
    if(is.null(current.varname) || current.varname == "None")
      return(NULL)
    #
    #print(current.varname)
    #print(current.gdata$dataset)
    width <- vector("numeric", length = nrow(current.gdata$dataset[current.varname]))
    width[1:length(width)] <- 1
    df <- data.frame(x = current.gdata$dataset[current.varname],
                     y =width)
    #print(df)
    p1 <- ggplot(df, aes_string(x = current.varname)) +
      geom_histogram() +
      theme_bw() +
      ggtitle("a) Histogram")
    if(input$eda.log == "Yes"){
      p1 <- p1 + scale_x_log10()
    }
    p2 <- ggplot(df, aes_string(x = "width", y = current.varname)) +
      geom_boxplot() +
      theme_bw() +
      ggtitle("a) Boxplot")
    if(input$eda.log  == "Yes"){
      p2 <- p2 + scale_y_log10()
    }
    p3 <- ggplot(df, aes_string(x = current.varname)) +
      stat_ecdf(geom = "step") +
      theme_bw() +
      ggtitle("c) ECDF")
    if(input$eda.log  == "Yes"){
      p3 <- p3 + scale_x_log10()
    }
    p4 <- ggplot(df, aes_string(sample = current.varname)) +
      stat_qq() +
      stat_qq_line() +
      theme_bw() +
      ggtitle("d) QQ Plot")
    if(input$eda.log == "Yes"){
      p4 <- p4 + scale_y_log10()
    }
    pdef <- grid.arrange(p1, p2, p3, p4,
                 ncol = 2)
    server.env$current.plot <- pdef
    return(pdef)
  })
  #
  output$eda.downloadPlot <- downloadHandler(
      filename <- function() {
        paste('plot1', 'png', sep = ".")
      },
      content <- function(file) {
        png(file)
        plot <- server.env$current.plot
        print(plot)
        dev.off()
      },
      contentType = "image/png"
    )

  #######################################################################################
  #                               Hydrogeochemical Plots tab
  #######################################################################################
  output$hplot.tselector <- renderUI({
    res <- NULL
    current.table <- server.env$current.table
    if(is.null(current.table)){
      return(NULL)
    }
    else{
      res <- selectInput(inputId = "hplot.tselector0", label = "Dataset type",
                         choices = c("None", "Single Sample", "All Samples"),
                         selected = "None")
    }
    return(res)
  })
  #
  output$hplot.tselector1 <- renderUI({
    res <- NULL
    current.table <- server.env$current.table
    if(is.null(current.table)){
      return(NULL)
    }
    current.hplot <- input$hplot.tselector0
    if(is.null(current.hplot))
      return(NULL)
    if(input$hplot.tselector0 == "Single Sample"){
      res <- selectInput(inputId = "hplot.tselector1a", label = "Plot type",
                         choices = plot.types.single, selected = "None")
    }
    else{
      res <- selectInput(inputId = "hplot.tselector1a", label = "Plot type",
                         choices = plot.types.all, selected = "None")
    }
    return(res)
  })
  #
  output$hplot.option1 <- renderUI({
    res <- NULL
    current.table <- server.env$current.table
    if(is.null(current.table)){
      return(NULL)
    }
    #
    current.hplot <- input$hplot.tselector0
    if(is.null(current.hplot))
      return(NULL)
    #
    current.tplot <- input$hplot.tselector1a
    if(is.null(current.tplot) || current.tplot == "None")
      return(NULL)
    #
    if(current.hplot == "Single Sample"){
      # all.variables <- server.env$current.names
      # res <- selectInput(inputId = "hplot.variables",
      #                    label = "Choose variables",
      #                    choice = all.variables,
      #                    selected = "Ca",
      #                    multiple = TRUE)
      all.samples <- c("None", as.character( seq(1, nrow(current.table), by = 1)))
      res <- selectInput(inputId = "hplot.samples",
                         label = "Choose Sample",
                         choices = all.samples,
                         selected = "None",
                         multiple = FALSE)
    }
    else if(current.hplot == "All Samples"){
      if(current.tplot == "ternary"){
        all.variables <- server.env$current.names
        res <- selectInput(inputId = "hplot.variables",
                           label = "Choose variables",
                           choice = all.variables,
                           selected = "Ca",
                           multiple = TRUE)
      }
    }
    #
    return(res)
  })
  #
  output$hplot.option2 <- renderUI({
    res <- NULL
    current.table <- server.env$current.table
    if(is.null(current.table)){
      return(NULL)
    }
    #
    current.tplot <- input$hplot.tselector1a
    if(is.null(current.tplot) || current.tplot == "None")
      return(NULL)
    #
    all.variables <- c("None", server.env$current.names)
    res <- selectInput(inputId = "hplot.color", label = "Choose variable to color",
                       choices = all.variables, selected = "None", multiple = FALSE)
    return(res)
  })
  #
  output$hplot.option3 <- renderUI({
    res <- NULL
    current.table <- server.env$current.table
    if(is.null(current.table)){
      return(NULL)
    }
    current.tplot <- input$hplot.tselector1a
    if(is.null(current.tplot) || current.tplot == "None")
      return(NULL)
    #
    res <- selectInput(inputId = "hplot.measure", label = "Measure to use",
                       choices = measures.type, selected = "None", multiple = FALSE)
    return(res)
  })
  #
  output$hplot.option4 <- renderUI({
    res <- NULL
    current.table <- server.env$current.table
    if(is.null(current.table)){
      return(NULL)
    }
    #
    current.tplot <- input$hplot.tselector1a
    if(is.null(current.tplot) || current.tplot == "None")
      return(NULL)
    #
    all.variables <- c("None", server.env$current.names)
    res <- selectInput(inputId = "hplot.size", label = "Choose variable to Size",
                       choices = all.variables, selected = "None", multiple = FALSE)
    return(res)
  })
  #
  output$hplot <- renderPlot({
    current.table <- server.env$current.table
    if(is.null(current.table)){
      return(NULL)
    }
    #
    current.hplot <- input$hplot.tselector0
    if(is.null(current.hplot) || current.hplot == "None")
      return(NULL)
    #
    current.tplot <- input$hplot.tselector1a
    if(is.null(current.tplot) || current.tplot == "None")
      return(NULL)
    #
    current.vars <- input$hplot.variables
    if(current.tplot == "ternary"){
      if(is.null(current.vars) || current.vars == "None")
        return(NULL)
    }
    #
    #
    #
    current.color <- input$hplot.color
    if(is.null(current.color))
      return(NULL)
    if(current.color == "None")
      current.color <- NULL
    #
    current.size <- input$hplot.size
    if(is.null(current.size))
      return(NULL)
    if(current.size == "None")
      current.size <- NULL
    #
    current.measure <- input$hplot.measure
    if(is.null(current.measure) || current.measure == "None")
      return(NULL)
    p1 <- NULL
    current.gdata <- server.env$current.gdata
    if(current.hplot == "All Samples"){
      if(current.tplot == "piper"){
        p1 <- plot(current.gdata, type = "piper",
                   measure = current.measure,
                   color = current.color,
                   Size = current.size)
      }
      else if(current.tplot == "modified_piper"){
        p1 <- plot(current.gdata, type = "modified_piper")
      }
      else if(current.tplot == "durov"){
        p1 <- plot(current.gdata, type = "durov",
                   measure = current.measure,
                   vars = current.vars,
                   color = current.color,
                   Size = current.size)
      }
      else if(current.tplot == "schoeller"){
        p1 <- plot(current.gdata, type = "schoeller",
                             measure = current.measure,
                             color = current.color)
      }
      else if(current.tplot == "multirectangular"){
        p1 <- plot(current.gdata, type = "multirectangular",
                   measure = current.measure,
                   vars = current.vars,
                   color = current.color,
                   Size = current.size)
      }
      else if(current.tplot == "ternary"){
        p1 <- plot(current.gdata, type = "ternary",
                   measure = current.measure,
                   vars = current.vars,
                   color = current.color,
                   Size = current.size)
      }
    }
    else if(current.hplot == "Single Sample"){
      current.sample <- as.integer(input$hplot.samples)
      if(current.tplot == "stiff"){
        p1 <- plot(current.gdata[current.sample,], type = "stiff",
                   measure = current.measure)
      }
      else if(current.tplot == "radial"){
        p1 <- plot(current.gdata[current.sample], type = "radial",
                   measure = current.measure)
      }
    }
    print(p1)
    return(p1)
  })

})
