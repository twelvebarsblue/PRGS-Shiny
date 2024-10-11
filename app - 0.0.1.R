library(shiny)
library(mrgsolve)
library(tidyverse)
library(gridExtra)
library(shinydashboard)
library(RCurl)
library(plotly)



####################################
# Essential functions              #
####################################
mdrd <- 
  function(serum, age, gender) {
    if (gender == "Male") {
      return(186 * ((serum/88.4)^-1.154) * (age^-0.203))}
    if (gender == "Female") {
      return(186 * ((serum/88.4)^-1.154) * (age^-0.203) * 0.742)}
  }

cockcroft <-
  function(serum, age, weight, gender) {
    if (gender == "Male") {
      return((1.23 * (140-age) * weight) / serum)}
    if (gender == "Female") {
      return((1.04 * (140-age) * weight) / serum)}
  }

####################################
# User interface                   #
####################################


ui <- dashboardPage(skin = "purple",
                    dashboardHeader(title = "Vancomycin Dosing Guide", titleWidth = 126),
                    dashboardSidebar(disable = TRUE),
                    dashboardBody(tags$head(
                      tags$style(HTML('.value-box {font-size: 11px;} 
                                      .info-box {min-height: 45px; margin-top: 10px; padding: 0px; font-size: 11px;} 
                                      .info-box-icon {height: 45px; line-height: 0px; padding-top: 12px; margin-right: 6px; width: 20%} 
                                      .info-box-content {margin: 0px; padding-top: 0px; padding-bottom: 0px; text-align: left;font-size: 8px;}'))),
                      column(6, style = "padding:0px",
                             box(width = 12, #height = "610px",
                                 tabBox(id = "Steps", width = 12,
                                        tabPanel("Model",
                                                 fluidRow(
                                                   column(3, 
                                                          selectInput(
                                                            "model", 
                                                            "PkPop Model", 
                                                            choices = list(
                                                              "Model 1" = 1, 
                                                              "Model 2" = 2, 
                                                              "Model 3" = 3, 
                                                              "Model 4" = 4,
                                                              "Model 5" = 5),
                                                            selected = 1)),
                                                   column(width = 9,
                                                          column(3, selectInput("gender", "Gender",choices = list("Male" = "Male","Female" = "Female"),selected = "Female")),
                                                          column(3, numericInput("age", "Age",50, min = 17,max = 100)),
                                                          column(3, selectInput("predCrcl", "CRCL model", choices = list("Literature" = 1, "MDRD" = 2, "Cockcroft-Gault" = 3), selected = 1)),
                                                          conditionalPanel(condition = "input.predCrcl == 2 | input.predCrcl == 3", column(3, numericInput("creatinine", "Cr (mg/dL)", 50, min = 0.1))),
                                                          conditionalPanel(condition = "input.predCrcl == 3", column(3, numericInput("weight", "TBW (Kg)",70,min = 20,max = 150))),
                                                          conditionalPanel(condition = "input.model == 2", column(3, numericInput("weight", "TBW (Kg)",70,min = 20,max = 150)))))),
                                        tabPanel("Dosage Screening",
                                                 fluidRow(
                                                   column(4, numericInput("min_dose", "Minimum Dose (mg)", value = 100, min = 1, step=1)),
                                                   column(4, numericInput("max_dose", "Maximal Dose (mg)", value = 2000, min = 1, step=1)),
                                                   column(4, numericInput("resolution", "Resolution", value = 100, min = 1, step=1)),
                                                   column(4, numericInput("ii", "Interdose interval (h)", value = 12, min = 1, max = 48, step=1)),
                                                   column(12, verbatimTextOutput("output_text")),
                                                   ))),
                                 tabBox(id = "PkPop", width = 12,
                                        tabPanel("Pharmacokinetics",
                                                 fluidRow(style = "padding: 0px;", 
                                                          valueBoxOutput("TFG", width = 3),
                                                          conditionalPanel(
                                                            condition = "input.model == 1",
                                                            valueBoxOutput("Cl", width = 2),
                                                            valueBoxOutput("Vc", width = 2),
                                                            valueBoxOutput("Vp", width = 2),
                                                            valueBoxOutput("Q", width = 2)),
                                                          conditionalPanel(
                                                            condition = "output.Revilla",
                                                            valueBoxOutput("ClR", width = 2),
                                                            valueBoxOutput("V", width = 2)),
                                                          conditionalPanel(
                                                            condition = "output.Goti",
                                                            valueBoxOutput("ClG", width = 2),
                                                            valueBoxOutput("VcG", width = 2),
                                                            column(2, numericInput("V2", "Vp", value = 38.4, min = 0, max = 10, step=0.01)),
                                                            column(2, numericInput("Q", "Q", value = 6.5, min = 0, max = 20, step=0.01))))),
                                        tabPanel("Coefficients",
                                                 fluidRow(
                                                   conditionalPanel(
                                                     condition = "input.model == 1",
                                                     column(2, numericInput("teta1", "Θ_1", value = 0.034, min = 0, max = 10, step=0.01)),
                                                     column(2, numericInput("teta2", "Θ_2", value = 0.015, min = 0, max = 10, step=0.01)),
                                                     column(2, numericInput("teta3", "Θ_3", value = 0.414, min = 0, max = 10, step=0.01)),
                                                     column(2, numericInput("teta4", "Θ_4", value = 7.48, min = 0, max = 10, step=0.01)),
                                                     column(2, numericInput("teta5", "Θ_5", value = 1.32, min = 0, max = 20, step=0.01))),
                                                   conditionalPanel(
                                                     condition = "output.Revilla",
                                                     column(2, numericInput("teta11", "Θ_1", value = 0.67, min = 0, max = 10, step=0.01)),
                                                     column(2, numericInput("teta22", "Θ_2", value = -0.24, min = -1, max = 10, step=0.01)),
                                                     column(2, numericInput("teta33", "Θ_3", value = 0.82, min = 0, max = 10, step=0.01)),
                                                     column(2, numericInput("teta44", "Θ_4", value = 2.49, min = 0, max = 10, step=0.01))),
                                                   conditionalPanel(
                                                     condition = "output.Goti",
                                                     column(2, numericInput("teta1G", "Θ_1", value = 4.5, min = 0, max = 10, step=0.01)),
                                                     column(2, numericInput("teta2G", "Θ_2", value = 0.8, min = -1, max = 10, step=0.01)),
                                                     column(2, numericInput("teta3G", "Θ_3", value = 0.7, min = 0, max = 10, step=0.01)),
                                                     column(2, numericInput("teta4G", "Θ_4", value = 58.4, min = 0, max = 10, step=0.01)),
                                                     column(2, numericInput("teta5G", "Θ_5", value = 0.5, min = 0, max = 20, step=0.01))))),
                                        tabPanel("Equations",
                                                 fluidRow(
                                                   conditionalPanel(
                                                     condition = "input.model == 1",
                                                     withMathJax(),
                                                     box(width = 4, helpText('$$Cl = \\theta_1 \\cdot Cl_{cr} + \\theta_2 \\cdot TBW $$', style = "font-size:10px",)),
                                                     box(width = 3, helpText('$$V_c = \\theta_3 \\cdot TBW $$', style = "font-size:10px",)),
                                                     box(width = 3, helpText('$$V_p = \\theta_5 \\cdot TBW $$', style = "font-size:10px",)),
                                                     box(width = 2, helpText('$$Q = \\theta_4 $$', style = "font-size:10px"))),
                                                   conditionalPanel(
                                                     condition = "output.Revilla",
                                                     withMathJax(),
                                                     box(width = 4, helpText('$$Cl = \\theta_1 \\cdot Cl_{cr} +  TBW^{\\theta_2} $$', style = "font-size:10px")),
                                                     box(width = 3, helpText('$$V = \\theta_3 \\cdot \\theta_4^{A} $$', style = "font-size:10px")),
                                                     box(width = 5, helpText('$$A = 0 \\text{ if } Cr < 1 \\text{; } A = 1 \\text{ if } Cr > 1$$', style = "font-size:10px"))),
                                                   conditionalPanel(
                                                     condition = "output.Goti",
                                                     withMathJax(),
                                                     box(width = 4, helpText('$$Cl = \\theta_1 \\cdot \\left(\\frac{Cl_{cr}}{120}\\right)^{\\theta_2} \\cdot \\theta_3^{DIAL} $$', style = "font-size:10px")),
                                                     box(width = 4, helpText('$$V_c = \\theta_4 \\cdot \\left(\\frac{TBW}{70}\\right) \\cdot \\theta_5^{DIAL} $$', style = "font-size:10px")),
                                                     box(width = 4, helpText('Clcr = 150 mL/min if Clcr > 150 mL/min. SCr = 1 mg/dL if SCr < 1 mg/dL and age > 65 years.', style = "font-size:10px")))))),
                                 tabBox(id = "Configuration", width = 12,
                                        tabPanel("Tolerance",
                                                 fluidRow(
                                                   column(2, numericInput("Llim", "Lower (mg/L)", value = 15, min = 1, step=1)),
                                                   column(2, numericInput("Ulim", "Upper (mg/L)", value = 50, min = 1, step=1))
                                                   )),
                                        tabPanel("Settings",
                                                 fluidRow(
                                                   column(3,numericInput("nsamples", "Number of samples for simulation", value = 1000, min=0)),
                                                   column(3,numericInput("simtime", "Simulation time", value = 36, min=0)),
                                                   column(3,numericInput("minprob", "Minimal probability", value = 0.25, min=0, max=1)),
                                                   column(3,numericInput("maxprob", label='Maximal probability', value = 0.75,min=0,max=1))
                                                   ))),
                                 column(12, style="margin-top: 2px; text-align:right", 
                                        column(10, style="margin-top: 10px", conditionalPanel(condition = "input.go > 0", downloadButton(style = "color:black;background-color:lightgray;font-size:12px;", "downloadData", "Download"))),
                                        column(2, style="margin-top: 10px", actionButton(style = "color:white; background-color:green;font-size:12px;", "go", strong(" Simulate"), icon = icon("play")))))),
                      
   
                      column(6, style = "padding:0px",
                             box(width = 12, #height = "610px",
                                 tabBox(id = "Output", width = 12,
                                        tabPanel("Plot",
                                                 fluidRow(
                                                   column(3, plotlyOutput(width = "600px", height = "400px", outputId="simPlot")))),
                                        tabPanel("Prediction",
                                                 fluidRow(
                                                   column(9, 
                                                          infoBoxOutput("clearance", width = 6),
                                                          infoBoxOutput("initDose", width = 6),
                                                          infoBoxOutput("dose2", width = 6),
                                                          infoBoxOutput("dose3", width = 6),
                                                          infoBoxOutput("dose4", width = 6),
                                                          ))))))
                      ))
                      
                      
                      
                                 


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  txt <- reactive({
    
    runSim <- function(min_dose, max_dose, resolution) {
      
      # Create combinations of dose
      studied_doses <- seq(min_dose, max_dose, resolution)
      
      combination <- as.data.frame(expand.grid(DOSE_1 = studied_doses, DOSE_2 = studied_doses, DOSE_3 = studied_doses, DOSE_4 = studied_doses))
      
      nrow(combination)
      
    }
    
    paste0("The total number of combinations is ", runSim(input$min_dose, input$max_dose, input$resolution))
    
    }) 
  
  output$output_text <- renderText({
      txt()
    })
  
  doseSim <- eventReactive(input$go, {
    
    #########################################
    # Download models
    linkMod1 <- "https://github.com/twelvebarsblue/PRGS-Shiny/raw/refs/heads/main/run001.cpp?raw=true"
    linkCRCL_reference <- "https://raw.githubusercontent.com/twelvebarsblue/PRGS-Shiny/refs/heads/main/CRCL_reference.csv"
    
    download.file(linkMod1, destfile = "run001.cpp")
    download.file(linkCRCL_reference, destfile = "CRCL_reference.csv", method = "curl")
    
    #########################################
    # Load models
    
    if (input$model == 1) {
      mod <- mread("./run001.cpp") }
    else if (input$model == 2) {
      mod <- mread("./run002.cpp") }
    else if (input$model == 3) {
      mod <- mread("./run003.cpp") }
    else if (input$model == 4) {
      mod <- mread("./run004.cpp") }
    else if (input$model == 5) {
      mod <- mread("./run005.cpp") }
    
    
    ################################################
    # Step 1: Obtaining creatinine clearance values
    
    # Option 1:
    
    # Create a model based on CrCL reference values
    # Male: 
    # Rowe, John W., et al. "The effect of age on creatinine clearance in men: a cross-sectional and longitudinal study." 
    # Journal of gerontology 31.2 (1976): 155-163.
    
    # Female
    # Sokoll, Lori J., et al. "Establishment of creatinine clearance reference values for older women." 
    # Clinical chemistry 40.12 (1994): 2276-2281.
    
    # Option 2/3:
    # MDRD/Cockcroft-Gault
    
    if (input$predCrcl == 1) {
      
      refRange <- read.csv("./CRCL_reference.csv")
      
      # LR model
      lm_crcl <- lm(CrCL ~ Gender + Age, refRange)
      new_data <- data_frame(Gender = input$gender, Age = input$age)
      
      crcl <- predict(lm_crcl, newdata = new_data)}
    
    else if (input$predCrcl == 2) {
      crcl <- mdrd(input$creatinine, input$age, input$gender)}
    
    else if (input$predCrcl == 3) {
      crcl <- cockcroft(input$creatinine, input$age, input$weight, input$gender)}
      
    #####################################################
    # Run simulation
    
    runSim <- function(nsamples, minprob, maxprob, crcl, age, sim_time, ii, min_dose, max_dose, resolution, llim, ulim) {
      
      # Create combinations of dose
      studied_doses <- seq(min_dose, max_dose, resolution)
      
      combination <- as.data.frame(expand.grid(CRCL = crcl, AGE = age, 
                                               DOSE_1 = studied_doses, DOSE_2 = studied_doses, DOSE_3 = studied_doses, DOSE_4 = studied_doses, 
                                               UNTIL = sim_time, TIME_INT = ii, INTERVAL = ii))
      
      combination$TAG = seq(nrow(combination))
    
      TIME_1 = 0 + ii
      TIME_2 <- TIME_1 + ii
      TIME_3 <- TIME_2 + ii
      
      # At each combination, simulate for x times
      withProgress(message = "Simulating ", value = 0, {
        
        n = nrow(combination)
        
        for (i in 1:n) {
          
          df <- filter(combination, TAG == i)
          
          repeated_df <- df[rep(rownames(df), each = nsamples), ]
          
          repeated_df$ID <- seq(1:nsamples)
          
          out <- mrgsim(mod, idata = repeated_df, end = sim_time)
          
          out <- as.data.frame(out)
          
          sum_stat <- out %>%
            group_by(time) %>%
            summarise(Median_C=median(DV),
                      Median_AUC=median(AUC),
                      LP_DV=quantile(DV,probs=minprob),
                      HP_DV=quantile(DV,probs=maxprob),
                      LP_AUC=quantile(AUC,probs=minprob),
                      HP_AUC=quantile(DV,probs=maxprob)
            )  
          
          
          # Check the concentration at 12 hours, if < 15 mg/mL, go next combination
          
          dv = filter(sum_stat, time %in% c(TIME_1, TIME_2, TIME_3))
          
          incProgress(1/n, detail = paste("Combination ", i)) # To check whether its running or not
          
          if (all(dv$Median_C > llim) & max(sum_stat$Median_C) < ulim) {
            # cat(paste0("Recommended loading dose: ", df$DOSE_1, " mg", "\n", 
            #            "Maintenance dose at 12 h: ", df$DOSE_2, " mg", "\n", 
            #            "Maintenance dose at 24 h: ", df$DOSE_3, " mg", "\n",
            #            "Maintenance dose at 36 h: ", df$DOSE_4, " mg", "\n"))
          
            
            plot_pk <- ggplot(sum_stat, aes(x=time,y=Median_C)) +
              
              ## Add ribbon for variability
              geom_ribbon(aes(ymin=LP_DV, ymax=HP_DV, x=time), alpha = 0.15, linetype=0)+
              
              geom_ribbon(aes(ymin=llim,ymax=ulim), fill="green", alpha=0.3) +
              
              ## Add median line
              geom_line(size=2) +
              
              # scale_y_log10()+
              
              # Set axis and theme
              ylab(paste("Concentration (mg/L)",sep=""))+
              xlab("Time after dose (h)")+
              theme_bw()+
              
              # Remove legend
              theme(legend.position="none")
            
            # print(plot_pk)
            
            crcl <<- crcl
            TIME_1 <<- TIME_1
            TIME_2 <<- TIME_2
            TIME_3 <<- TIME_3
            initDose <<- df$DOSE_1
            dose2 <<- df$DOSE_2
            dose3 <<- df$DOSE_3
            dose4 <<- df$DOSE_4
            plot_pk <<- plot_pk
            
            break
          }
        }
      })
    }
    
    runSim(nsamples = input$nsamples, minprob = input$minprob, maxprob = input$maxprob, 
           crcl = crcl, age = input$age, sim_time = input$simtime, ii = input$ii, 
           min_dose = input$min_dose, max_dose = input$max_dose, resolution = input$resolution,
           llim = input$Llim, ulim = input$Ulim)
    
    return(list(crcl = crcl, TIME_1 = TIME_1, TIME_2 = TIME_2, TIME_3 = TIME_3, initDose = initDose, dose2 = dose2, dose3 = dose3, dose4 = dose4, plot_pk = plot_pk))
  })
  
  
  output$downloadData <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "report.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      params <- list(n = input$slider)
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv()))})

                        

  observeEvent(doseSim(),{
    output$clearance <- renderInfoBox({
      infoBox(
        "Clearance", paste0(round(doseSim()$crcl), " mL/min"), icon = tags$i(class = "fas fa-hospital-user", style="font-size: 12px"),
        color = "purple"
      )
    })
    output$initDose <- renderInfoBox({
      infoBox(
        "Initial Dose", paste0(round(doseSim()$initDose), " mg/mL"), icon = tags$i(class = "fas fa-prescription", style="font-size: 12px"),
        color = "purple"
      )
    })
    output$dose2 <- renderInfoBox({
      infoBox(
        paste0("Dose at ", doseSim()$TIME_1, "h"), paste0(round(doseSim()$dose2), " mg/mL"), icon = tags$i(class = "fas fa-prescription", style="font-size: 12px"),
        color = "purple"
      )
    })
    output$dose3 <- renderInfoBox({
      infoBox(
        paste0("Dose at ", doseSim()$TIME_2, "h"), paste0(round(doseSim()$dose3), " mg/mL"), icon = tags$i(class = "fas fa-prescription", style="font-size: 12px"),
        color = "purple"
      )
    })
    output$dose4 <- renderInfoBox({
      infoBox(
        paste0("Dose at ", doseSim()$TIME_3, "h"), paste0(round(doseSim()$dose4), " mg/mL"), icon = tags$i(class = "fas fa-prescription", style="font-size: 12px"),
        color = "purple"
      )
    })
  })
   

  
  output$simPlot <- renderPlotly({
    withProgress(message = 'Creating plot', detail='Wait...', value = 0.1, expr={
      ggplotly(doseSim()$plot_pk) %>%
        config(
          displaylogo = FALSE,
          modeBarButtonsToRemove = c("lasso2d","zoomIn2d", "zoomOut2d","pan2d","autoScale2d","zoom2d","hoverClosestCartesian"))
    
    })
  }) 
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
