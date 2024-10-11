library(shiny)
library(mrgsolve)
library(tidyverse)
library(gridExtra)
library(shinydashboard)
library(plotly)
library(mapbayr)
library(PKNCA)

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
# Model description                #
####################################

desc_m1 <- "Pharmacokinetic/pharmacodynamic (PK/PD)-modeling links dose-concentration relationships (PK) and concentration-effect relationships (PD), thereby facilitating the description and prediction of the time course of drug effects resulting from a certain dosing regimen. PK/PD-modeling approaches can basically be distinguished by four major attributes."


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
                                 tabBox(id = "Model", width = 12,
                                        tabPanel("Model",
                                                 fluidRow(
                                                   column(6, 
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
                                                   column(6, 
                                                          selectInput(
                                                            "method", 
                                                            "Simulation Method", 
                                                            choices = list(
                                                              "Monte Carlo" = 1, 
                                                              "MAP Bayesian" = 2, 
                                                              "Manual" = 3),
                                                            selected = 1)))), 
                                        
                                        tabPanel("Description",
                                                 fluidRow(
                                                   conditionalPanel(
                                                     condition = "input.model == 1",
                                                     withMathJax(),
                                                     box(width = 12, helpText(desc_m1, style = "font-size:10px",))),
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
                                                     box(width = 4, helpText('Clcr = 150 mL/min if Clcr > 150 mL/min. SCr = 1 mg/dL if SCr < 1 mg/dL and age > 65 years.', style = "font-size:10px"))))),
                                        
                                        
                                        tabPanel("Equations",
                                                 fluidRow(
                                                   conditionalPanel(
                                                     condition = "input.model == 1",
                                                     withMathJax(),
                                                     box(width = 4, helpText('$$Cl = \\theta_1 \\cdot Cl_{cr} + \\theta_2 \\cdot TBW $$', style = "font-size:10px",)),
                                                     box(width = 3, helpText('$$V_c = \\theta_3 \\cdot TBW $$', style = "font-size:10px",)),
                                                     box(width = 3, helpText('$$V_p = \\theta_5 \\cdot TBW $$', style = "font-size:10px",)),
                                                     box(width = 2, helpText('$$Q = \\theta_4 $$', style = "font-size:10px")))
                                                   
                                                   
                                                   ))
                                        ),
                                        
                                 tabBox(id = "Characteristics", width = 12,
                                        tabPanel("Patient's Characteristics",
                                                 fluidRow(style = "padding: 0px;",
                                                          column(width = 9,
                                                                 column(4, selectInput("gender", "Gender",choices = list("Male" = "Male","Female" = "Female"),selected = "Female")),
                                                                 column(4, numericInput("age", "Age",50, min = 17,max = 100)),
                                                                 conditionalPanel(condition = "input.model == 2", column(4, numericInput("weight", "TBW (Kg)",70,min = 20,max = 150))),
                                                                 column(4, selectInput("predCrcl", "CRCL model", choices = list("Literature" = 1, "MDRD" = 2, "Cockcroft-Gault" = 3), selected = 1)),
                                                                 conditionalPanel(condition = "input.predCrcl == 2 | input.predCrcl == 3", column(4, numericInput("creatinine", "Cr (mg/dL)", 50, min = 0.1))),
                                                                 conditionalPanel(condition = "input.predCrcl == 3", column(4, numericInput("weight", "TBW (Kg)",70,min = 20,max = 150))))
                                                          )),
                                        tabPanel("Dose",
                                                 fluidRow(
                                                   conditionalPanel(
                                                     condition = "input.method == 1",
                                                     column(4, numericInput("min_dose", "Minimum Dose (mg)", value = 100, min = 1, step=1)),
                                                     column(4, numericInput("max_dose", "Maximal Dose (mg)", value = 2000, min = 1, step=1)),
                                                     column(4, numericInput("resolution", "Resolution", value = 100, min = 1, step=1)),
                                                     column(4, numericInput("ii", "Interdose interval (h)", value = 12, min = 1, max = 48, step=1)),
                                                     column(12, verbatimTextOutput("output_text"))),
                                                   
                                                   conditionalPanel(
                                                     condition = "input.method == 2",
                                                     column(12, numericInput("initDose", "Loading Dose (mg)", value = 1000, min = 1, step=1)),
                                                     column(6, 
                                                            numericInput("time1", "Time 1 (h)", value = 50.25),
                                                            numericInput("time2", "Time 2 (h)", value = 85.95),
                                                            numericInput("time3", "Time 3 (h)", value = 120.1)),
                                                     column(6, 
                                                            numericInput("dv1", "Conc 1 (mg/mL)", value = 13.58),
                                                            numericInput("dv2", "Conc 2 (mg/mL)", value = 40.18),
                                                            numericInput("dv3", "Conc 3 (mg/mL)", value = 20.29))),
                                                   
                                                   conditionalPanel(
                                                     condition = "input.method == 3",
                                                     column(12, numericInput("init_dose", "Loading Dose (mg)", value = 1000, min = 1, step=1)),
                                                     column(6, 
                                                            numericInput("time_1", "Time 1 (h)", value = 12, min = 1, step=1),
                                                            numericInput("time_2", "Time 2 (h)", value = 24, min = 1, step=1),
                                                            numericInput("time_3", "Time 3 (h)", value = 36, min = 1, step=1),
                                                            numericInput("time_4", "Time 4 (h)", value = 48, min = 1, step=1)),
                                                     column(6, 
                                                            numericInput("d1", "Dose 1 (mg/mL)", value = 1000, min = 1, step=1),
                                                            numericInput("d2", "Dose 2 (mg/mL)", value = 1000, min = 1, step=1),
                                                            numericInput("d3", "Dose 3 (mg/mL)", value = 1000, min = 1, step=1),
                                                            numericInput("d4", "Dose 4 (mg/mL)", value = 1000, min = 1, step=1)))
        
                                                   
                                                   ))),
                                        
                                 
                                 
                                 tabBox(id = "Configuration", width = 12,
                                        tabPanel("Therapeutic Index",
                                                 fluidRow(
                                                   column(2, numericInput("Llim", "Minimum Effective Concentration (mg/L)", value = 15, min = 1, step=1)),
                                                   column(2, numericInput("Ulim", "Minimum Toxic Concentration (mg/L)", value = 50, min = 1, step=1))
                                                   )),
                                        tabPanel("Settings",
                                                 fluidRow(
                                                   conditionalPanel(condition = "input.method == 1 | input.method == 3", column(3,numericInput("nsamples", "Number of samples for simulation", value = 1000, min=0))),
                                                   column(3,numericInput("simtime", "Simulation time", value = 60, min=0)),
                                                   conditionalPanel(condition = "input.method == 1 | input.method == 3", column(3,numericInput("minprob", "Minimal probability", value = 0.25, min=0, max=1))),
                                                   conditionalPanel(condition = "input.method == 1 | input.method == 3", column(3,numericInput("maxprob", label='Maximal probability', value = 0.75,min=0,max=1))),
                                                   conditionalPanel(
                                                     condition = "input.method == 2",
                                                     column(4, numericInput("min_dose_bayes", "Minimum Dose (mg)", value = 100, min = 1, step=1)),
                                                     column(4, numericInput("max_dose_bayes", "Maximal Dose (mg)", value = 2000, min = 1, step=1)),
                                                     column(4, numericInput("resolution_bayes", "Resolution", value = 100, min = 1, step=1)),
                                                     column(4, numericInput("ii_bayes", "Interdose interval (h)", value = 12, min = 1, max = 48, step=1)))
                                                   ))),
                                 column(12, style="margin-top: 2px; text-align:right", 
                                        column(10, style="position:absolute;top:-18px;right:7em", conditionalPanel(condition = "input.go1", downloadButton(style = "color:black;background-color:lightgray;font-size:12px;", "downloadData1", "Download"))),
                                        column(10, style="position:absolute;top:-18px;right:7em", conditionalPanel(condition = "input.go2", downloadButton(style = "color:black;background-color:lightgray;font-size:12px;", "downloadData2", "Download"))),
                                        column(10, style="position:absolute;top:-18px;right:7em", conditionalPanel(condition = "input.go3", downloadButton(style = "color:black;background-color:lightgray;font-size:12px;", "downloadData3", "Download"))),
                                        column(2, style="position:absolute;top:-18px;right:1em", conditionalPanel(condition = "input.method == 1", actionButton(style = "color:white; background-color:purple;font-size:12px;", "go1", strong(" Simulate"), icon = icon("play")))),
                                        column(2, style="position:absolute;top:-18px;right:1em", conditionalPanel(condition = "input.method == 2", actionButton(style = "color:white; background-color:purple;font-size:12px;", "go2", strong(" Simulate"), icon = icon("play")))),
                                        column(2, style="position:absolute;top:-18px;right:1em", conditionalPanel(condition = "input.method == 3", actionButton(style = "color:white; background-color:purple;font-size:12px;", "go3", strong(" Simulate"), icon = icon("play"))))
                                        ))),
                      
   
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
                                                          infoBoxOutput("dose5", width = 6),
                                                          infoBoxOutput("auc1", width = 6),
                                                          infoBoxOutput("auc2", width = 6),
                                                          infoBoxOutput("auc3", width = 6),
                                                          infoBoxOutput("auc4", width = 6),
                                                          infoBoxOutput("auc5", width = 6)
                                                          ))))))
                      ))
                      
                      
                      
                                 


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  txt <- reactive({
    
    tot_combination <- function(min_dose, max_dose, resolution) {
      
      # Create combinations of dose
      studied_doses <- seq(min_dose, max_dose, resolution)
      
      combination <- as.data.frame(expand.grid(DOSE_1 = studied_doses, DOSE_2 = studied_doses, DOSE_3 = studied_doses, DOSE_4 = studied_doses))
      
      nrow(combination)
      
    }
    
    paste0("Combinations: ", tot_combination(input$min_dose, input$max_dose, input$resolution))
    
    }) 
  
  output$output_text <- renderText({
      txt()
    })
  
  monte_carlo <- eventReactive(input$go1, {
    
    # #########################################
    # # Download models
    # linkMod1 <- "https://github.com/twelvebarsblue/PRGS-Shiny/raw/refs/heads/main/run001.cpp?raw=true"
    # linkCRCL_reference <- "https://raw.githubusercontent.com/twelvebarsblue/PRGS-Shiny/refs/heads/main/CRCL_reference.csv"
    # 
    # download.file(linkMod1, destfile = "run001.cpp")
    # download.file(linkCRCL_reference, destfile = "CRCL_reference.csv", method = "curl")
    # 
    
    #########################################
    # Load models
    
    if (input$model == 1) {
      mod <- mread("C:/Dropbox/R Shiny/model/run001.cpp")
      helpText = desc_m1 }
    else if (input$model == 2) {
      mod <- mread("./run002.cpp")
      helpText = desc_m2 }
    else if (input$model == 3) {
      mod <- mread("./run003.cpp")
      helpText = desc_m3 }
    else if (input$model == 4) {
      mod <- mread("./run004.cpp") 
      helpText = desc_m4 }
    else if (input$model == 5) {
      mod <- mread("./run005.cpp") 
      helpText = desc_m5 }
    
    
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
      
      refRange <- read_csv("./CRCL_reference.csv")
      
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
      
      
      TIME_1 <- ii
      TIME_2 <- TIME_1 + ii
      TIME_3 <- TIME_2 + ii
      TIME_4 <- TIME_3 + ii
    
      
      combination <- as.data.frame(expand.grid(CRCL = crcl, 
                                               AGE = age,
                                               DOSE_1 = studied_doses, DOSE_2 = studied_doses, DOSE_3 = studied_doses, DOSE_4 = studied_doses, DOSE_5 = studied_doses,
                                               II_1 = ii, II_2 = ii, II_3 = ii, II_4 = ii,
                                               TIME_1 = TIME_1, TIME_2 = TIME_2, TIME_3 = TIME_3, TIME_4 = TIME_4,
                                               UNTIL = sim_time))
      
      combination$TAG = seq(nrow(combination))
    

      
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
                      LP_DV=quantile(DV,probs=minprob),
                      HP_DV=quantile(DV,probs=maxprob)
            )  
          
          
          # Check the concentration at 12 hours, if < 15 mg/mL, go next combination
          
          dv = filter(sum_stat, time %in% c(TIME_1, TIME_2, TIME_3, TIME_4))
          
          incProgress(1/n, detail = paste("Combination ", i)) # To check whether its running or not
          
          if (all(dv$Median_C > llim) & max(sum_stat$Median_C) < ulim) {
            
            # cat(paste0("Recommended loading dose: ", df$DOSE_1, " mg", "\n",
            #            "Maintenance dose at 12 h: ", df$DOSE_2, " mg", "\n",
            #            "Maintenance dose at 24 h: ", df$DOSE_3, " mg", "\n",
            #            "Maintenance dose at 36 h: ", df$DOSE_4, " mg", "\n",
            #            "Maintenance dose at 48 h: ", df$DOSE_5, " mg", "\n"))
            
            plot_pk <- ggplot(sum_stat, aes(x=time,y=Median_C)) +
              
              ## Add ribbon for variability
              geom_ribbon(aes(ymin=LP_DV, ymax=HP_DV, x=time), alpha = 0.15, linetype=0)+
              
              geom_ribbon(aes(ymin=llim,ymax=ulim), fill="purple", alpha=0.3) +
              
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
            
            # Clean up the table
            
            cleaned_table <- sum_stat %>% select(time, Median_C)
            cleaned_table$Median_C <- abs(cleaned_table$Median_C)
            cleaned_table$Dose <- 0 # Initiate blank column
            cleaned_table$`Drug Concentration` <- 0
            cleaned_table$AUC <- 0 # Initiate blank column
            
            for (i in 1:nrow(cleaned_table)) {
              time = cleaned_table$time[i]
              
              if (time == 0) {
                cleaned_table$Dose[i] <- 1
                cleaned_table$`Drug Concentration`[i] <- df$DOSE_1
                cleaned_table$AUC[i] <- 0 
              } 
              if (time > 0 & time <= TIME_1) {
                cleaned_table$Dose[i] <- 1
                cleaned_table$`Drug Concentration`[i] <- df$DOSE_1
                cleaned_table$AUC[i] <- pk.calc.auc(conc = cleaned_table$Median_C, 
                                                    time = cleaned_table$time, 
                                                    interval = c(0, TIME_1))
                auc1 <- cleaned_table$AUC[i]
              }
              if (time > TIME_1 & time <= TIME_2) {
                cleaned_table$Dose[i] <- 2
                cleaned_table$`Drug Concentration`[i] <- df$DOSE_2
                cleaned_table$AUC[i] <- pk.calc.auc(conc = cleaned_table$Median_C, 
                                                    time = cleaned_table$time, 
                                                    interval = c(TIME_1, TIME_2))
                auc2 <- cleaned_table$AUC[i]
              }
              if (time > TIME_2 & time <= TIME_3) {
                cleaned_table$Dose[i] <- 3
                cleaned_table$`Drug Concentration`[i] <- df$DOSE_3
                cleaned_table$AUC[i] <- pk.calc.auc(conc = cleaned_table$Median_C, 
                                                    time = cleaned_table$time, 
                                                    interval = c(TIME_2, TIME_3))
                auc3 <- cleaned_table$AUC[i]
              }
              if (time > TIME_3 & time <= TIME_4) {
                cleaned_table$Dose[i] <- 4
                cleaned_table$`Drug Concentration`[i] <- df$DOSE_4
                cleaned_table$AUC[i] <- pk.calc.auc(conc = cleaned_table$Median_C, 
                                                    time = cleaned_table$time, 
                                                    interval = c(TIME_3, TIME_4))
                auc4 <- cleaned_table$AUC[i]
              }
              if (time > TIME_4) {
                cleaned_table$Dose[i] <- 5
                cleaned_table$`Drug Concentration`[i] <- df$DOSE_5
                cleaned_table$AUC[i] <- pk.calc.auc(conc = cleaned_table$Median_C, 
                                                    time = cleaned_table$time, 
                                                    interval = c(TIME_4, sim_time))
                auc5 <- cleaned_table$AUC[i]
              }
            }
            
            colnames(cleaned_table) <- c("Time (h)", "Trough Concentration (mg/mL)", "Regimen", "IV (mg/mL)", "AUC")
            
            cleaned_table <- cleaned_table[, c("Time (h)", "Regimen", "IV (mg/mL)", "Trough Concentration (mg/mL)", "AUC")]
            
            # print(cleaned_table)
            
            crcl <<- crcl
            TIME_1 <<- TIME_1
            TIME_2 <<- TIME_2
            TIME_3 <<- TIME_3
            TIME_4 <<- TIME_4
            initDose <<- df$DOSE_1
            dose2 <<- df$DOSE_2
            dose3 <<- df$DOSE_3
            dose4 <<- df$DOSE_4
            dose5 <<- df$DOSE_5
            auc1 <<- auc1
            auc2 <<- auc2
            auc3 <<- auc3
            auc4 <<- auc4
            auc5 <<- auc5
            plot_pk <<- plot_pk
            cleaned_table <<- cleaned_table
            
            break
          }
        }
      })
    }
    
    runSim(nsamples = input$nsamples, minprob = input$minprob, maxprob = input$maxprob, 
           crcl = crcl, age = input$age, sim_time = input$simtime, ii = input$ii, 
           min_dose = input$min_dose, max_dose = input$max_dose, resolution = input$resolution,
           llim = input$Llim, ulim = input$Ulim)
    
    return(list(helpText = helpText, crcl = crcl, 
                TIME_1 = TIME_1, TIME_2 = TIME_2, TIME_3 = TIME_3, TIME_4 = TIME_4, 
                initDose = initDose, dose2 = dose2, dose3 = dose3, dose4 = dose4, dose5 = dose5,
                auc1 = auc1, auc2 = auc2, auc3 = auc3, auc4 = auc4, auc5 = auc5,
                plot_pk = plot_pk, cleaned_table = cleaned_table
                ))
  })
  
  bayes <- eventReactive(input$go2, {
    
    # #########################################
    # # Download models
    # linkMod1 <- "https://github.com/twelvebarsblue/PRGS-Shiny/raw/refs/heads/main/run001.cpp?raw=true"
    # linkCRCL_reference <- "https://raw.githubusercontent.com/twelvebarsblue/PRGS-Shiny/refs/heads/main/CRCL_reference.csv"
    # 
    # download.file(linkMod1, destfile = "run001.cpp")
    # download.file(linkCRCL_reference, destfile = "CRCL_reference.csv", method = "curl")
    
    #########################################
    # Load models
    
    if (input$model == 1) {
      mod <- mread("C:/Dropbox/R Shiny/model/run001.cpp") }
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
      
      refRange <- read_csv("./CRCL_reference.csv")
      
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
    
    runSim <- function(nsamples, crcl, age, sim_time, ii, min_dose = 100, max_dose = 500, resolution = 100, llim, ulim) {
      
      # Update model with individual variability
      # Get ETA 1 and ETA 2
      my_est <- mod %>% 
        adm_rows(time = 0, amt = input$initDose, DV = NA) %>% 
        obs_rows(time = input$time1, DV = input$dv1) %>% 
        obs_rows(time = input$time2, DV = input$dv2) %>% 
        obs_rows(time = input$time3, DV = input$dv3) %>% 
        add_covariates(CRCL = crcl, AGE = age) %>% 
        mapbayest() %>%
        as.data.frame()
      
      eta1 <- unique(my_est$ETA1)
      eta2 <- unique(my_est$ETA2)
      
      # Create combinations of dose
      studied_doses <- seq(min_dose, max_dose, resolution)
      
      TIME_1 <- ii
      TIME_2 <- TIME_1 + ii
      TIME_3 <- TIME_2 + ii
      TIME_4 <- TIME_3 + ii
      
      combination <- data_frame(expand.grid(CRCL = crcl, AGE = age, 
                                               DOSE_1 = studied_doses, DOSE_2 = studied_doses, DOSE_3 = studied_doses, DOSE_4 = studied_doses, DOSE_5 = studied_doses,
                                               TIME_1 = TIME_1, TIME_2 = TIME_2, TIME_3 = TIME_3, TIME_4 = TIME_4,
                                               UNTIL = sim_time))
      
      
      combination$TAG = seq(nrow(combination))
      
      withProgress(message = "Simulating ", value = 0, {
        
        n = nrow(combination)
      
        for (i in 1:nrow(combination)) {
          
          df <- filter(combination, TAG == i)
          
          out <- mod %>% 
            zero_re() %>% ## Run a simulation, with the new ETA1 and ETA2 values but without any other level of variability
            within(mod, {ETA1 <- eta1}) %>%
            within(mod, {ETA1 <- eta2}) %>%
            mrgsim(idata = df, end = sim_time) %>%
            as.data.frame()
          
          # Check the concentration at 12 hours, if < 15 mg/mL, go next combination
          
          dv = filter(out, time %in% c(TIME_1, TIME_2, TIME_3, TIME_4))
          
          incProgress(1/n, detail = paste("Combination ", i)) # To check whether its running or not
          
          if (all(dv$DV > llim) & max(out$DV) < ulim) {
            # cat(paste0("Recommended loading dose: ", df$DOSE_1, " mg", "\n", 
            #            "Maintenance dose at 12 h: ", df$DOSE_2, " mg", "\n", 
            #            "Maintenance dose at 24 h: ", df$DOSE_3, " mg", "\n",
            #            "Maintenance dose at 36 h: ", df$DOSE_4, " mg", "\n"))
            
            print(dv)
            
            plot_pk <- ggplot(out, aes(x=time,y=DV)) +
              
              geom_ribbon(aes(ymin=llim,ymax=ulim), fill="purple", alpha=0.3) +
              
              ## Add median line
              geom_line(size=2) +
              
              # scale_y_log10()+
              
              # Set axis and theme
              ylab(paste("Concentration",sep=""))+
              xlab("Time after dose (h)")+
              theme_bw()+
              
              # Remove legend
              theme(legend.position="none")
            
            # print(plot_pk)
            
            # Clean up the table
            
            cleaned_table <- out %>% select(time, DV)
            cleaned_table$Dose <- 0 # Initiate blank column
            cleaned_table$`Drug Concentration` <- 0
            cleaned_table$AUC <- 0 # Initiate blank column
            
            for (i in 1:nrow(cleaned_table)) {
              time = cleaned_table$time[i]
              
              if (time == 0) {
                cleaned_table$Dose[i] <- 1
                cleaned_table$`Drug Concentration`[i] <- df$DOSE_1
                cleaned_table$AUC[i] <- 0 
              } 
              if (time > 0 & time <= TIME_1) {
                cleaned_table$Dose[i] <- 1
                cleaned_table$`Drug Concentration`[i] <- df$DOSE_1
                cleaned_table$AUC[i] <- pk.calc.auc(conc = cleaned_table$DV, 
                                                    time = cleaned_table$time, 
                                                    interval = c(0, TIME_1))
                auc1 <- cleaned_table$AUC[i]
              }
              if (time > TIME_1 & time <= TIME_2) {
                cleaned_table$Dose[i] <- 2
                cleaned_table$`Drug Concentration`[i] <- df$DOSE_2
                cleaned_table$AUC[i] <- pk.calc.auc(conc = cleaned_table$DV, 
                                                    time = cleaned_table$time, 
                                                    interval = c(TIME_1, TIME_2))
                auc2 <- cleaned_table$AUC[i]
              }
              if (time > TIME_2 & time <= TIME_3) {
                cleaned_table$Dose[i] <- 3
                cleaned_table$`Drug Concentration`[i] <- df$DOSE_3
                cleaned_table$AUC[i] <- pk.calc.auc(conc = cleaned_table$DV, 
                                                    time = cleaned_table$time, 
                                                    interval = c(TIME_2, TIME_3))
                auc3 <- cleaned_table$AUC[i]
              }
              if (time > TIME_3 & time <= TIME_4) {
                cleaned_table$Dose[i] <- 4
                cleaned_table$`Drug Concentration`[i] <- df$DOSE_4
                cleaned_table$AUC[i] <- pk.calc.auc(conc = cleaned_table$DV, 
                                                    time = cleaned_table$time, 
                                                    interval = c(TIME_3, TIME_4))
                auc4 <- cleaned_table$AUC[i]
              }
              if (time > TIME_4) {
                cleaned_table$Dose[i] <- 5
                cleaned_table$`Drug Concentration`[i] <- df$DOSE_5
                cleaned_table$AUC[i] <- pk.calc.auc(conc = cleaned_table$DV, 
                                                    time = cleaned_table$time, 
                                                    interval = c(TIME_4, sim_time))
                auc5 <- cleaned_table$AUC[i]
              }
            }
            
            colnames(cleaned_table) <- c("Time (h)", "Trough Concentration (mg/mL)", "Regimen", "IV (mg/mL)", "AUC")
            
            cleaned_table <- cleaned_table[, c("Time (h)", "Regimen", "IV (mg/mL)", "Trough Concentration (mg/mL)", "AUC")]
            
            print(cleaned_table)
            
            crcl <<- crcl
            TIME_1 <<- TIME_1
            TIME_2 <<- TIME_2
            TIME_3 <<- TIME_3
            TIME_4 <<- TIME_4
            initDose <<- df$DOSE_1
            dose2 <<- df$DOSE_2
            dose3 <<- df$DOSE_3
            dose4 <<- df$DOSE_4
            dose5 <<- df$DOSE_5
            auc1 <<- auc1
            auc2 <<- auc2
            auc3 <<- auc3
            auc4 <<- auc4
            auc5 <<- auc5
            plot_pk <<- plot_pk
            cleaned_table <<- cleaned_table
            
            break
          }
        }
      })
    }
  
    
    runSim(nsamples = input$nsamples, crcl = crcl, age = input$age,
           sim_time = input$simtime, ii = input$ii_bayes, 
           min_dose = input$min_dose_bayes, max_dose = input$max_dose_bayes, resolution = input$resolution_bayes,
           llim = input$Llim, ulim = input$Ulim)
    
    return(list(crcl = crcl, TIME_1 = TIME_1, TIME_2 = TIME_2, TIME_3 = TIME_3, TIME_4 = TIME_4, 
                initDose = initDose, dose2 = dose2, dose3 = dose3, dose4 = dose4, dose5 = dose5, 
                auc1 = auc1, auc2 = auc2, auc3 = auc3, auc4 = auc4, auc5 = auc5, 
                plot_pk = plot_pk,
                cleaned_table = cleaned_table))
    
  })
  
  manual <- eventReactive(input$go3, {
    
    # #########################################
    # # Download models
    # linkMod1 <- "https://github.com/twelvebarsblue/PRGS-Shiny/raw/refs/heads/main/run001.cpp?raw=true"
    # linkCRCL_reference <- "https://raw.githubusercontent.com/twelvebarsblue/PRGS-Shiny/refs/heads/main/CRCL_reference.csv"
    # 
    # download.file(linkMod1, destfile = "run001.cpp")
    # download.file(linkCRCL_reference, destfile = "CRCL_reference.csv", method = "curl")
    # 
    
    #########################################
    # Load models
    
    if (input$model == 1) {
      mod <- mread("C:/Dropbox/R Shiny/model/run001.cpp") }
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
      
      refRange <- read_csv("./CRCL_reference.csv")
      
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
    
    runSim <- function(nsamples, minprob, maxprob, crcl, age, d1, d2, d3, d4, d5, t1, t2, t3, t4, sim_time, minpro, llim, ulim) {
      
      if (t2 %% t1 == 0 && t3 %% t1 == 0 && t4 %% t1 == 0) {
        ii_1 <- t1
        ii_2 <- t1
        ii_3 <- t1
        ii_4 <- t1 }
      else {
        ii_1 <- t1
        ii_2 <- t2-t1
        ii_3 <- t3-t2
        ii_4 <- t4-t2 }
      
      df <- data_frame(CRCL = crcl, AGE = age, 
                          DOSE_1 = d1, DOSE_2 = d2, DOSE_3 = d3, DOSE_4 = d4, DOSE_5 = d5,
                          II_1 = ii_1, II_2 = ii_2, II_3 = ii_3, II_4 = ii_4,
                          TIME_1 = t1, TIME_2 = t2, TIME_3 = t3, TIME_4 = t4,
                          UNTIL = sim_time)
      
      print(df)
      
      repeated_df <- df[rep(rownames(df), each = nsamples), ]
      
      repeated_df$ID <- seq(1:nsamples)
      
      out <- mrgsim(mod, idata = repeated_df, end = sim_time) %>% as.data.frame()
      
      sum_stat <- out %>%
        group_by(time) %>%
        summarise(Median_C=median(DV),
                  LP_DV=quantile(DV,probs=minprob),
                  HP_DV=quantile(DV,probs=maxprob)
        )  
      
      plot_pk <- ggplot(sum_stat, aes(x=time,y=Median_C)) +
        
        ## Add ribbon for variability
        geom_ribbon(aes(ymin=LP_DV, ymax=HP_DV, x=time), alpha = 0.15, linetype=0)+
        
        geom_ribbon(aes(ymin=llim,ymax=ulim), fill="purple", alpha=0.3) +
        
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
      
      # Clean up the table
      
      cleaned_table <- sum_stat %>% select(time, Median_C)
      cleaned_table$Median_C <- abs(cleaned_table$Median_C)
      cleaned_table$Dose <- 0 # Initiate blank column
      cleaned_table$`Drug Concentration` <- 0
      cleaned_table$AUC <- 0 # Initiate blank column
      
      for (i in 1:nrow(cleaned_table)) {
        time = cleaned_table$time[i]
        
        if (time == 0) {
          cleaned_table$Dose[i] <- 1
          cleaned_table$`Drug Concentration`[i] <- df$DOSE_1
          cleaned_table$AUC[i] <- 0 
        } 
        if (time > 0 & time <= t1) {
          cleaned_table$Dose[i] <- 1
          cleaned_table$`Drug Concentration`[i] <- df$DOSE_1
          cleaned_table$AUC[i] <- pk.calc.auc(conc = cleaned_table$Median_C, 
                                              time = cleaned_table$time, 
                                              interval = c(0, t1))
          auc1 <- cleaned_table$AUC[i]
        }
        if (time > t1 & time <= t2) {
          cleaned_table$Dose[i] <- 2
          cleaned_table$`Drug Concentration`[i] <- df$DOSE_2
          cleaned_table$AUC[i] <- pk.calc.auc(conc = cleaned_table$Median_C, 
                                              time = cleaned_table$time, 
                                              interval = c(t1, t2))
          auc2 <- cleaned_table$AUC[i]
        }
        if (time > t2 & time <= t3) {
          cleaned_table$Dose[i] <- 3
          cleaned_table$`Drug Concentration`[i] <- df$DOSE_3
          cleaned_table$AUC[i] <- pk.calc.auc(conc = cleaned_table$Median_C, 
                                              time = cleaned_table$time, 
                                              interval = c(t2, t3))
          auc3 <- cleaned_table$AUC[i]
        }
        if (time > t3 & time <= t4) {
          cleaned_table$Dose[i] <- 4
          cleaned_table$`Drug Concentration`[i] <- df$DOSE_4
          cleaned_table$AUC[i] <- pk.calc.auc(conc = cleaned_table$Median_C, 
                                              time = cleaned_table$time, 
                                              interval = c(t3, t4))
          auc4 <- cleaned_table$AUC[i]
        }
        if (time > t4) {
          cleaned_table$Dose[i] <- 5
          cleaned_table$`Drug Concentration`[i] <- df$DOSE_5
          cleaned_table$AUC[i] <- pk.calc.auc(conc = cleaned_table$Median_C, 
                                              time = cleaned_table$time, 
                                              interval = c(t4, sim_time))
          auc5 <- cleaned_table$AUC[i]
        }
      }
      
      colnames(cleaned_table) <- c("Time (h)", "Trough Concentration (mg/mL)", "Regimen", "IV (mg/mL)", "AUC")
      
      cleaned_table <- cleaned_table[, c("Time (h)", "Regimen", "IV (mg/mL)", "Trough Concentration (mg/mL)", "AUC")]
      
      crcl <<- crcl
      TIME_1 <<- t1
      TIME_2 <<- t2
      TIME_3 <<- t3
      TIME_4 <<- t4
      initDose <<- df$DOSE_1
      dose2 <<- df$DOSE_2
      dose3 <<- df$DOSE_3
      dose4 <<- df$DOSE_4
      dose5 <<- df$DOSE_5
      auc1 <<- auc1
      auc2 <<- auc2
      auc3 <<- auc3
      auc4 <<- auc4
      auc5 <<- auc5
      plot_pk <<- plot_pk
      cleaned_table <<- cleaned_table
        
    }
    
    runSim(nsamples = input$nsamples, minprob = input$minprob, maxprob = input$maxprob, 
           crcl = crcl, age = input$age,
           d1 = input$init_dose, d2 = input$d1, d3 = input$d2, d4 = input$d3, d5 = input$d4,
           t1 = input$time_1, t2 = input$time_2, t3 = input$time_3, t4 = input$time_4, 
           sim_time = input$simtime,
           llim = input$Llim, ulim = input$Ulim
          )
    
    return(list(crcl = crcl, 
                TIME_1 = TIME_1, TIME_2 = TIME_2, TIME_3 = TIME_3, TIME_4 = TIME_4, 
                initDose = initDose, dose2 = dose2, dose3 = dose3, dose4 = dose4, dose5 = dose5, 
                auc1 = auc1, auc2 = auc2, auc3 = auc3, auc4 = auc4, auc5 = auc5, 
                plot_pk = plot_pk, cleaned_table = cleaned_table))
    
  })
  
  
  #########################################################################################################
  # PDF Output 
  #########################################################################################################
  
  output$downloadData1 <- downloadHandler(
    # For PDF output, change this to "*.pdf"
    filename = "Monte-Carlo.pdf",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "Monte-Carlo.Rmd")
      tempLogo <- file.path(tempdir(), "logo.png")
      tempCover <- file.path(tempdir(), "cover.png")
      file.copy("./Report/Monte-Carlo.Rmd", tempReport, overwrite = TRUE)
      file.copy("./Report/logo.png", tempLogo, overwrite = TRUE)
      file.copy("./Report/cover.png", tempCover, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      params <- list(
        gender = input$gender,
        age = input$age,
        monte_carlo = monte_carlo())
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv()))})
  
  output$downloadData2 <- downloadHandler(
    # For PDF output, change this to "*.pdf"
    filename = "MAP-Bayesian.pdf",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "MAP-Bayesian.Rmd")
      tempLogo <- file.path(tempdir(), "logo.png")
      tempCover <- file.path(tempdir(), "cover.png")
      file.copy("./Report/MAP-Bayesian.Rmd", tempReport, overwrite = TRUE)
      file.copy("./Report/logo.png", tempLogo, overwrite = TRUE)
      file.copy("./Report/cover.png", tempCover, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      params <- list(
        gender = input$gender,
        age = input$age,
        bayes = bayes())
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv()))})
  
  output$downloadData3 <- downloadHandler(
    # For PDF output, change this to "*.pdf"
    filename = "Manual.pdf",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "Manual.Rmd")
      tempLogo <- file.path(tempdir(), "logo.png")
      tempCover <- file.path(tempdir(), "cover.png")
      file.copy("./Report/Manual.Rmd", tempReport, overwrite = TRUE)
      file.copy("./Report/logo.png", tempLogo, overwrite = TRUE)
      file.copy("./Report/cover.png", tempCover, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      params <- list(
        gender = input$gender,
        age = input$age,
        manual = manual())
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv()))})
  
  #########################################################################################################
  # Prediction Tab
  #########################################################################################################
  
  ###############################
  # Monte-Carlo
  ###############################
  observeEvent(monte_carlo(),{
    output$clearance <- renderInfoBox({
      infoBox(
        "Clearance", paste0(round(monte_carlo()$crcl), " mL/min"), icon = tags$i(class = "fas fa-clock-rotate-left", style="font-size: 12px"),
        color = "teal"
      )
    })
    output$initDose <- renderInfoBox({
      infoBox(
        "Initial Dose", paste0(round(monte_carlo()$initDose), " mg"), icon = tags$i(class = "fas fa-syringe", style="font-size: 12px"),
        color = "purple"
      )
    })
    output$dose2 <- renderInfoBox({
      infoBox(
        paste0("Dose at ", monte_carlo()$TIME_1, "h"), paste0(round(monte_carlo()$dose2), " mg"), icon = tags$i(class = "fas fa-syringe", style="font-size: 12px"),
        color = "purple"
      )
    })
    output$dose3 <- renderInfoBox({
      infoBox(
        paste0("Dose at ", monte_carlo()$TIME_2, "h"), paste0(round(monte_carlo()$dose3), " mg"), icon = tags$i(class = "fas fa-syringe", style="font-size: 12px"),
        color = "purple"
      )
    })
    output$dose4 <- renderInfoBox({
      infoBox(
        paste0("Dose at ", monte_carlo()$TIME_3, "h"), paste0(round(monte_carlo()$dose4), " mg"), icon = tags$i(class = "fas fa-syringe", style="font-size: 12px"),
        color = "purple"
      )
    })
    output$dose5 <- renderInfoBox({
      infoBox(
        paste0("Dose at ", monte_carlo()$TIME_4, "h"), paste0(round(monte_carlo()$dose5), " mg"), icon = tags$i(class = "fas fa-syringe", style="font-size: 12px"),
        color = "purple"
      )
    })
    output$auc1 <- renderInfoBox({
      infoBox(
        paste0("AUC", monte_carlo()$TIME_1), paste0(round(monte_carlo()$auc1, 2)), icon = tags$i(class = "fas fa-prescription", style="font-size: 12px"),
        color = "red"
      )
    })
    output$auc2 <- renderInfoBox({
      infoBox(
        paste0("AUC", monte_carlo()$TIME_2), paste0(round(monte_carlo()$auc2, 2)), icon = tags$i(class = "fas fa-prescription", style="font-size: 12px"),
        color = "red"
      )
    })
    output$auc3 <- renderInfoBox({
      infoBox(
        paste0("AUC", monte_carlo()$TIME_3), paste0(round(monte_carlo()$auc3, 2)), icon = tags$i(class = "fas fa-prescription", style="font-size: 12px"),
        color = "red"
      )
    })
    output$auc4 <- renderInfoBox({
      infoBox(
        paste0("AUC", monte_carlo()$TIME_4), paste0(round(monte_carlo()$auc4, 2)), icon = tags$i(class = "fas fa-prescription", style="font-size: 12px"),
        color = "red"
      )
    })
    output$auc5 <- renderInfoBox({
      infoBox(
        paste0("AUC", input$simtime), paste0(round(monte_carlo()$auc5, 2)), icon = tags$i(class = "fas fa-prescription", style="font-size: 12px"),
        color = "red"
      )
    })
  })
  
  ###############################
  # MAP Bayesian
  ###############################
  observeEvent(bayes(),{
    output$clearance <- renderInfoBox({
      infoBox(
        "Clearance", paste0(round(bayes()$crcl), " mL/min"), icon = tags$i(class = "fas fa-clock-rotate-left", style="font-size: 12px"),
        color = "teal"
      )
    })
    output$initDose <- renderInfoBox({
      infoBox(
        "Initial Dose", paste0(round(bayes()$initDose), " mg"), icon = tags$i(class = "fas fa-syringe", style="font-size: 12px"),
        color = "purple"
      )
    })
    output$dose2 <- renderInfoBox({
      infoBox(
        paste0("Dose at ", bayes()$TIME_1, "h"), paste0(round(bayes()$dose2), " mg"), icon = tags$i(class = "fas fa-syringe", style="font-size: 12px"),
        color = "purple"
      )
    })
    output$dose3 <- renderInfoBox({
      infoBox(
        paste0("Dose at ", bayes()$TIME_2, "h"), paste0(round(bayes()$dose3), " mg"), icon = tags$i(class = "fas fa-syringe", style="font-size: 12px"),
        color = "purple"
      )
    })
    output$dose4 <- renderInfoBox({
      infoBox(
        paste0("Dose at ", bayes()$TIME_3, "h"), paste0(round(bayes()$dose4), " mg"), icon = tags$i(class = "fas fa-syringe", style="font-size: 12px"),
        color = "purple"
      )
    })
    output$dose5 <- renderInfoBox({
      infoBox(
        paste0("Dose at ", bayes()$TIME_4, "h"), paste0(round(bayes()$dose5), " mg"), icon = tags$i(class = "fas fa-syringe", style="font-size: 12px"),
        color = "purple"
      )
    })
    output$auc1 <- renderInfoBox({
      infoBox(
        paste0("AUC", bayes()$TIME_1), paste0(round(bayes()$auc1, 2)), icon = tags$i(class = "fas fa-prescription", style="font-size: 12px"),
        color = "red"
      )
    })
    output$auc2 <- renderInfoBox({
      infoBox(
        paste0("AUC", bayes()$TIME_2), paste0(round(bayes()$auc2, 2)), icon = tags$i(class = "fas fa-prescription", style="font-size: 12px"),
        color = "red"
      )
    })
    output$auc3 <- renderInfoBox({
      infoBox(
        paste0("AUC", bayes()$TIME_3), paste0(round(bayes()$auc3, 2)), icon = tags$i(class = "fas fa-prescription", style="font-size: 12px"),
        color = "red"
      )
    })
    output$auc4 <- renderInfoBox({
      infoBox(
        paste0("AUC", bayes()$TIME_4), paste0(round(bayes()$auc4, 2)), icon = tags$i(class = "fas fa-prescription", style="font-size: 12px"),
        color = "red"
      )
    })
    output$auc5 <- renderInfoBox({
      infoBox(
        paste0("AUC", input$simtime), paste0(round(bayes()$auc5, 2)), icon = tags$i(class = "fas fa-prescription", style="font-size: 12px"),
        color = "red"
      )
    })
  })
    
  ###############################
  # Manual
  ###############################
  observeEvent(manual(),{
    output$clearance <- renderInfoBox({
      infoBox(
        "Clearance", paste0(round(manual()$crcl), " mL/min"), icon = tags$i(class = "fas fa-clock-rotate-left", style="font-size: 12px"),
        color = "teal"
      )
    })
    output$initDose <- renderInfoBox({
      infoBox(
        "Initial Dose", paste0(round(manual()$initDose), " mg"), icon = tags$i(class = "fas fa-syringe", style="font-size: 12px"),
        color = "purple"
      )
    })
    output$dose2 <- renderInfoBox({
      infoBox(
        paste0("Dose at ", manual()$TIME_1, "h"), paste0(round(manual()$dose2), " mg"), icon = tags$i(class = "fas fa-syringe", style="font-size: 12px"),
        color = "purple"
      )
    })
    output$dose3 <- renderInfoBox({
      infoBox(
        paste0("Dose at ", manual()$TIME_2, "h"), paste0(round(manual()$dose3), " mg"), icon = tags$i(class = "fas fa-syringe", style="font-size: 12px"),
        color = "purple"
      )
    })
    output$dose4 <- renderInfoBox({
      infoBox(
        paste0("Dose at ", manual()$TIME_3, "h"), paste0(round(manual()$dose4), " mg"), icon = tags$i(class = "fas fa-syringe", style="font-size: 12px"),
        color = "purple"
      )
    })
    output$dose5 <- renderInfoBox({
      infoBox(
        paste0("Dose at ", manual()$TIME_4, "h"), paste0(round(manual()$dose5), " mg"), icon = tags$i(class = "fas fa-syringe", style="font-size: 12px"),
        color = "purple"
      )
    })
    output$auc1 <- renderInfoBox({
      infoBox(
        paste0("AUC", manual()$TIME_1), paste0(round(manual()$auc1, 2)), icon = tags$i(class = "fas fa-prescription", style="font-size: 12px"),
        color = "red"
      )
    })
    output$auc2 <- renderInfoBox({
      infoBox(
        paste0("AUC", manual()$TIME_2), paste0(round(manual()$auc2, 2)), icon = tags$i(class = "fas fa-prescription", style="font-size: 12px"),
        color = "red"
      )
    })
    output$auc3 <- renderInfoBox({
      infoBox(
        paste0("AUC", manual()$TIME_3), paste0(round(manual()$auc3, 2)), icon = tags$i(class = "fas fa-prescription", style="font-size: 12px"),
        color = "red"
      )
    })
    output$auc4 <- renderInfoBox({
      infoBox(
        paste0("AUC", manual()$TIME_4), paste0(round(manual()$auc4, 2)), icon = tags$i(class = "fas fa-prescription", style="font-size: 12px"),
        color = "red"
      )
    })
    output$auc5 <- renderInfoBox({
      infoBox(
        paste0("AUC", input$simtime), paste0(round(manual()$auc5, 2)), icon = tags$i(class = "fas fa-prescription", style="font-size: 12px"),
        color = "red"
      )
    })
  })
  
  #########################################################################################################
  # Plot
  #########################################################################################################
  
  observeEvent(monte_carlo(),{  
    output$simPlot <- renderPlotly({
      withProgress(message = 'Creating plot', detail='Wait...', value = 0.1, expr={
        ggplotly(monte_carlo()$plot_pk) %>%
          config(
            displaylogo = FALSE,
            modeBarButtonsToRemove = c("lasso2d","zoomIn2d", "zoomOut2d","pan2d","autoScale2d","zoom2d","hoverClosestCartesian"))
      
      })
    }) 
  })
  
  observeEvent(bayes(),{  
    output$simPlot <- renderPlotly({
      withProgress(message = 'Creating plot', detail='Wait...', value = 0.1, expr={
        ggplotly(bayes()$plot_pk) %>%
          config(
            displaylogo = FALSE,
            modeBarButtonsToRemove = c("lasso2d","zoomIn2d", "zoomOut2d","pan2d","autoScale2d","zoom2d","hoverClosestCartesian"))
        
      })
    }) 
  })
  
  observeEvent(manual(),{  
    output$simPlot <- renderPlotly({
      withProgress(message = 'Creating plot', detail='Wait...', value = 0.1, expr={
        ggplotly(manual()$plot_pk) %>%
          config(
            displaylogo = FALSE,
            modeBarButtonsToRemove = c("lasso2d","zoomIn2d", "zoomOut2d","pan2d","autoScale2d","zoom2d","hoverClosestCartesian"))
        
      })
    }) 
  })
    
  }

# Run the application 
shinyApp(ui = ui, server = server)
