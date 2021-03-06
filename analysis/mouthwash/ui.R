
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(
    
  # Application title
    titlePanel("MOUTHWASH: Maximizing Over Unobservables To Help With Adaptive SHrinkage"),
    
  # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            h2("Data"),
            fileInput(inputId = "YinFile", label = "Choose Matrix of Responses (rds file)"),
            fileInput(inputId = "XinFile", label = "Choose Matrix of Covariates (rds file)"),
            h2("MOUTHWASH Settings"),
            numericInput(inputId = "numsv", label = "Number of Confounders", value = -1),
            sliderInput("lambda0",
                        "Regularization Parameter",
                        min = 1,
                        max = 50,
                        value = 10),
            selectInput(inputId = "famethod", label = "Choose the Factor Analysis",
                        choices = c("Heteroscedastic PCA",
                                    "Shrunken Variance PCA",
                                    "Homoscedastic PCA",
                                    "Quasi MLE (with cate)",
                                    "Regularized MLE",
                                    "Moderated Factor Analysis",
                                    "Homoscedastic FLASH",
                                    "Heteroscedastic FLASH",
                                    "No FA Heteroscedastic",
                                    "No FA Shrunken Var",
                                    "No FA Homoscedastic")),
            selectInput(inputId = "likelihood", label = "Choose your likelihood",
                        choices = c("Normal", "t")),
            selectInput(inputId = "mixcompdist", label = "Choose your mixing distribution",
                        choices = c("Normal", "Uniform"))     
        ),
        
                                        # Show a plot of the generated distribution
        mainPanel(
            div(actionButton(inputId = "runnow", label = "Wash Dat Mouth!"), align = "center"),
            textOutput("message"),
            plotOutput("hist_plot")
        )
    )
))
