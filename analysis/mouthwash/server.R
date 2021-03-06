
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(succotashr)
library(ggplot2)


shinyServer(function(input, output) {
    
    Y <- reactive({
      if(!is.null(input$YinFile)){
        as.matrix(readRDS(input$YinFile$datapath))
      }
    })
    
    X <- reactive({
      if(!is.null(input$XinFile)){
        as.matrix(readRDS(input$XinFile$datapath))
      }
    })
    
    num_sv <- reactive({
        if (input$numsv < 0) {
            sva::num.sv(dat = t(Y()), mod = model.matrix(~X()))
        } else {
            input$numsv
        }
    })
    
    lambda0 <- reactive({
        input$lambda0
    })
    
    fa_method <- reactive({
        switch(input$famethod,
               "Heteroscedastic PCA"       = "pca",
               "Shrunken Variance PCA"     = "pca_shrinkvar",
               "Homoscedastic PCA"         = "homoPCA",
               "Quasi MLE (with cate)"     = "quasi_mle",
               "Regularized MLE"           = "reg_mle",
               "Moderated Factor Analysis" = "mod_fa",
               "Homoscedastic FLASH"       = "flash",
               "Heteroscedastic FLASH"     = "flash_hetero",
               "No FA Heteroscedastic"     = "non_hetero",
               "No FA Shrunken Var"        = "non_shrinkvar",
               "No FA Homoscedastic"       = "non_homo"
               )
    })
    
    likelihood <- reactive({
        switch(input$likelihood,
               "Normal" = "normal",
               "t"      = "t")
    })
    
    mixcompdist <- reactive({
        switch(input$mixcompdist,
               "Normal"  = "normal",
               "Uniform" = "uniform")
    })
    
    
    
    observeEvent(input$runnow, {
        if(!is.null(Y()) & !is.null(X())){
            suc_out <- succotash(Y = Y(), X = X(), k = num_sv(), fa_method = fa_method(),
                      likelihood = likelihood(), lambda0 = lambda0())
        
        output$hist_plot <- renderPlot(
            qplot(suc_out$betahat, bins = 20, fill = I("skyblue"), color = I("black"), xlab = "betahat")
        )
          } else {
              output$message <- renderText("No mouth to wash")
          }
    })
    
    
    
    
    
})
