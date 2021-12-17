###########################################################################
# Know your Distributions
###########################################################################

# Libraries ---------------------------------------------------------------
library(shiny) # shiny

# UI ----------------------------------------------------------------------
ui <- fluidPage(
  
  selectInput(inputId = "Distribution",
              label = "Choose a Distribution", 
              choices = list(normal = "normal", 
                             beta = "beta"
                             ,gamma = "gamma",
                             `inverse gamma` = "inverse gamma"
              )
              
              
  ),
  
  # conditional on distribution chose, make input of parameters
  conditionalPanel(condition = "input.Distribution == 'normal'",
                   
                   numericInput("mean",
                                value = 0,
                                label = ("Choose the mean")),
                   numericInput("sd",
                                value = 1,
                                label = ("Choose the Standard Deviation"))),
  
  
  conditionalPanel(condition = "input.Distribution == 'beta'",
                   
                   numericInput("shape1",
                                value = 1,
                                label = ("Choose the shape1 parameter")),
                   numericInput("shape2",
                                value = 1,
                                label = ("Choose the shape2 parameter"))),
  
  conditionalPanel(condition = "input.Distribution == 'gamma'",
                   
                   numericInput("shape",
                                value = 1,
                                label = ("Choose the shape parameter")),
                   numericInput("rate",
                                value = 1,
                                label = ("Choose the rate parameter"))),
  
  conditionalPanel(condition = "input.Distribution == 'inverse gamma'",
                   
                   numericInput("shape_inv",
                                value = 1,
                                label = ("Choose the shape parameter")),
                   numericInput("rate_inv",
                                value = 1,
                                label = ("Choose the rate parameter"))),
  
  # output plot            
  plotOutput("plot")
  #align = "left"# plot
  
  
  
)
# Server ------------------------------------------------------------------
server <- function(input, output) {
  
  # use renderPlot() function to make output$plot
  output$plot <- renderPlot({
    
    # conditional on input
    if (input$Distribution == "normal"){
      
      plot(density(rnorm(2000000, input$mean, input$sd)), main = "", xlab = "")
      
      
    }else if(input$Distribution == "beta"){
      
      plot(density(rbeta(2000000, input$shape1, input$shape2)), main = "", xlab = "")
      
    }
    
    else if(input$Distribution == "gamma"){
      
      plot(density(rgamma(2000000, input$shape, input$rate)), main = "", xlab = "")
      
      
    }else if(input$Distribution == "inverse gamma"){
      
      plot(density(invgamma::rinvgamma(2000000, input$shape_inv, input$rate_inv)),
           main = "", xlab = "", form = 1, to = 2)
      
    }
    
  })
  
}
# Run App -----------------------------------------------------------------
shinyApp(ui, server)
