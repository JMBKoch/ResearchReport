###########################################################################
# Shiny App for Visualizing Regularization Priors
###########################################################################

# Libraries ---------------------------------------------------------------
library(shiny) # shiny
library(tidyverse)
library(papaja)

# Preperations ------------------------------------------------------------
# 10000 draws
ndraws <- 5e+05
# sample Small Variance Normal Prior
smallVar <- rnorm(ndraws, mean = 0, sd = 0.1)

# UI ----------------------------------------------------------------------
ui <- fluidPage(
  
  selectInput(inputId = "Distribution",
              label = "Choose a Distribution", 
              choices = list(SVNP = "normal", 
                             RHP = "beta"
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
  

  
  # output plot            
  plotOutput("plot")
  #align = "left"# plot
  
  
  
)
# Server ------------------------------------------------------------------
server <- function(input, output) {
  
  # use renderPlot() function to make output$plot
  output$plot <- renderPlot({
    
 
    
  })
  
}
# Run App -----------------------------------------------------------------
shinyApp(ui, server)







# sample regularized horseshoe prior
regHs <- rep(NA, ndraws)
for(i in 1:ndraws){
  c2 <- LaplacesDemon::rinvgamma(1, shape= nu/2, scale= nu*s2/2)
  lambda <- t(1, scale=1)
  tau <- LaplacesDemon::rhalfcauchy(1, scale=1)
  lambda2_tilde <- c2 * lambda^2/(c2 + tau^2*lambda^2)
  regHs[i] <- rnorm(1, 0, sqrt(tau^2*lambda2_tilde))
}

# make plot
data.frame(dens = c(smallVar, regHs), 
           prior = as.factor(rep(c("Small Variance Normal Prior (\u03c3\u00B2 = 0.01)", "Regularized Horseshoe Prior"), each = ndraws)),
           asymp = rep(0, ndraws)) %>% 
  ggplot(aes(x = dens, fill = prior, linetype = prior)) + 
  geom_density(alpha = .5)+
  geom_vline(aes(xintercept = asymp), linetype = "dashed") +
  xlim(-.5, .5)+
  labs(x = "Size Cross-Loading", title = NULL)+
  theme_apa(legend.pos = "bottom")
