library(shiny)
library(rstan)

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      numericInput("N", label = "N", value = 10)
    ),
    mainPanel(
      plotOutput("posteriors")
    )
  )
)

server <- function(input, output, session) {
  ## compile stan model
  model <- stan_model(file = "~/1vs2StepBayesianRegSEM/shiny/lm.stan")
  ## draw samples
  draws <- reactive({
    N <- input$N
    sampling(
      object = model,
      data = list(N = N, x = seq_len(N), y = rnorm(N, seq_len(N), 0.1)),
      chains = 2,
      iter = 1000
    )
  })
  ## plot histograms
  output$posteriors <- renderPlot({
    req(draws())
    op <- par(mfrow = c(1, 2), cex = 1.25)
    hist(extract(draws(), "alpha")[[1]], main = bquote("Posterior samples"~alpha), xlab = expression(alpha))
    hist(extract(draws(), "beta")[[1]], main = bquote("Posterior samples"~beta), xlab = expression(beta))
    par(op)
  })
  
}

shinyApp(ui = ui, server = server)