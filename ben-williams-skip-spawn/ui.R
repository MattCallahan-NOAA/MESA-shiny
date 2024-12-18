shinyUI(fluidPage(
  tags$head(includeHTML("google-analytics.html")),

  navbarPage("Skip spawning simulator",
             navbarMenu("Instructions",
                        tabPanel("Setup",
                                 mainPanel("This Shiny app is designed to explore different estimation models for fish that exhibit skip spawning." ,
                                           br(),
                                           "These models are then evaluated on how well they accurately reflect the true functional maturity and how any errors propagate through to an estimate of spawning stock biomass (SSB).",
                                           br(),
                                           br(),
                                           "The", em("Input"), "and", em("Models"), "tabs have parameters that can be adjusted in order to explore different maturation rates, variable skip spawning etc.",
                                           br(),
                                           "These simulations are based upon having complete sample sizes for all ages",
                                           br(),
                                           "If no skip spawning is present then a GLM and GAM will produce very similar, if not identical, results."
                                 )
                        )
             ),
             navbarMenu("Inputs",
                        tabPanel("Biological", fluid = TRUE,
                                 sidebarLayout(
                                   sidebarPanel(
                                     tags$div(HTML("<script type='text/x-mathjax-config' >
                                            MathJax.Hub.Config({
                                            tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
                                            });
                                            </script >
                                            ")),
                                     sliderInput("ages", "Max age:", 10, 100, 50),
                                     br(),
                                     "age-length von Bertalanffy parameters",
                                     sliderInput("linf", withMathJax("$L_{\\infty}$"), 10, 100, 80.2),
                                     sliderInput("k", withMathJax("$\\kappa$"), 0.01, 0.4, 0.22, 0.05),
                                     sliderInput("t0", withMathJax("$t_0$"), -5, 2, -0.2, 0.1 ),
                                     "length-weight parameters",
                                     sliderInput("wa1", "Growth parameter a:", min = 0.000001, max = 0.02, value = 0.006),
                                     sliderInput("wa2", "Growth parameter b:", min = 1, max = 4, value = 3.1, 0.1)),
                                   mainPanel(plotOutput("alplot"),
                                             plotOutput("wtplot"))
                                 )
                        ),
                        tabPanel("Sample sizes", fluid = TRUE,
                                 sidebarLayout(
                                   sidebarPanel(
                                     "The number of samples by age.",
                                     sliderInput("nsamp",
                                                 "Sample size:",
                                                 min = 10,
                                                 max = 1000,
                                                 value = 100),
                                     "The number of simulations",
                                     sliderInput("nsim",
                                                 "Number of simulations:",
                                                 min = 10,
                                                 max = 1000,
                                                 value = 200)),
                                   mainPanel(
                                     # dt_output("Samples", "sample_tbl")
                                     )
                                 )
                        ),
                        tabPanel("Skip spawning", fluid = TRUE,
                                 sidebarLayout(
                                   sidebarPanel(
                                     "Skip spawning parameters",
                                     sliderInput("skip", "Percent skip:", min = 0, max = 30, value = 5),
                                     sliderInput("minskip", "Min skip age:", min = 0, max = 100, value = 5),
                                     sliderInput("maxskip", "Max skip age:", min = 0, max = 100, value = 50),
                                     br(),
                                     h4("Only adjust this if a gradient in skip spawning is desired"),
                                     sliderInput("shigh", "% for gradient age:", min = 0, max = 30, value = 5)),

                                   mainPanel(plotOutput("skipplot"))
                                 )
                        ),
                        tabPanel("Population", fluid = TRUE,
                                 sidebarLayout(
                                   sidebarPanel(
                                     sliderInput("mean_pop", "Mean pop (log):", min = 2, max = 15, value =6.5, step = 0.5),
                                     sliderInput("sigr", "Pop variance:", min = 0.1, max = 2, value = 1, step = 0.1),
                                     sliderInput("bin", "Plus age group:", 10, 100, 40),
                                     sliderInput("m", "Natural mortality:", min = 0.01, max = 0.9, value = 0.08),
                                     sliderInput("fish", "Fishing mortality:", min = 0.01, max = 0.9, value = 0.08),
                                     sliderInput("s50", "Age 50% selectivity:", min = 1, max = 15, value = 4.0, step = 0.1)),
                                   mainPanel(plotOutput("popplot")%>% withSpinner(color="#0dc5c1"),
                                             plotOutput("slxplot"))
                                 )
                        )
             ),

             navbarMenu("Models",
                        tabPanel("Age", fluid = TRUE,
                                 sidebarLayout(
                                   sidebarPanel(
                                     "Age-based maturity parameters can be adjusted to reflect a known maturity curve. Additionally the ",
                                     sliderInput("b0", "slope:", min = 0, max = 2, value = 0.84, step = 0.05),
                                     sliderInput("b1", "intercept:", min = 0, max = 15, value = 5.5, step = 0.25),
                                     sliderInput("mat_bin", "Maturity plus age group:", 10, 100, 40),
                                     sliderInput("a_k", "GAM knots:", min = 0, max = 10, value = 10, step = 1)),

                                   mainPanel(plotOutput("matplot") %>% withSpinner(color="#0dc5c1"))
                                 )
                        ),
                        tabPanel("Length", fluid = TRUE,
                                 sidebarLayout(
                                   sidebarPanel(
                                     "Length-based maturity parameters",
                                     sliderInput("b0_l", "slope:", min = 0, max = 2, value = 0.4, step = 0.05),
                                     sliderInput("b1_l", "intercept:", min = 0, max = 100, value = 55, step = 1),
                                     sliderInput("l_k", "GAM knots:", min = 0, max = 10, value = 10, step = 1)),

                                   mainPanel(plotOutput("matplotl") %>% withSpinner(color="#0dc5c1"))
                                 )
                        ),
                        tabPanel("Age-length", fluid = TRUE,
                                 sidebarLayout(
                                   sidebarPanel(
                                     "Age/length-based maturity parameters",
                                     sliderInput("aal_k", "GAM knots - age:", min = 0, max = 10, value = 10, step = 1),
                                     sliderInput("al_k", "GAM knots - length:", min = 0, max = 10, value = 10, step = 1)),

                                   mainPanel(plotOutput("matplotal") %>% withSpinner(color="#0dc5c1"))
                                 )
                        )
                        ,
                        # tabPanel("Weighted-Age", fluid = TRUE,
                        #          sidebarLayout(
                        #            sidebarPanel(
                        #              "Sometimes not all ages have been sampled, this table can be updated with different sample sizes. The associated figure will change along with the table.
                        #       Additionally, data weights can be applied that affect both the GLM and GAM results."),
                        #
                        #            mainPanel(plotOutput("matplotwt") %>% withSpinner(color="#0dc5c1"))
                        #          )
                        # )
             ),
             navbarMenu("Results",
                        tabPanel("Age", fluid = TRUE,
                                 sidebarLayout(
                                   sidebarPanel(
                                     "SSB estimates from (age-based) maturity parameters",
                                     br(),
                                     "If the population or maturity estimator has a 'plus group' then there will be four figures displayed.
                                     On the left is SSB all ages present, on the right is the SSB based upon using a plu group.
                                     Below each figure is a plot of the percent difference from the true SSB ((dotted~line)) compared to the modeled estimates."),
                                   mainPanel(plotOutput("ssb_plot_a") %>% withSpinner(color="#0dc5c1"),
                                             plotOutput("ssb_apd")  %>% withSpinner(color="#0dc5c1"))
                                 )
                        ),
                        tabPanel("Length", fluid = TRUE,
                                 sidebarLayout(
                                   sidebarPanel(
                                     "SSB estimates from (length-based) maturity parameters",
                                     br(),
                                     "If the population or maturity estimator has a 'plus group' then there will be four figures displayed.
                                     On the left is SSB all ages present, on the right is the SSB based upon using a plu group.
                                     Below each figure is a plot of the percent difference from the true SSB ((dotted~line)) compared to the modeled estimates."),
                                   mainPanel(plotOutput("ssb_plot_l") %>% withSpinner(color="#0dc5c1"),
                                             plotOutput("ssb_lpd")  %>% withSpinner(color="#0dc5c1"))
                                 )
                        ),
                        tabPanel("Age/length", fluid = TRUE,
                                 sidebarLayout(
                                   sidebarPanel(
                                     "SSB estimates from (age/length-based) maturity parameters",
                                     br(),
                                     "If the population or maturity estimator has a 'plus group' then there will be four figures displayed.
                                     On the left is SSB all ages present, on the right is the SSB based upon using a plu group.
                                     Below each figure is a plot of the percent difference from the true SSB ((dotted~line)) compared to the modeled estimates."),
                                   mainPanel(plotOutput("ssb_plot_al") %>% withSpinner(color="#0dc5c1"),
                                             plotOutput("ssb_alpd")  %>% withSpinner(color="#0dc5c1"))
                                 )
                        )
                        # ,
                        # tabPanel("Weighted age", fluid = TRUE,
                                 # sidebarLayout(
                                   # sidebarPanel(
                                     # "SSB estimates from weighted age-based maturity parameters"),
                                   # mainPanel(plotOutput("ssb_plot_awt") %>% withSpinner(color="#0dc5c1"),
                                             # plotOutput("ssb_awtpd"))
                                 # )
                        # )
             )
  )
))
