library(shiny)
library(tidyverse) 
library(rsconnect)
source('utils.R') 

ui <- fluidPage(

  # App title ----
  titlePanel("Effective SARS-CoV-2 test sensitivity"),

  fluidRow(
    column(12, 
      h4("Please wait a few moments for plots to load. Refresh page to restore defaults."),
      h4("Further details are available in the manuscript ", tags$a(href="https://www.medrxiv.org/content/10.1101/2020.10.21.20217042v1", tags$em("Viral dynamics of SARS-CoV-2 infection and the predictive value of repeat testing (doi: 10.1101/2020.10.21.20217042)"), target="_blank")),
      h5("This page is currently in beta; it is updated frequently and some features may be inaccurate. Please direct feedback to:",style="color:red"),
      h5("Stephen Kissler (",tags$a(href="mailto:skissler@hsph.harvard.edu","skissler@hsph.harvard.edu")," | ",tags$a(href="https://twitter.com/StephenKissler","@StephenKissler",target="_blank"),")."),
      h5("Yonatan Grad (",tags$a(href="mailto:ygrad@hsph.harvard.edu","ygrad@hsph.harvard.edu")," | ",tags$a(href="https://twitter.com/yhgrad","@yhgrad",target="_blank"),")."))),

  fluidRow(

    column(3,
      wellPanel(
        h2("Test characteristics"), 
        sliderInput(inputId = "lod",
                  label = "Test limit of detection (Ct)",
                  min = 0,
                  max = 40,
                  value = 40,
                  step=1),
        sliderInput(inputId = "se",
                  label = "Test sensitivity (for viral concentrations above the limit of detection)",
                  min = 0,
                  max = 1,
                  value = 0.99,
                  step=0.01),
        sliderInput(inputId = "inf_ct",
                  label = "Infectiousness threshold (Ct)",
                  min = 0,
                  max = 40,
                  value = 30,
                  step=1),
        )),

    column(3,
      wellPanel(
        h2("Population characteristics"),
        numericInput(inputId = "prev",
                    label="Prevalence of infectious individuals in the population (0-1):",
                    min=0,
                    max=1,
                    value=0.05,
                    step=0.01
          ),
        # sliderInput(inputId = "prev",
        #             label = "Prevalence of infectious individuals in the population:",
        #             min = 0,
        #             max = 1,
        #             value = 0.05,
        #             step=0.01), 
        numericInput(inputId = "n_attendees",
                    label = "Number of event attendees:",
                    min = 0,
                    value = 500,
                    step=1), 
        sliderInput(inputId = "event_duration",
                    label = "Event duration (hours):",
                    min = 0,
                    max = 24,
                    value = 3,
                    step=0.25), 
        )),

    column(2,
      wellPanel(
        h2("Proliferation parameters"),  
        sliderInput(inputId = "wpmean",
                    label = "Proliferation time mean (days):",
                    min = 0,
                    max = 7,
                    value = 3.2,
                    step=0.1), 

        sliderInput(inputId = "wpsd",
                    label = "Proliferation time std. dev. (days):",
                    min = 0,
                    max = 7,
                    value = 2.1,
                    step=0.1), 
        )),

    column(2,
      wellPanel(
        h2("Clearance parameters"),
        sliderInput(inputId = "wrmean",
                    label = "Clearance time mean (days):",
                    min = 0,
                    max = 30,
                    value = 8.0,
                    step=0.1), 

        sliderInput(inputId = "wrsd",
                    label = "Clearance time std. dev. (days):",
                    min = 0,
                    max = 30,
                    value = 5.1,
                    step=0.1),
        )),

    column(2,
      wellPanel(
        h2("Peak Ct parameters"),
        sliderInput(inputId = "peakctmean",
                    label = "Peak Ct mean:",
                    min = 10,
                    max = 40,
                    value = 22.3,
                    step=0.1), 

        sliderInput(inputId = "peakctsd",
                    label = "Peak Ct std. dev.:",
                    min = 0,
                    max = 40,
                    value = 4.2,
                    step=0.1),
        )),
    ),

  fluidRow(

    column(4,
      # h2("Effective sensitivity"),
      plotOutput(outputId = "distPlot")
      ),

    column(4,
      # h2("Expected number of infecttious attendees"), 
      plotOutput(outputId = "ninfPlot")
      ),

    column(4,
      # h2("Distribution of Ct trajectories")
      plotOutput(outputId = "trajectoryPlot")
      ),
    ),

  fluidRow(

    column(12,
      hr(), 
      h2("About this tool"), 
      tags$p(
        "This tool addresses the question: how many infectious people are likely to show up to an event, given a screening test administered ", tags$em("n"), "days prior to the event?",
        style="font-size:18px"
        ),
      tags$p(
        "The effective sensitivity (left-hand plot) is the probability that an individual who would have been infectious at the event gets screened by the pre-event test. It's best to test as close as possible to the event, but high-sensitivity tests (e.g. RT-qPCR) often have a longer turnaround time than lower-sensitivity tests (e.g. rapid antigen tests). This leads to an inherent trade-off between test sensitivity and turnaround time.",
        style="font-size:18px"
        ),
      tags$p(
        "The center plot gives the estimated number of people who will arrive at the event given a test ", tags$em("n"), "days before the event. It accounts for the event size and prevalence of infectious individuals in the population.",
        style="font-size:18px"
        ),
      tags$p(
        "The right-hand plot depicts the typical viral concentration trajectory specified by the inputs (proliferation duration, clearance duration, and peak Ct). It is expressed in terms of the RT-qPCR cycle threshold (Ct). The default values are extracted from the manuscript 'Viral dynamics of SARS-CoV-2 infection and the predictive value of repeat testing' (doi: 10.1101/2020.10.21.20217042) based on viral concentration trajectories from the NBA. These should not be changed unless there is good reason to!",
        style="font-size:18px"
        ),
      tags$p("Some caveats:",style="font-size:18px"),
      tags$ul(
        tags$li("This assumes that all tests are administered at precisely the same time. In reality, tests will be administered within some window of time, leading to a blending of effective sensitivity values and predicted number of infectious individuals at the event."), 
        tags$li("Individuals are assumed to become infectious immediately upon passing the infectiousness threshold. We also use the simplifying assumption that the Ct threshold is equivalent in the proliferation and clearance phases."), 
        tags$li("We assume that prevalence has held roughly steady during the days prior to the event. If prevalence is increasing, there will be more individuals near the start of their infections at the time of the event, leading to reduced effective sensitivity. The opposite will be true if prevalence is decreasing prior to the event."),
        style="font-size:18px"
        ),
      tags$p("This tool is maintained by", tags$b("Stephen Kissler "), "(skissler@hsph.harvard.edu). Ideas belong to the team below, but any mistakes are his. If you find any, please let him know!",
        style="font-size:18px"),
    tags$p("It was developed by:",style="font-size:18px"),
    tags$ul(
        tags$li("Stephen Kissler"), 
        tags$li("Scott Olesen"), 
        tags$li("Joseph Fauver"),
        tags$li("Christina Mack"),
        tags$li("Caroline Tai"),
        tags$li("Kristin Shiue"),
        tags$li("Chaney Kalinich"),
        tags$li("Sarah Jednak"),
        tags$li("Isabel Ott"),
        tags$li("Chantal Vogels"),
        tags$li("Jay Wohlgemuth"),
        tags$li("James Weisberger"),
        tags$li("John DiFiori"),
        tags$li("Deverick Anderson"),
        tags$li("Jimmie Mancell"),
        tags$li("David Ho"),
        tags$li("Nathan D. Grubaugh"),
        tags$li("Yonatan H. Grad"),
        style="font-size:18px"
        ),
    tags$p("For more information, see 'Viral dynamics of SARS-CoV-2 infection and the predictive value of repeat testing' on medRxiv (doi:10.1101/2020.10.21.20217042)",style="font-size:18px"),

      ), 

    )




  # Sidebar layout with input and output definitions ----
  # sidebarLayout(

  #   # Sidebar panel for inputs ----
  #   sidebarPanel(

  #     h2("Test characteristics"), 

  #     sliderInput(inputId = "lod",
  #                 label = "Test limit of detection (Ct)",
  #                 min = 0,
  #                 max = 40,
  #                 value = 40,
  #                 step=1),

  #     sliderInput(inputId = "se",
  #                 label = "Test sensitivity (for viral concentrations above the limit of detection)",
  #                 min = 0,
  #                 max = 1,
  #                 value = 0.99,
  #                 step=0.01),

  #     sliderInput(inputId = "inf_ct",
  #                 label = "Infectiousness threshold (Ct)",
  #                 min = 0,
  #                 max = 40,
  #                 value = 30,
  #                 step=1),

  #     HTML('<hr style="border: 1px solid black;">'),

  #     h2("Proliferation parameters"), 

  #     sliderInput(inputId = "wpmean",
  #                 label = "Proliferation time mean (days):",
  #                 min = 0,
  #                 max = 7,
  #                 value = 3.2,
  #                 step=0.1), 

  #     sliderInput(inputId = "wpsd",
  #                 label = "Proliferation time std. dev. (days):",
  #                 min = 0,
  #                 max = 7,
  #                 value = 2.1,
  #                 step=0.1), 

  #     HTML('<hr style="border: 1px solid black;">'),

  #     h2("Clearance parameters"), 

  #     sliderInput(inputId = "wrmean",
  #                 label = "Clearance time mean (days):",
  #                 min = 0,
  #                 max = 30,
  #                 value = 8.0,
  #                 step=0.1), 

  #     sliderInput(inputId = "wrsd",
  #                 label = "Clearance time std. dev. (days):",
  #                 min = 0,
  #                 max = 30,
  #                 value = 5.1,
  #                 step=0.1),

  #     HTML('<hr style="border: 1px solid black;">'),

  #     h2("Peak Ct parameters"), 

  #     sliderInput(inputId = "peakctmean",
  #                 label = "Peak Ct mean:",
  #                 min = 0,
  #                 max = 40,
  #                 value = 22.3,
  #                 step=0.1), 

  #     sliderInput(inputId = "peakctsd",
  #                 label = "Peak Ct std. dev.:",
  #                 min = 0,
  #                 max = 40,
  #                 value = 4.2,
  #                 step=0.1),



  #   width=3),

  #   # Main panel for displaying outputs ----
  #   mainPanel(

  #     # plotOutput(outputID = "tempPlot"),

  #     # Output: Histogram ----
  #     plotOutput(outputId = "distPlot"),

  #   width=4)
  # )






)
