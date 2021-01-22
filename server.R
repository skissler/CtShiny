library(shiny)
library(tidyverse) 
library(rsconnect)
source('utils.R') 

server <- function(input, output) {

  t_test <- seq(-3,0,1/24)


  output$distPlot <- renderPlot({

    trajectory_pars <- c(
      wpmean = input$wpmean,     # Mean last-neg-before-first-pos to peak time
      wpsd = input$wpsd, # SD of last-neg-before-first-pos to peak time
      wrmean = input$wrmean,     # Mean peak to first-neg-after-last-pos time
      wrsd = input$wrsd, # SD of peak to first-neg-after-last-pos time
      dpmean = 40-input$peakctmean,  # Mean Ct difference from lod
      dpsd = input$peakctsd,   # SD of Ct difference from lod
      inf_ct = input$inf_ct,  # Ct threshold for infectiousness 
      maxct = 40    # Maxumum possible Ct, corresponding e.g. to qpcr lod
      )

    eff_se <- unlist(lapply(t_test, get_effective_sensitivity, 
        lod=input$lod, se=input$se, trajectory_pars=trajectory_pars, event_duration=(input$event_duration)/24))

    ggplot(data=tibble(x=-t_test,y=eff_se), aes(x=x, y=y)) + 
      geom_point(size=0.1, alpha=0) + 
      geom_line(stat="smooth", method="loess", span=0.6) + 
      scale_y_continuous(limits=c(0,1)) + 
      scale_x_reverse() + 
      labs(title="Effective sensitivity", x="Days prior to event", y="Effective sensitivity") + 
      theme_minimal() + 
      theme(text=element_text(size=18))

    })

    output$ninfPlot <- renderPlot({

      trajectory_pars <- c(
        wpmean = input$wpmean,    # Mean last-neg-before-first-pos to peak time
        wpsd = input$wpsd, # SD of last-neg-before-first-pos to peak time
        wrmean = input$wrmean,     # Mean peak to first-neg-after-last-pos time
        wrsd = input$wrsd, # SD of peak to first-neg-after-last-pos time
        dpmean = 40-input$peakctmean,  # Mean Ct difference from lod
        dpsd = input$peakctsd,   # SD of Ct difference from lod
        inf_ct = input$inf_ct,  # Ct threshold for infectiousness 
        maxct = 40    # Maxumum possible Ct, corresponding e.g. to qpcr lod
        )

      pop_pars <- c(
        n_attendees=input$n_attendees,
        prev=input$prev
        )

      ninf <- reduce(lapply(t_test, get_n_infectious, 
        lod=input$lod, se=input$se, trajectory_pars=trajectory_pars, pop_pars=pop_pars, event_duration=(input$event_duration)/24), bind_rows) %>% 
      pivot_wider(names_from="statistic", values_from=c("value")) %>%
      mutate(mean_smooth=predict(loess(mean~t, data=., span=0.6))) %>%
      mutate(lwr_smooth=predict(loess(lwr~t, data=., span=0.6))) %>%
      mutate(upr_smooth=predict(loess(upr~t, data=., span=0.6))) %>%
      mutate(t=-t)

      ggplot() + 
        geom_ribbon(
          data=ninf, 
          aes(x=t, ymin=lwr_smooth, ymax=upr_smooth), alpha=0.2, fill="grey") + 
        geom_point(data=ninf, aes(x=t, y=lwr), size=0.1, alpha=0) +
        geom_point(data=ninf, aes(x=t, y=upr), size=0.1, alpha=0) +
        geom_point(data=ninf, aes(x=t, y=mean), size=0.1, alpha=0) + 
        geom_line(data=ninf, aes(x=t, y=mean), stat="smooth", method="loess", span=0.6) + 
        coord_cartesian(ylim=c(0,max(ninf$upr)), expand=FALSE) + 
        scale_x_reverse() + 
        theme_minimal() + 
        theme(text=element_text(size=18)) + 
        labs(title="Number infectious at event", subtitle="(90% pred. interval)", x="Days prior to event", y="Number infectious at event")

    })

    output$trajectoryPlot <- renderPlot({

      trajectory_pars <- c(
        wpmean = input$wpmean,   # Mean last-neg-before-first-pos to peak time
        wpsd = input$wpsd, # SD of last-neg-before-first-pos to peak time
        wrmean = input$wrmean,     # Mean peak to first-neg-after-last-pos time
        wrsd = input$wrsd, # SD of peak to first-neg-after-last-pos time
        dpmean = 40-input$peakctmean,  # Mean Ct difference from lod
        dpsd = input$peakctsd,   # SD of Ct difference from lod
        inf_ct = input$inf_ct,  # Ct threshold for infectiousness 
        maxct = 40,    # Maxumum possible Ct, corresponding e.g. to qpcr lod
        lod=input$lod  # Test limit of detecttion
        )

      make_sample_trajectory(trajectory_pars)
    })

}
