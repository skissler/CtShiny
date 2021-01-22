library(tidyverse) 
library(purrr)
library(truncnorm)

# Vectorize the uniform draw function:
runif_V = Vectorize(runif)

# Draw gamma-distributed values given a mean and sd; 
rgamma_musd <- function(n, mean, sd){
	shape <- (mean^2)/(sd^2)
	rate <- mean/(sd^2)
	out <- rgamma(n, shape=shape, rate=rate) 
	return(out) 
}

qgamma_musd <- function(p, mean, sd){
	shape <- (mean^2)/(sd^2)
	rate <- mean/(sd^2)
	out <- qgamma(p, shape=shape, rate=rate) 
	return(out) 
}

# Randomly sample trajectories (output is a dataframe with n rows)
rtrajectory <- function(n, trajectory_pars, event_duration){
	with(as.list(c(trajectory_pars)),{
 
		# Assign trajectory values to all individuals 
		wpvec <- rgamma_musd(n, wpmean, wpsd)
		wrvec <- rgamma_musd(n, wrmean, wrsd)
		dpvec <- rtruncnorm(n, a=maxct-inf_ct, b=maxct, dpmean, dpsd)

		trajectory_df <- tibble(
			wp=wpvec,
			wr=wrvec,
			dp=dpvec)

		trajectory_df <- trajectory_df %>% 
			mutate(to_lwr = (wr/dp)*(maxct-dp-inf_ct) - wp ) %>% 
			mutate(to_upr = -(wp/dp)*(maxct-inf_ct) + event_duration) %>%
			mutate(to=runif_V(1,to_lwr,to_upr)) %>%
			select(-to_lwr, -to_upr) %>%
			mutate(id=1:n()) %>%
			select(id, to, wp, wr, dp)
	
		return(trajectory_df)	
	})
}

getct <- function(trajectory_df, t, trajectory_pars){

	with(as.list(trajectory_pars),{

	out <- trajectory_df %>% 
		mutate(t=t) %>% 
		mutate(ct=case_when(
			(t>to) & (t<=(to+wp)) ~ maxct - (dp/wp)*(t-to),
			(t>(to+wp)) & (t<=(to+wp+wr)) ~ maxct - (dp-(dp/wr)*(t-(to+wp))),
			TRUE ~ maxct
			))

	return(out) 

	})
}

runtest <- function(ct_df, lod, se){
	out <- ct_df %>% 
		mutate(falseneg=rbinom(n(),1,1-se)) %>% 
		mutate(screened=case_when(
			ct<lod & falseneg==0 ~ 1,
			TRUE ~ 0
		)) %>%
		select(-falseneg)

	return(out)
}

get_effective_sensitivity <- function(t_test, lod, se, trajectory_pars, event_duration=3/24, ndraws=1000){

	trajectory_df <- rtrajectory(ndraws, trajectory_pars, event_duration)

	out <- trajectory_df %>% 
		getct(t_test, trajectory_pars) %>% 
		runtest(lod=lod, se=se) %>%
		summarise(sensitivity=sum(screened)/n()) %>% 
		pull(sensitivity)

	return(out)

}

get_n_infectious <- function(t_test, lod, se, trajectory_pars, pop_pars, event_duration=3/24, ngames=1000, siglevel=0.9){

	with(as.list(c(trajectory_pars, pop_pars)),{

		n_infectious_baseline <- rbinom(ngames, n_attendees, prev)
		eff_se <- get_effective_sensitivity(t_test, lod, se, trajectory_pars, event_duration, ndraws=5000)
		n_infectious <- unlist(lapply(n_infectious_baseline, function(x) rbinom(1, x, 1-eff_se)))
		out <- tibble(
			mean=mean(n_infectious),
			lwr=quantile(n_infectious, (1-siglevel)/2), 
			upr=quantile(n_infectious, 1-(1-siglevel)/2)) %>%
			pivot_longer(everything(), names_to="statistic", values_to="value") %>%
			mutate(t=t_test)

		# n_infectious_baseline <- rbinom(ngames, n_attendees, prev)
		# trajectory_df <- rtrajectory(sum(n_infectious_baseline), trajectory_pars, event_duration)
		# trajectory_df <- mutate(trajectory_df, game=unlist(imap(n_infectious_baseline, ~rep(.y,.x))))

		# out <- trajectory_df %>% 
		# 	getct(t_test, trajectory_pars) %>%
		# 	runtest(lod=lod, se=se) %>% 
		# 	group_by(game) %>% 
		# 	summarise(n_infectious=n()-sum(screened)) %>% 
		# 	ungroup() %>% 
		# 	summarise(
		# 		mean=mean(n_infectious), 
		# 		lwr=quantile(n_infectious, (1-siglevel)/2), 
		# 		upr=quantile(n_infectious, 1-(1-siglevel)/2)) %>%
		# 	pivot_longer(everything(), names_to="statistic", values_to="value") %>% 
		# 	mutate(t=t_test)

		return(out)

		})

}

srise <- function(x, dp, wp){
	out <- dp*(1+x/wp)
	return(out)
}
sfall <- function(x, dp, wr){
	out <- dp*(1-x/wr)
	return(out)
}

make_sample_trajectory <- function(trajectory_pars, siglevel=0.9){
	with(as.list(trajectory_pars),{

	wp_mean <- wpmean
	wp_lwr <- qgamma_musd((1-siglevel)/2, wpmean, wpsd)
	wp_upr <- qgamma_musd(1-(1-siglevel)/2, wpmean, wpsd)

	wr_mean <- wrmean
	wr_lwr <- qgamma_musd((1-siglevel)/2, wrmean, wrsd)
	wr_upr <- qgamma_musd(1-(1-siglevel)/2, wrmean, wrsd)

	dp_mean <- dpmean
	dp_lwr <- qnorm((1-siglevel)/2, dpmean, dpsd)
	dp_upr <- qnorm(1-(1-siglevel)/2, dpmean, dpsd)

	xvals_proliferation <- seq(from=-wp_upr, 0, length.out=500)
	xvals_clearance <- seq(from=0, wr_upr, length.out=500)

	yvals_upr_proliferation <- unlist(lapply(xvals_proliferation, 
	function(x) quantile(srise(x, rnorm(1000,dpmean,dpsd), rgamma_musd(1000,wpmean,wpsd)),1-(1-siglevel)/2)))
	yvals_upr_proliferation_smooth <- predict(loess(yvals~xvals, data=tibble(yvals=yvals_upr_proliferation, xvals=xvals_proliferation), span=0.5))
	yvals_lwr_proliferation <- unlist(lapply(xvals_proliferation, 
	function(x) quantile(srise(x, rnorm(1000,dpmean,dpsd), rgamma_musd(1000,wpmean,wpsd)),(1-siglevel)/2)))
	yvals_lwr_proliferation_smooth <- predict(loess(yvals~xvals, data=tibble(yvals=yvals_lwr_proliferation, xvals=xvals_proliferation), span=0.5))
	
	yvals_upr_clearance <- unlist(lapply(xvals_clearance, 
	function(x) quantile(sfall(x, rnorm(1000,dpmean,dpsd), rgamma_musd(1000,wrmean,wrsd)),1-(1-siglevel)/2)))
	yvals_upr_clearance_smooth <- predict(loess(yvals~xvals, data=tibble(yvals=yvals_upr_clearance, xvals=xvals_clearance), span=0.5))
	yvals_lwr_clearance <- unlist(lapply(xvals_clearance, 
	function(x) quantile(sfall(x, rnorm(1000,dpmean,dpsd), rgamma_musd(1000,wrmean,wrsd)),(1-siglevel)/2)))
	yvals_lwr_clearance_smooth <- predict(loess(yvals~xvals, data=tibble(yvals=yvals_lwr_clearance, xvals=xvals_clearance), span=0.5))

	out <- ggplot() + 
		geom_ribbon(
			data=tibble(
				xvals=xvals_proliferation, 
				yvals_lwr=yvals_lwr_proliferation_smooth, 
				yvals_upr=yvals_upr_proliferation_smooth),
			aes(x=xvals, ymin=maxct-yvals_lwr, ymax=maxct-yvals_upr), alpha=0.2, fill="grey") + 
		geom_segment(aes(x=-wp_mean,xend=0,y=maxct,yend=maxct-dp_mean),col="black") + 
		geom_ribbon(
			data=tibble(
				xvals=xvals_clearance, 
				yvals_lwr=yvals_lwr_clearance_smooth, 
				yvals_upr=yvals_upr_clearance_smooth),
			aes(x=xvals, ymin=maxct-yvals_lwr, ymax=maxct-yvals_upr), alpha=0.2, fill="grey") + 
		geom_segment(aes(x=0,xend=wr_mean,y=maxct-dp_mean,yend=maxct),col="black") + 
		geom_hline(aes(yintercept=inf_ct), linetype="dashed", size=0.2) + 
		geom_hline(aes(yintercept=lod), size=0.2) + 
		geom_text(aes(y=inf_ct),x=0.5,label="Infectiousness threshold",hjust=0,vjust=0) +
		geom_text(aes(y=lod),x=-0.5,label="Limit of detection",hjust=1,vjust=0) +
		coord_cartesian(ylim=c(40,0), expand=FALSE) + 
		theme_minimal() + 
		labs(title="Ct trajectory distribution", subtitle="(90% pred. interval)", x="Days from peak", y="Ct") + 
		scale_y_reverse() + 
		theme(text=element_text(size=18))

	return(out)

	})

}

