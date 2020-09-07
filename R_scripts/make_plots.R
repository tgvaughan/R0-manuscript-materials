## Generate figures for R0 paper
##
## Requires ./extract_sequence_dates.sh be run first and for
## the results of the BEAST analyses to be in ../Results.
##
## Also, be sure that Johns_Hopkins contains the repository described
## in Johns_Hopkins/README, and that PHW contains the CSV file described
## in PHW/README.
##
## Before running, create a directory figures/ to hold the generated
## figures.`


library(tidyverse)
library(lubridate)

outbreaks <- c("Australia",
               "China",
               "The Netherlands (1)",
               "The Netherlands (2)",
               "France (1)",
               "France (2)",
               "Iceland (1)",
               "Iceland (2)",
               "Italy",
               "Spain",
               "WA State (1)",
               "WA State (2)",
               "Iran",
               "Wales",
               "Diamond Princess")

populations <- c("Australia",
                 "China",
                 "The Netherlands",
                 "France",
                 "Iceland",
                 "Italy",
                 "Spain",
                 "WA State",
                 "Iran",
                 "Wales",
                 "Diamond Princess")

outbreak_beast_codes <- c("Australia" = "australia",
                           "China" = "china",
                           "The Netherlands (1)" = "dutch1",
                           "The Netherlands (2)" = "dutch2",
                           "France (1)" = "french1",
                           "France (2)" = "french2",
                           "Iceland (1)" = "iceland1",
                           "Iceland (2)" = "iceland2",
                           "Italy" = "italy",
                           "Spain" = "spain",
                           "WA State (1)" = "washington1",
                           "WA State (2)" = "washington2",
                           "Iran" = "iran",
                           "Wales" = "welsh",
                           "Diamond Princess" = "diamond_princess")

outbreak_population <- c("Australia" = "Australia",
                         "China" = "China",
                         "The Netherlands (1)" = "The Netherlands",
                         "The Netherlands (2)" = "The Netherlands",
                         "France (1)" = "France",
                         "France (2)" = "France",
                         "Iceland (1)" = "Iceland",
                         "Iceland (2)" = "Iceland",
                         "Italy" = "Italy",
                         "Spain" = "Spain",
                         "WA State (1)" = "WA State",
                         "WA State (2)" = "WA State",
                         "Iran" = "Iran",
                         "Wales" = "Wales",
                         "Diamond Princess" = "Diamond Princess")

### Load confirmed case data:

csse_repo_dir <- "Johns_Hopkins/repository/csse_covid_19_data/csse_covid_19_time_series/"

raw_data <- t(read_csv(paste0(csse_repo_dir,"time_series_covid19_confirmed_global.csv"),
                       col_names=FALSE))

confirmed_cases <- NULL
for (population in c("Australia",
                     "China",
                     "France",
                     "Italy",
                     "Spain",
                     "Diamond Princess",
                     "Iceland",
                     "Netherlands",
                     "Iran")) {
    print(population)
    popCol <- which(raw_data[2,] == population)
    print(popCol)
    cases <- raw_data[-(1:4),popCol]
    if (length(dim(cases))>1) {
        mode(cases) <- "numeric"
        cases <- rowSums(cases)
    } else {
        cases <- as.numeric(cases)
    }
    if (population == "Netherlands")
        population <- "The Netherlands"
    if (population == "Australia")
        print(cases)
    confirmed_cases <- bind_rows(confirmed_cases,
                                 tibble(population=population,
                                        date=mdy(raw_data[-(1:4),1]),
                                        confirmed_cumulativeInfections=cases))
}

raw_data <- t(read_csv(paste0(csse_repo_dir, "time_series_covid19_confirmed_US.csv"),
                       col_names=FALSE))
cols <- which(raw_data[7,] == "Washington")
cases <- raw_data[-(1:11),cols]
mode(cases) <- "numeric"
cases <- rowSums(cases)
confirmed_cases <- bind_rows(confirmed_cases,
                             tibble(population="WA State",
                                    date=mdy(raw_data[-(1:11),1]),
                                    confirmed_cumulativeInfections=cases))

welsh_data_file <- "PHW/Wales_cumulative_cases.csv"
confirmed_cases <- bind_rows(confirmed_cases,
                             read_csv(welsh_data_file, col_types="cciinii") %>%
                             mutate(date=dmy(`Specimen date`)) %>%
                             group_by(date) %>%
                             summarize(confirmed_cumulativeInfections=sum(`Cumulative cases`)) %>%
                             mutate(population="Wales"))


# Load parameter logs

load_and_combine_logs <- function(logfiles, burnin=0.1) {

    print("Loading logfiles:")
    print(logfiles)

    df <- NULL

    for (logfile in logfiles) {

        thisdf <- read_tsv(logfile, comment = "#", col_types = cols(Sample = "i"))
        N <- dim(thisdf)[1]
        thisdf <- thisdf[-(1:round(burnin*N)),]

        if (is.null(df)) {
            df <- thisdf
        } else {
            df <- bind_rows(df, thisdf)
        }
    }

    return(df)
}

process_results <- function(logfile_pattern, rateshift="none", burnin=0.1,
                            rR0prior=function(n) {rlnorm(n, 0.8, 0.5)}) {

    logfiles <- dir("Results", logfile_pattern, full.names=TRUE)

    if (length(logfiles) == 0)
        return(NULL)

    df <- load_and_combine_logs(logfiles, burnin=burnin)
    
    out_df <- NULL
    for (outbreak in outbreaks) {

        if (rateshift != "none" && outbreak == "Diamond Princess")
            next

        print(outbreak)

        switch(rateshift,
               none = {
                   Re <- df[[paste0("Re_", outbreak_beast_codes[[outbreak]])]]

                   if (outbreak == "Diamond Princess")
                       sProp <- df[[paste0("sampProp_", outbreak_beast_codes[[outbreak]])]]
                   else
                       sProp <- df[[paste0("sampProp_", outbreak_beast_codes[[outbreak]], "2")]]
                   out_df <- rbind(out_df,
                                   data.frame(Outbreak=outbreak,
                                              Population=outbreak_population[outbreak],
                                              Re=Re,
                                              samplingProp=sProp))
               },

               R0 = {
                   Re <- df[[paste0("Re_", outbreak_beast_codes[[outbreak]], "[0]")]]
                   Re_post <- df[[paste0("Re_", outbreak_beast_codes[[outbreak]], "[1]")]]
                   Re_splitTime <- df[[paste0("ReSplitTime_", outbreak_beast_codes[[outbreak]])]]
                   origin <- df[[paste0("origin_", outbreak_beast_codes[[outbreak]])]]

                   sProp <- df[[paste0("sampProp_", outbreak_beast_codes[[outbreak]], "2")]]
                   out_df <- rbind(out_df,
                                   data.frame(Outbreak=outbreak,
                                              Population=outbreak_population[outbreak],
                                              Re=Re,
                                              Re_post=Re_post,
                                              Re_splitTime=Re_splitTime,
                                              Re_splitTimeRel=Re_splitTime/origin,
                                              samplingProp=sProp))
                   
               },
               
               sampProp = {
                   Re <- df[[paste0("Re_", outbreak_beast_codes[[outbreak]])]]
                   sProp <-  df[[paste0("sampProp_", outbreak_beast_codes[[outbreak]], "[1]")]]
                   sProp_post <-  df[[paste0("sampProp_", outbreak_beast_codes[[outbreak]], "[2]")]]
                   sProp_splitTime <-  df[[paste0("sampChangeTimes_", outbreak_beast_codes[[outbreak]], "2")]]
                   sProp_splitTime_max <-  df[[paste0("sampChangeTimes_", outbreak_beast_codes[[outbreak]], "3")]]

                   out_df <- rbind(out_df,
                                   data.frame(Outbreak=outbreak,
                                              Population=outbreak_population[outbreak],
                                              Re=Re,
                                              samplingProp=sProp,
                                              samplingProp_post=sProp_post,
                                              samplingProp_splitTime=sProp_splitTime,
                                              samplingProp_splitTimeRel=sProp_splitTime/sProp_splitTime_max))

               },

               both = {
                   Re_d <- df[[paste0("Re_", outbreak_beast_codes[[outbreak]], "1")]]
                   Re_post_d <- df[[paste0("Re_", outbreak_beast_codes[[outbreak]], "2")]]
                   sampProp_d <- df[[paste0("sampProp_", outbreak_beast_codes[[outbreak]], "2")]]
                   sampProp_post_d <- df[[paste0("sampProp_", outbreak_beast_codes[[outbreak]], "3")]]
                   out_df <- rbind(out_df,
                                   data.frame(Outbreak=outbreak,
                                              Population=outbreak_population[outbreak],
                                              Re=Re_d,
                                              Re_post=Re_post_d,
                                              samplingProp=sampProp_d,
                                              samplingProp_post=sampProp_post_d))
               }
               )
    }

    print("Prior")

    switch(rateshift,
           none = {
               out_df <- rbind(out_df,
                               data.frame(Outbreak="Prior",
                                          Population=NA,
                                          Re=rR0prior(10000),
                                          samplingProp=rbeta(10000, 1, 4)))

           },
           R0 = {
               base_prior <- rR0prior(5000)
               out_df <- rbind(out_df,
                               data.frame(Outbreak="Prior",
                                          Population=NA,
                                          Re=c(base_prior, base_prior),
                                          Re_post=c(base_prior, rR0prior(5000)),
                                          Re_splitTime=NA,
                                          Re_splitTimeRel=runif(10000, 0, 1),
                                          samplingProp=rbeta(10000, 1, 4)))

           },
           sampProp = {

               base_prior <- rbeta(5000, 1, 4)
               out_df <- rbind(out_df,
                               data.frame(Outbreak="Prior",
                                          Population=NA,
                                          Re=rR0prior(10000),
                                          samplingProp=c(base_prior, base_prior),
                                          samplingProp_post=c(base_prior, rbeta(5000, 1, 4)),
                                          samplingProp_splitTime=NA,
                                          samplingProp_splitTimeRel=runif(10000, 0, 1)))

           },

           both = {
               out_df <- rbind(out_df,
                               data.frame(Outbreak="Prior",
                                          Population=NA,
                                          Re=rR0prior(10000),
                                          Re_post=rR0prior(10000),
                                          samplingProp=rbeta(10000, 1, 4),
                                          samplingProp_post=rbeta(10000, 1, 4)))

           }
           )
    
    out_df$Outbreak <- factor(out_df$Outbreak, levels=c(outbreaks, "Prior"))

    return(out_df)
}


results_indep <- process_results("BD_indep.clock_8e-4.bu_36.5.[1-5].log")

results_altPrior_indep <- process_results("BD_altPrior_indep.clock_8e-4.bu_36.5.[1-5].log",
                                          rR0prior=function(n) { runif(n, 0, 10) })

results_fixedrateshift_indep <- process_results("BD_fixedrateshift_indep.[1-5].log",
                                                rateshift="both")

results_linked <- process_results("BD_linked.clock_8e-4.bu_36.5.[1-5].log")

results_indep_noseq <- process_results("BD_indep_noseq.clock_8e-4.bu_36.5.[1-5].log")


## Estimate R0 values from sample times

outbreak_Re_estimates <- NULL
outbreak_Re_stderrs <- NULL

df_Re <- NULL
for (outbreak in outbreaks) {

    datefile <- paste0(outbreak_beast_codes[[outbreak]],"_dates.txt")
    dates <- sort(ymd(read.table(datefile, header=F)[[1]]))

    x <- dates-dates[1]
    y <- log(1:length(dates))

    res <- lm(y~x)
    lambda_minus_mu <- res$coefficients[2]
    lambda_minus_mu_stderr <- summary(res)$coefficients[2,2]


    bu <- 1/10 # per day

    R0 <- lambda_minus_mu/bu + 1
    R0_stderr <- lambda_minus_mu_stderr/bu + 1

    outbreak_Re_estimates[[outbreak]] <- R0
    outbreak_Re_stderrs[[outbreak]] <- R0_stderr

    df_Re <- rbind(df_Re,
                   data.frame(Outbreak=outbreak,Re_est=R0, Re_est_err=R0_stderr))
}


## Include R0 values from EpiEstim

df_Re_epiEstim <- read_csv("R_estimates_clusters.csv",
                            col_types="cnnn",
                            quote = "") %>%
     group_by(cluster_id) %>%
     mutate(Outbreak=names(outbreak_beast_codes[outbreak_beast_codes==cluster_id]))


## Plot unique Re count derived from DPP analysis:

p <- ggplot(results_linked_unprocessed) +
    geom_histogram(aes(x=uniqueReCount, y=..density..),
                   breaks=seq(1.5,8.5,by=1),
                   col="black",
                   fill="orange") +
    theme_bw() +
    xlab(expression(paste("Unique ",R[e], " values"))) + ylab("Probability")
ggsave("figures/unique_Re_count.png", p, width=15, height=15, units="cm")


## Sample date distribution plot ##

outbreak_dates <- NULL
for (outbreak in outbreaks) {
    datefile <- paste0("sequence_dates/",outbreak_beast_codes[[outbreak]],"_dates.txt")
    dates <- sort(ymd(read.table(datefile, header=F)[[1]]))
    outbreak_dates <- rbind(outbreak_dates,
                            data.frame(Outbreak = outbreak,
                                       Date=dates,
                                       Offset=(dates - min(dates))))
}

p <- ggplot(outbreak_dates) +
    stat_summary(geom="linerange",aes(x=Offset, y=Outbreak, colour=Outbreak),
                 fun.min=min, fun.max=max, fun=mean) +
    geom_point(aes(x=Offset, y=Outbreak, color=Outbreak),
               position=position_jitter(width=0, height=0.2), alpha=0.5) +
    guides(colour=FALSE) +
    xlab("Days following first sampled sequence") +
    ylab("Outbreak")
ggsave("figures/sample_ranges.png", p, width=25, height=15, units="cm")


## Sample count table ##
sink("Sequence_dates_and_counts.txt")
outbreak_dates %>%
    group_by(Outbreak) %>%
    summarize("Final date"=max(Date), "Sequence count"=length(Date))
sink()


## Plot marginal R0 posteriors

## Indep
p <- ggplot(results_indep) +
    geom_violin(aes(Outbreak, Re, fill=Population),
                draw_quantiles = c(0.025, 0.5, 0.975)) + 
    coord_cartesian(ylim=c(0,10)) +
    geom_hline(yintercept=1, linetype="dashed") +
    ylab(expression(R[0])) +
    theme_bw() + guides(fill=FALSE) +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave(paste0("figures/Re_comparison_indep.png"), p,
       width=15, height=10, units="cm")

## Linked
p <- ggplot(results_linked) +
    geom_violin(aes(Outbreak, Re, fill=Population),
                draw_quantiles = c(0.025, 0.5, 0.975)) + 
    coord_cartesian(ylim=c(0,10)) +
    geom_hline(yintercept=1, linetype="dashed") +
    ylab(expression(R[0])) +
    theme_bw() + guides(fill=FALSE) +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave(paste0("figures/Re_comparison_linked.png"), p,
       width=15, height=10, units="cm")

## Alternate R0 prior
p <- ggplot(results_altPrior_indep) +
    geom_violin(aes(Outbreak, Re, fill=Population),
                draw_quantiles = c(0.025, 0.5, 0.975)) + 
    coord_cartesian(ylim=c(0,10)) +
    geom_hline(yintercept=1, linetype="dashed") +
    ylab(expression(R[0])) +
    theme_bw() + guides(fill=FALSE) +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave(paste0("figures/Re_comparison_altPrior_indep.png"), p,
       width=15, height=10, units="cm")

# Sample and R0 rate shift (fixed change time)

p <- ggplot(results_fixedrateshift_indep) +
    geom_violin(aes(Outbreak, Re, fill=Population),
                draw_quantiles = c(0.025, 0.5, 0.975)) + 
    coord_cartesian(ylim=c(0,10)) +
    geom_hline(yintercept=1, linetype="dashed") +
    ylab(expression(R[0])) +
    theme_bw() + guides(fill=FALSE) +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0)) #+
    ## ggtitle("Fixed time rate shift (s and R0)")
ggsave(paste0("figures/Re_comparison_fixedrateshift.png"), p,
       width=15, height=10, units="cm")


## Plot sampling proportions

## Indep
p <- ggplot(results_indep %>% filter(Outbreak != "Prior")) +
    geom_violin(aes(Outbreak, samplingProp, fill=Population),
                scale="width",
                draw_quantiles = c(0.025, 0.5, 0.975)) +
    ylab("Sampling proportion") +
    guides(fill=FALSE) + theme_bw() +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave(paste0("figures/samplingProp_comparison_indep.png"), p,
       width=15, height=10, units="cm")


## Linked
p <- ggplot(results_linked %>% filter(Outbreak != "Prior")) +
    geom_violin(aes(Outbreak, samplingProp, fill=Population),
                scale="width",
                draw_quantiles = c(0.025, 0.5, 0.975)) +
    ylab("Sampling proportion") +
    guides(fill=FALSE) + theme_bw() +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave(paste0("figures/samplingProp_comparison_linked.png"), p,
       width=15, height=10, units="cm")


## Plot marginal R0 posteriors before and after DP quarantine:

df <- load_and_combine_logs(dir("../Results", "BD_indep.clock_8e-4.bu_36.5.[1-5].log",
                                full.names=TRUE))
dfRe <- bind_rows(df %>% select(Re_diamond_princess) %>%
                  rename(Re=Re_diamond_princess) %>%
                  mutate(Phase = "Pre-quarantine"),
                  df %>% select(Re_diamond_princess_postQ) %>%
                  rename(Re=Re_diamond_princess_postQ) %>%
                  mutate(Phase = "Post-quarantine"))
p <- ggplot(dfRe) +
    geom_violin(aes(factor(Phase, levels=c("Pre-quarantine", "Post-quarantine")), Re)) +
    ylab("Reproductive number") + xlab(element_blank())
ggsave(paste0("figures/Re_comparison_DP_pre_post_indep.png"), p,
       width=10, height=10, units="cm")


## Write a table containing R0 and sampling proportion quantiles of indep analysis

sink("Re_and_samplingProp_indep.txt")
results_indep %>%
    group_by(Outbreak) %>%
    summarize(Re_med = median(Re), Re_low = quantile(Re, 0.025), Re_high = quantile(Re, 0.975),
              s_med = median(samplingProp), s_low = quantile(samplingProp, 0.025),
              s_high = quantile(samplingProp, 0.975))
sink()



## Supplemental plot showing marginal Re posteriors
## alongside basic regression estimates and epiestim results

p <- ggplot(results_indep_noseq) +
    geom_violin(aes(Outbreak, Re, fill="BDSKY (no sequences)"),
                draw_quantiles = c(0.025, 0.5, 0.975),
                position = position_nudge(x=-0.2)) + 
    ## scale_y_continuous(limits=c(0,10)) +
    geom_pointrange(aes(Outbreak, Re_est,
                        ymin=Re_est-2*Re_est_err,
                        ymax=Re_est+2*Re_est_err,
                        colour="Regression"),
                    df_Re %>% filter(Outbreak != "Diamond Princess"),
                    position = position_nudge(x=+0.0)) +
    geom_pointrange(aes(Outbreak, R_mean,
                        ymin=R_lowQt,
                        ymax=R_highQt, colour="EpiEstim"),
                    df_Re_epiEstim %>% filter(Outbreak != "Diamond Princess"),
                    position = position_nudge(x=+0.2)) +
    geom_hline(yintercept=1, linetype="dashed") +
    coord_cartesian(ylim=c(0,10)) +
    ylab(expression(R[0])) +
    scale_colour_manual(name=element_blank(),
                        values=c("Regression"="blue",
                                 "EpiEstim"="red")) +
    scale_fill_manual(name=element_blank(), values=c("BDSKY (no sequences)"="orange")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave("figures/Re_comparison_indep_noseq_regression.png", p,
       width=20, height=10, units="cm")



## Case count marginal posteriors (from EpiInf)

source("trajDataTools.R")

casecountdf <- NULL
for (population in unique(outbreak_population)) {
    print(population)

    popdf <- NULL
    
    for (outbreak in names(outbreak_population[outbreak_population==population])) {

        print(outbreak)

        popdf <- bind_rows(popdf,
                           loadTrajectories(paste0("../Results/BD_indep.",
                                                   outbreak_beast_codes[outbreak],
                                                   ".clock_8e-4.bu_36.5.1.traj")) %>%
                           subset(age==0) %>%
                           distinct(traj, .keep_all=TRUE))
    }

    casecountdf <- bind_rows(casecountdf,
                             popdf %>%
                             group_by(traj) %>%
                             summarize(cumulativeInfections=sum(cumulativeInfections)) %>%
                             mutate(population=population))
}

casecountdf$population <- factor(casecountdf$population, levels=populations)

final_confirmed_cases <- NULL
final_confirmed_cases_offset <- NULL
for (population in unique(outbreak_population)) {
    outbreak <- names(outbreak_population[outbreak_population==population])[1]
    final_sample_date <- max((outbreak_dates %>% filter(Outbreak == outbreak))$Date)
    final_confirmed_cases <- bind_rows(final_confirmed_cases,
                                       confirmed_cases %>%
                                       filter(population == !!population) %>%
                                       filter(date == final_sample_date))
    final_confirmed_cases_offset <- bind_rows(final_confirmed_cases_offset,
                                       confirmed_cases %>%
                                       filter(population == !!population) %>%
                                       filter(date == final_sample_date+10))
}

p <- ggplot(casecountdf) +
    geom_violin(aes(population, cumulativeInfections, fill=population),
                draw_quantiles=c(0.025,0.5,0.975)) +
    geom_point(aes(population, confirmed_cumulativeInfections),
               final_confirmed_cases, shape=5, size=2.5, stroke=2, color="white") +
    geom_point(aes(population, confirmed_cumulativeInfections),
               final_confirmed_cases, shape=5, size=3.0, stroke=1, color="red") +
    scale_y_log10() +
    ylab("Total cases") +
    xlab("Population") +
    guides(fill=FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave("figures/caseCounts_comparison_indep.png", p,
       width=15, height=10, units="cm")

p <- ggplot(casecountdf) +
    geom_violin(aes(population, cumulativeInfections, fill=population),
                draw_quantiles=c(0.025,0.5,0.975)) +
    geom_point(aes(population, confirmed_cumulativeInfections),
               final_confirmed_cases_offset, shape=5, size=2.5, stroke=2, color="white") +
    geom_point(aes(population, confirmed_cumulativeInfections),
               final_confirmed_cases_offset, shape=5, size=3.0, stroke=1, color="red") +
    scale_y_log10() +
    ylab("Total cases") +
    xlab("Population") +
    guides(fill=FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave("figures/caseCounts_comparison_indep_offset.png", p,
       width=15, height=10, units="cm")


## Write a table containing case count estimates

sink("case_count_indep.txt")
casecountdf %>%
    filter(indep==TRUE) %>%
    group_by(population) %>%
    summarize(median = median(cumulativeInfections),
              low = quantile(cumulativeInfections, 0.025),
              high = quantile(cumulativeInfections, 0.975)) %>%
    left_join(final_confirmed_cases %>% rename(confirmed=confirmed_cumulativeInfections),
              by="population") %>%
    select(population, median, low, high, confirmed, date) %>%
    mutate(median/confirmed)
sink()

sink("case_count_indep_offset.txt")
casecountdf %>%
    filter(indep==TRUE) %>%
    group_by(population) %>%
    summarize(median = median(cumulativeInfections),
              low = quantile(cumulativeInfections, 0.025),
              high = quantile(cumulativeInfections, 0.975)) %>%
    left_join(final_confirmed_cases_offset %>% rename(confirmed=confirmed_cumulativeInfections),
              by="population") %>%
    select(population, median, low, high, confirmed, date) %>%
    mutate(median/confirmed)
sink()


### Trajectory plots

source("trajDataTools.R")

makeTrajPlot <- function(filenames, final_sample_date, casecount_data=NA, offset=0) {
    if (length(filenames)>1) {

        df <- NULL
        for (thisfname in filenames) {
            df <- bind_rows(df, loadTrajectories(thisfname) %>% mutate(f=thisfname))
        }

        max_age <- max(df$age)

        df_grid <- NULL
        for (thisfname in filenames) {
            df_grid <- bind_rows(df_grid,
                                 gridTrajectories(subset(df,f==thisfname),
                                                  seq(0,max_age,length.out=201)))
        }
        
        df_grid <- df_grid %>%
            group_by(age, traj) %>%
            summarize(cumulativeInfections = sum(cumulativeInfections), .groups="drop_last")

        df <- df_grid

    } else {

        df <- loadTrajectories(filenames)
        df_grid <- gridTrajectories(df, seq(0, max(df$age), length.out=101))
    }

    df$date <- final_sample_date - df$age*366
    df_grid$date <- final_sample_date - df_grid$age*366
    p <- ggplot(df, aes(date, cumulativeInfections)) +
        geom_step(aes(group=traj), alpha=0.05) +
        geom_ribbon(aes(date,cimin, ymin=cimin, ymax=cimax),
                    df_grid %>%
                    group_by(date) %>%
                    summarize(cimin=quantile(cumulativeInfections, 0.025),
                              cimax=quantile(cumulativeInfections, 0.975), .groups="drop_last"),
                    alpha=0.25, fill="red")
    if (length(casecount_data)>0) {
        p <- p + geom_point(aes(date+offset, confirmed_cumulativeInfections), casecount_data,
                            shape=5, size=1, stroke=2, color="white") +
            geom_point(aes(date+offset, confirmed_cumulativeInfections), casecount_data,
                       shape=5, size=1.5, stroke=1, color="red")
    }
    p <- p + geom_line(aes(date,median),
                       df_grid %>%
                       group_by(date) %>%
                       summarize(median=median(cumulativeInfections), .groups="drop_last"),
                       size=1.0, color="orange") +
        xlab(element_blank()) + ylab("Cumulative Infections") + 
        theme_bw()
    return(p)
}

ggsave("figures/trajectories_australia.png",
       makeTrajPlot("../Results/BD_indep.australia.clock_8e-4.bu_36.5.1.traj",
                    ymd("2020-03-11"),
                    confirmed_cases %>% filter(population=="Australia",
                                               confirmed_cumulativeInfections>0)) +
       ggtitle("Australia") +
       scale_x_date(limits = c(ymd("2020-01-01"),ymd("2020-03-11"))) +
       scale_y_log10(limits=c(1,1e4)),
       width=15, height=10, units="cm")

ggsave("figures/trajectories_china.png",
       makeTrajPlot("../Results/BD_indep.china.clock_8e-4.bu_36.5.1.traj",
                    ymd("2020-01-23"),
                    confirmed_cases %>% filter(population=="China",
                                               confirmed_cumulativeInfections>0)) +
       ggtitle("China") +
       scale_x_date(limits = c(ymd("2019-11-15"),ymd("2020-01-23"))) +
       scale_y_log10(limits = c(1,1e5)),
       width=15, height=10, units="cm")

ggsave("figures/trajectories_netherlands.png",
       makeTrajPlot(c("../Results/BD_indep.dutch1.clock_8e-4.bu_36.5.1.traj",
                      "../Results/BD_indep.dutch2.clock_8e-4.bu_36.5.1.traj"),
                    ymd("2020-03-12"),
                    confirmed_cases %>% filter(population=="The Netherlands",
                                               confirmed_cumulativeInfections>0)) + 
       ggtitle("The Netherlands") +
       scale_x_date(limits = c(ymd("2020-01-01"),ymd("2020-03-12"))) +
       scale_y_log10(limits = c(1,1e4)),
       width=15, height=10, units="cm")

ggsave("figures/trajectories_france.png",
       makeTrajPlot(c("../Results/BD_indep.french1.clock_8e-4.bu_36.5.1.traj",
                      "../Results/BD_indep.french2.clock_8e-4.bu_36.5.1.traj"),
                    ymd("2020-03-16"),
                    confirmed_cases %>% filter(population=="France",
                                               confirmed_cumulativeInfections>0)) + 
       ggtitle("France") +
       scale_x_date(limits = c(ymd("2020-01-1"),ymd("2020-03-16"))) +
       scale_y_log10(limits = c(1,1e4)),
       width=15, height=10, units="cm")

ggsave("figures/trajectories_iceland.png",
       makeTrajPlot(c("../Results/BD_indep.iceland1.clock_8e-4.bu_36.5.1.traj",
                      "../Results/BD_indep.iceland2.clock_8e-4.bu_36.5.1.traj"),
                    ymd("2020-03-18"),
                    confirmed_cases %>% filter(population=="Iceland",
                                               confirmed_cumulativeInfections>0)) + 
       ggtitle("Iceland") +
       scale_x_date(limits = c(ymd("2020-02-15"),ymd("2020-03-18"))) +
       scale_y_log10(limits = c(1,1e4)),
       width=15, height=10, units="cm")

ggsave("figures/trajectories_iran.png",
       makeTrajPlot("../Results/BD_indep.iran.clock_8e-4.bu_36.5.1.traj",
                    ymd("2020-03-04"),
                    confirmed_cases %>% filter(population=="Iran",
                                               confirmed_cumulativeInfections>0)) + 
       ggtitle("Iran") +
       scale_x_date(limits = c(ymd("2019-12-01"),ymd("2020-03-04"))) +
       scale_y_log10(limits = c(1,1e5)),
       width=15, height=10, units="cm")

ggsave("figures/trajectories_italy.png",
       makeTrajPlot("results/BD_indep.italy.clock_8e-4.bu_36.5.1.traj",
                    ymd("2020-03-08"),
                    confirmed_cases %>% filter(population=="Italy",
                                               confirmed_cumulativeInfections>0)) + 
       ggtitle("Italy") +
       scale_x_date(limits = c(ymd("2020-01-01"),ymd("2020-03-08"))) +
       scale_y_log10(limits = c(1,1e4)),
       width=15, height=10, units="cm")

ggsave("figures/trajectories_spain.png",
       makeTrajPlot("results/BD_indep.spain.clock_8e-4.bu_36.5.1.traj",
                    ymd("2020-03-12"),
                    confirmed_cases %>% filter(population=="Spain",
                                               confirmed_cumulativeInfections>0)) + 
       ggtitle("Spain") +
       scale_x_date(limits = c(ymd("2019-12-15"),ymd("2020-03-12"))) +
       scale_y_log10(limits = c(1,2e3)),
       width=15, height=10, units="cm")

ggsave("figures/trajectories_washington.png",
       makeTrajPlot(c("results/BD_indep.washington1.clock_8e-4.bu_36.5.1.traj",
                      "results/BD_indep.washington2.clock_8e-4.bu_36.5.1.traj"),
                    ymd("2020-03-11"),
                    confirmed_cases %>% filter(population=="WA State",
                                               confirmed_cumulativeInfections>0)) + 
       ggtitle("WA State (USA)") +
       scale_x_date(limits = c(ymd("2020-01-15"),ymd("2020-03-11"))) +
       scale_y_log10(limits = c(1,1e4)),
       width=15, height=10, units="cm")

ggsave("figures/trajectories_wales.png",
       makeTrajPlot("results/BD_indep.welsh.clock_8e-4.bu_36.5.1.traj",
                    ymd("2020-03-16"),
                    confirmed_cases %>% filter(population=="Wales",
                                               confirmed_cumulativeInfections>0)) + 
       ggtitle("Wales") +
       scale_x_date(limits = c(ymd("2020-02-24"),ymd("2020-03-16"))) +
       scale_y_log10(limits = c(1,1e4)),
       width=15, height=10, units="cm")

ggsave("figures/trajectories_diamond_princess.png",
       makeTrajPlot("results/BD_indep.diamond_princess.clock_8e-4.bu_36.5.1.traj",
                    ymd("2020-02-25"),
                    confirmed_cases %>% filter(population=="Diamond Princess",
                                               confirmed_cumulativeInfections>0)) +
       ggtitle("Diamond Princess") +
       scale_x_date(limits = c(ymd("2020-01-20"),ymd("2020-02-25"))) +
       scale_y_log10(limits = c(1,1e4)),
       width=15, height=10, units="cm")
