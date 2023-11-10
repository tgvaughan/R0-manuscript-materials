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

    logfiles <- dir("../Results", logfile_pattern, full.names=TRUE)

    if (length(logfiles) == 0)
        return(NULL)

    df <- load_and_combine_logs(logfiles, burnin=burnin)

    maybeDot <- ifelse("sampProp_australia.1" %in% names(df), ".", "")

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
                       sProp <- df[[paste0("sampProp_", outbreak_beast_codes[[outbreak]],
                                           maybeDot, "2")]]
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

                   sProp <- df[[paste0("sampProp_", outbreak_beast_codes[[outbreak]],
                                       maybeDot, "2")]]
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
                   sProp_splitTime <-  df[[paste0("sampChangeTimes_", outbreak_beast_codes[[outbreak]], maybeDot, "2")]]
                   sProp_splitTime_max <-  df[[paste0("sampChangeTimes_", outbreak_beast_codes[[outbreak]], maybeDot, "3")]]

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
                   sampProp_d <- df[[paste0("sampProp_", outbreak_beast_codes[[outbreak]], maybeDot, "2")]]
                   sampProp_post_d <- df[[paste0("sampProp_", outbreak_beast_codes[[outbreak]], maybeDot, "3")]]
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

    return(as_tibble(out_df))
}


results_indep <- process_results("BD_indep.clock_8e-4.bu_36.5.[1-5].log")

results_altPrior_indep <- process_results("BD_altPrior_indep.clock_8e-4.bu_36.5.[1-5].log",
                                          rR0prior=function(n) { runif(n, 0, 10) })

results_rateshift_indep <- process_results("BD_rateshift_indep.[1-5].log",
                                           rateshift="R0")

results_samprateshift_indep <- process_results("BD_samprateshift_indep.[1-5].log",
                                               rateshift="sampProp")

results_fixedrateshift_indep <- process_results("BD_fixedrateshift_indep.[1-5].log",
                                                rateshift="both")

results_linked <- process_results("BD_linked.clock_8e-4.bu_36.5.[1-5].log")

results_indep_noseq <- process_results("BD_indep_noseq.clock_8e-4.bu_36.5.[1-5].log")

results_linked_unprocessed <-
    load_and_combine_logs(dir("../Results",
                              "BD_linked.clock_8e-4.bu_36.5.[1-5].log",
                              full.names=TRUE))
                              

## Estimate R0 values from sample times

outbreak_Re_estimates <- NULL
outbreak_Re_stderrs <- NULL

df_Re <- NULL
for (outbreak in outbreaks) {

    datefile <- paste0("sequence_dates/",outbreak_beast_codes[[outbreak]],"_dates.txt")
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

df_Re_epiEstim <- read_csv("cluster_analysis_with_EpiEstim/R_estimates_clusters.csv",
                            col_types="cnnn",
                            quote = "") %>%
     group_by(cluster_id) %>%
     mutate(Outbreak=names(outbreak_beast_codes[outbreak_beast_codes==cluster_id]))


## Plot unique Re count derived from DPP analysis:

p <- ggplot(results_linked_unprocessed) +
    geom_histogram(aes(x=uniqueReCount, y=..density..),
                   breaks=seq(-.5,15.5,by=1),
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
sink("tables/Sequence_dates_and_counts.txt")
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
    theme_bw() + guides(fill="none") +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave(paste0("figures/Re_comparison_indep.pdf"), p,
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

## Rateshift model selection estimates
p <- ggplot(results_rateshift_indep) +
    geom_violin(aes(Outbreak, Re, fill=Population),
                draw_quantiles = c(0.025, 0.5, 0.975)) + 
    scale_y_continuous(limits=c(0,10)) +
    geom_hline(yintercept=1, linetype="dashed") +
    ylab(expression(R[0])) +
    theme_bw() + guides(fill=FALSE) +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave(paste0("figures/Re_comparison_rateshift_modelsel.png"), p,
       width=15, height=10, units="cm")

## Rateshift (conditioned on rate shift)
p <- ggplot(results_rateshift_indep %>% filter(Re_post != Re)) +
    geom_violin(aes(Outbreak, Re, fill=Population),
                draw_quantiles = c(0.025, 0.5, 0.975)) + 
    scale_y_continuous(limits=c(0,10)) +
    geom_hline(yintercept=1, linetype="dashed") +
    ylab(expression(R[0])) +
    theme_bw() + guides(fill=FALSE) +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave(paste0("figures/Re_comparison_rateshift_cond_for.png"), p,
       width=15, height=10, units="cm")

## Rateshift (conditioned on no shift)
p <- ggplot(results_rateshift_indep %>% filter(Re_post == Re)) +
    geom_violin(aes(Outbreak, Re, fill=Population),
                draw_quantiles = c(0.025, 0.5, 0.975)) + 
    scale_y_continuous(limits=c(0,10)) +
    geom_hline(yintercept=1, linetype="dashed") +
    ylab(expression(R[0])) +
    theme_bw() + guides(fill=FALSE) +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave(paste0("figures/Re_comparison_rateshift_cond_against.png"), p,
       width=15, height=10, units="cm")


# Rateshift post
p <- ggplot(results_rateshift_indep %>% filter(Re_post != Re)) +
    geom_violin(aes(Outbreak, Re_post, fill=Population),
                draw_quantiles = c(0.025, 0.5, 0.975)) + 
    scale_y_continuous(limits=c(0,10)) +
    geom_hline(yintercept=1, linetype="dashed") +
    ylab(expression(R[0])) +
    theme_bw() + guides(fill=FALSE) +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave(paste0("figures/Re_comparison_rateshift_post.png"), p,
       width=15, height=10, units="cm")

# Rateshift times (absolute)
p <- ggplot(results_rateshift_indep %>% filter(Re_post != Re, Outbreak != "Prior")) +
    geom_violin(aes(Outbreak, Re_splitTime, fill=Population),
                draw_quantiles = c(0.025, 0.5, 0.975)) + 
    scale_y_continuous(limits=c(0,1)) +
    ## ylab(expression(R[0])) +
    theme_bw() + guides(fill=FALSE) +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave(paste0("figures/Re_comparison_rateshift_times_abs.png"), p,
       width=15, height=10, units="cm")

# Rateshift times (relative)
p <- ggplot(results_rateshift_indep %>% filter(Re_post != Re)) +
    geom_violin(aes(Outbreak, Re_splitTimeRel, fill=Population),
                draw_quantiles = c(0.025, 0.5, 0.975)) + 
    scale_y_continuous(limits=c(0,1)) +
    ## ylab(expression(R[0])) +
    theme_bw() + guides(fill=FALSE) +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave(paste0("figures/Re_comparison_rateshift_times_rel.png"), p,
       width=15, height=10, units="cm")

# Scatter plots
p <- ggplot(results_rateshift_indep) +
    geom_point(aes(Re, Re_post, color=Population), alpha=0.05) +
    facet_wrap(vars(Outbreak), scales="free") +
    xlab(expression(R[1])) +
    ylab(expression(R[2])) +
    theme_bw() + guides(fill=FALSE) + guides(color=FALSE)
    ## theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave(paste0("figures/Re_pre_post.png"), p,
       width=20, height=20, units="cm")

# Sample rate shift model selection estimates

p <- ggplot(results_samprateshift_indep) +
    geom_violin(aes(Outbreak, Re, fill=Population),
                draw_quantiles = c(0.025, 0.5, 0.975)) + 
    scale_y_continuous(limits=c(0,10)) +
    geom_hline(yintercept=1, linetype="dashed") +
    ylab(expression(R[0])) +
    theme_bw() + guides(fill=FALSE) +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave(paste0("figures/Re_comparison_samprateshift_modelsel.png"), p,
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

sink("tables/Re_and_samplingProp_indep.txt")
results_indep %>%
    group_by(Outbreak) %>%
    summarize(Re_med = median(Re), Re_low = quantile(Re, 0.025), Re_high = quantile(Re, 0.975),
              s_med = median(samplingProp), s_low = quantile(samplingProp, 0.025),
              s_high = quantile(samplingProp, 0.975))
sink()



## Supplemental plot showing marginal Re posteriors
## alongside basic regression estimates and epiestim results

p <- ggplot() +
    geom_violin(aes("BDSKY (with sequences)", Re),
                results_indep %>% filter (Outbreak != "Prior"),
                draw_quantiles = c(0.025, 0.5, 0.975),
                fill="purple") + 
    geom_violin(aes("BDSKY (without sequences)", Re),
                results_indep_noseq %>% filter(Outbreak != "Prior"),
                draw_quantiles = c(0.025, 0.5, 0.975),
                fill="orange") + 
    ## scale_y_continuous(limits=c(0,10)) +
    geom_pointrange(aes("Regression", Re_est,
                        ymin=Re_est-2*Re_est_err,
                        ymax=Re_est+2*Re_est_err),
                    df_Re %>% filter(Outbreak != "Diamond Princess"),
                    col="blue") +
    geom_pointrange(aes("EpiEstim", R_mean,
                        ymin=R_lowQt,
                        ymax=R_highQt),
                    df_Re_epiEstim %>% filter(Outbreak != "Diamond Princess"),
                    col="red") +
    geom_hline(yintercept=1, linetype="dashed") +
    coord_cartesian(ylim=c(0,10)) +
    ylab(expression(R[0])) +
    ## scale_colour_manual(name=element_blank(),
    ##                     values=c("Regression"="blue",
    ##                              "EpiEstim"="red")) +
    ## scale_fill_manual(name=element_blank(),
    ##                   values=c("BDSKY (with sequences)"="purple",
    ##                            "BDSKY (no sequences)"="orange")) +
    facet_wrap(facets=vars(Outbreak), ncol=5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave("figures/Re_comparison_indep_noseq_regression.png", p,
       width=20, height=20, units="cm")



## Case count marginal posteriors (from EpiInf)

source("~/code/beast_and_friends/EpiInf/scripts/trajDataTools.R")

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

scientific_10 <- function(x) {
    parse(text=gsub("1e\\+*", "10^", scales::label_scientific()(x)))
}

p <- ggplot(casecountdf) +
    geom_violin(aes(population, cumulativeInfections, fill=population),
                draw_quantiles=c(0.025,0.5,0.975)) +
    geom_point(aes(population, confirmed_cumulativeInfections),
               final_confirmed_cases_offset, shape=5, size=2.5, stroke=2, color="white") +
    geom_point(aes(population, confirmed_cumulativeInfections),
               final_confirmed_cases_offset, shape=5, size=3.0, stroke=1, color="red") +
    ## scale_y_log10() +
    scale_y_log10(label=scientific_10) +
    ylab("Total cases") +
    xlab("Population") +
    guides(fill="none") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave("figures/caseCounts_comparison_indep_offset.pdf", p,
       width=15, height=10, units="cm")


## Write a table containing case count estimates

sink("tables/case_count_indep.txt")
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

sink("tables/case_count_indep_offset.txt")
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


### Estimating total case counts and cryptic genetic diversity from infection fatality rate

# Load death data:

non_us_raw_data <- read_csv("Johns_Hopkins/repository/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv") %>%
    pivot_longer(cols=-(1:4), names_to="date") %>%
    rename(country="Country/Region",state="Province/State")

non_us_data <- non_us_raw_data %>%
    filter(country %in% c("Australia",
                          "China",
                          "France",
                          "Iceland",
                          "Iran",
                          "Italy",
                          "Netherlands",
                          "Spain",
                          "Diamond Princess"
                          )) %>%
    mutate(date=mdy(date)) %>%
    mutate(country=ifelse(country=="Netherlands","The Netherlands",country)) %>%
    group_by(country,date) %>%
    summarize(value=sum(value)) %>%
    rename(Population=country)

us_raw_data <- read_csv("Johns_Hopkins/repository/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_us.csv") %>%
    pivot_longer(cols=-(1:12), names_to="date") %>%
    rename(state=Province_State)

us_data <- us_raw_data %>%
    filter(state == "Washington") %>%
    mutate(date=mdy(date), state="WA State") %>%
    group_by(state, date) %>%
    summarize(value=sum(value)) %>%
    rename(Population=state)

death_data <- bind_rows(non_us_data, us_data)

final_IFR_cases_offset <- outbreak_dates %>%
    left_join(tibble(Population=outbreak_population, Outbreak=names(outbreak_population))) %>%
    group_by(Population) %>%
    summarize(date=max(Date)) %>%
    inner_join(death_data %>% mutate(date=date-18)) %>% #Offset + 8 to account for test->death lag
    transmute(population=Population,date,IFR_cumulativeInfections=value/0.0064,
              IFR_cumulativeInfections_high=value/0.0038,
              IFR_cumulativeInfections_low=value/0.0098)

# Case count plot

p <- ggplot(casecountdf) +
    geom_violin(aes(population, cumulativeInfections, fill=population),
                draw_quantiles=c(0.025,0.5,0.975)) +
    geom_point(aes(population, confirmed_cumulativeInfections),
               final_confirmed_cases_offset, shape=5, size=2.5, stroke=2, color="white") +
    geom_point(aes(population, confirmed_cumulativeInfections),
               final_confirmed_cases_offset, shape=5, size=3.0, stroke=1, color="red") +
    ## geom_point(aes(population, IFR_cumulativeInfections),
    ##            final_IFR_cases_offset, shape=1, size=2.5, stroke=2, color="white") +
    ## geom_point(aes(population, IFR_cumulativeInfections),
    ##            final_IFR_cases_offset, shape=1, size=3.0, stroke=1, color="purple") +
    geom_errorbar(aes(population, ymin=IFR_cumulativeInfections_low,
                        ymax=IFR_cumulativeInfections_high),
                    final_IFR_cases_offset, linewidth=2, width=0.25, color="white") +
    geom_point(aes(population, IFR_cumulativeInfections),
                    final_IFR_cases_offset, stroke=2, color="white") +
    geom_errorbar(aes(population, ymin=IFR_cumulativeInfections_low,
                        ymax=IFR_cumulativeInfections_high),
                    final_IFR_cases_offset, linewidth=1, width=0.2, color="purple") +
    geom_point(aes(population, IFR_cumulativeInfections),
                    final_IFR_cases_offset, stroke=1, color="purple") +
    scale_y_log10() +
    ylab("Total cases") +
    xlab("Population") +
    guides(fill="none") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave("figures/caseCountsIFR_comparison_indep_offset.png", p,
       width=15, height=10, units="cm")

# Cryptic genetic diversity plot

sampleddiversitydf <- casecountdf %>%
    left_join(final_IFR_cases_offset, by="population") %>%
    mutate(sampled_diversity=cumulativeInfections/IFR_cumulativeInfections)

sampleddiversitydf$population <- factor(sampleddiversitydf$population, levels=populations)

p <- ggplot(sampleddiversitydf %>% filter(population!="Wales")) +
    geom_violin(aes(population, sampled_diversity, fill=population),
                draw_quantiles=c(0.025,0.5,0.975)) +
    scale_y_log10() +
    ylab("Sampled diversity ratio") +
    xlab("Population") +
    guides(fill="none") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
ggsave("figures/sampled_diversity_indep_offset.png", p,
       width=15, height=10, units="cm")


### Empirical sampling fraction plots

sequence_dates <- NULL
for (outbreak in outbreaks) {
    datefile <- paste0("sequence_dates/",outbreak_beast_codes[[outbreak]],"_dates.txt")
    sequence_dates <- bind_rows(sequence_dates,
                                read_tsv(datefile, col_names=c("Date")) %>%
                                mutate(Outbreak=outbreak,
                                       Population=outbreak_population[[outbreak]]))
}
sequence_dates <- sequence_dates %>%
    arrange(by_group=Date) %>%
    group_by(Population) %>%
    mutate(Count=row_number())

# Calculate normalized cumulative sample counts
total_sequence_counts <- sequence_dates %>%
    group_by(Population) %>%
    summarize(totalCount=max(Count))
sequence_dates_ecdf <- 
    sequence_dates %>%
    left_join(total_sequence_counts) %>%
    filter(Population != "China") %>%
    mutate(normCount=Count/totalCount)

# Calculate normalized cumulative confirmed cases
final_sample_dates <- sequence_dates %>%
    group_by(Population) %>%
    summarize(finalSampleDate=max(Date), initialSampleDate=min(Date))
total_confirmed_cases <- confirmed_cases %>%
    transmute(Population=population, Date=date, Count=confirmed_cumulativeInfections) %>%
    left_join(final_sample_dates) %>%
    filter(Date<=finalSampleDate) %>%
    filter(Date>=initialSampleDate) %>%
    group_by(Population) %>%
    summarize(initialCount=min(Count), finalCount=max(Count))
confirmed_cases_ecdf <- confirmed_cases %>%
     transmute(Population=population, Date=date, Count=confirmed_cumulativeInfections) %>%
    left_join(total_confirmed_cases) %>%
    left_join(final_sample_dates) %>%
    filter(Date<=finalSampleDate) %>%
    filter(Date>=initialSampleDate) %>%
    filter(Population != "China") %>%
    mutate(normCount=(Count-initialCount)/(finalCount-initialCount))
   
# Combined plot of empirical CDFs
p <- ggplot(confirmed_cases_ecdf, aes(Date, normCount, col=Population)) +
    geom_line() +
    facet_wrap(facets=vars(Population),ncol=4, scales="free") +
    geom_line(data=sequence_dates_ecdf, linetype="dashed") +
    guides(col="none") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=-45, vjust=0.5, hjust=0)) +
    xlab(NULL) + ylab("Cumulative proportion")
 ggsave("figures/cumulative_sample_count_comparison.png", p,
       width=15, height=15, units="cm")
   

### Trajectory plots

source("~/code/beast_and_friends/EpiInf/scripts/trajDataTools.R")

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
       makeTrajPlot("../Results/BD_indep.italy.clock_8e-4.bu_36.5.1.traj",
                    ymd("2020-03-08"),
                    confirmed_cases %>% filter(population=="Italy",
                                               confirmed_cumulativeInfections>0)) + 
       ggtitle("Italy") +
       scale_x_date(limits = c(ymd("2020-01-01"),ymd("2020-03-08"))) +
       scale_y_log10(limits = c(1,1e4)),
       width=15, height=10, units="cm")

ggsave("figures/trajectories_spain.png",
       makeTrajPlot("../Results/BD_indep.spain.clock_8e-4.bu_36.5.1.traj",
                    ymd("2020-03-12"),
                    confirmed_cases %>% filter(population=="Spain",
                                               confirmed_cumulativeInfections>0)) + 
       ggtitle("Spain") +
       scale_x_date(limits = c(ymd("2019-12-15"),ymd("2020-03-12"))) +
       scale_y_log10(limits = c(1,2e3)),
       width=15, height=10, units="cm")

ggsave("figures/trajectories_washington.png",
       makeTrajPlot(c("../Results/BD_indep.washington1.clock_8e-4.bu_36.5.1.traj",
                      "../Results/BD_indep.washington2.clock_8e-4.bu_36.5.1.traj"),
                    ymd("2020-03-11"),
                    confirmed_cases %>% filter(population=="WA State",
                                               confirmed_cumulativeInfections>0)) + 
       ggtitle("WA State (USA)") +
       scale_x_date(limits = c(ymd("2020-01-15"),ymd("2020-03-11"))) +
       scale_y_log10(limits = c(1,1e4)),
       width=15, height=10, units="cm")

ggsave("figures/trajectories_wales.png",
       makeTrajPlot("../Results/BD_indep.welsh.clock_8e-4.bu_36.5.1.traj",
                    ymd("2020-03-16"),
                    confirmed_cases %>% filter(population=="Wales",
                                               confirmed_cumulativeInfections>0)) + 
       ggtitle("Wales") +
       scale_x_date(limits = c(ymd("2020-02-24"),ymd("2020-03-16"))) +
       scale_y_log10(limits = c(1,1e4)),
       width=15, height=10, units="cm")

ggsave("figures/trajectories_diamond_princess.png",
       makeTrajPlot("../Results/BD_indep.diamond_princess.clock_8e-4.bu_36.5.1.traj",
                    ymd("2020-02-25"),
                    confirmed_cases %>% filter(population=="Diamond Princess",
                                               confirmed_cumulativeInfections>0)) +
       ggtitle("Diamond Princess") +
       scale_x_date(limits = c(ymd("2020-01-20"),ymd("2020-02-25"))) +
       scale_y_log10(limits = c(1,1e4)),
       width=15, height=10, units="cm")


### Topological ESS plots

library(rwty)

ess_plots <- NULL
for (outbreak in outbreaks) {
    cat(paste("Processing outbreak", outbreak, "...\n"))

    pathPattern <- paste0("BD_indep.",
                          outbreak_beast_codes[outbreak],
                          ".clock_8e-4.bu_36.5.[1,2].trees")

    files <- dir(path="../Results",
                 pattern=pathPattern,
                 full.names=TRUE)

    chains <- lapply(files,
                     function (f) {
                         print(paste("Loading trees from",f));
                         return(load.trees(f, format="beast"))
                     })
    
    ess_plots[[outbreak]] <- makeplot.topology(chains, burnin=2000)$density.plot
}

for (outbreak in outbreaks) {
    ggsave(paste0("figures/topology_distr_",outbreak_beast_codes[outbreak],".png"),
           ess_plots[[outbreak]] + ggtitle(outbreak) + xlab("Topological Distance"),
           width=8, height=8, units="cm")
}


### Sequence identity plots

library(seqinr)

aliFiles <- dir("sequences", pattern="*.masked", full.names=T)

alignments <- list()
for (i in 1:length(aliFiles)) 
    alignments[[i]] <- read.alignment(aliFiles[i], "fasta")

get_pairwise_identity <- function(seq1, seq2) {
    L <- nchar(seq1)

    nIdent <- 0
    for (i in 1:L) {
        c1 <- substr(seq1,i,i)
        c2 <- substr(seq2,i,i)
        if (c1!="n" && c2!="n" && c1!=c2)
            return(0)
    }

    return(1)
}

get_average_pairwise_identity_within <- function(ali) {
    pairs <- 0
    totalIdent <- 0
    for (i in 1:(ali$nb-1)) {
        for (j in (i+1):ali$nb) {
            totalIdent <- totalIdent + get_pairwise_identity(ali$seq[i], ali$seq[j])
            pairs <- pairs + 1
        }
    }

    return(totalIdent/pairs)
}

get_average_pairwise_identity_between <- function(ali1, ali2) {
    pairs <- 0
    totalIdent <- 0
    for (i in 1:ali1$nb) {
        for (j in 1:ali2$nb) {
            totalIdent <- totalIdent + get_pairwise_identity(ali1$seq[i], ali2$seq[j])
            pairs <- pairs + 1
        }
    }

    return(totalIdent/pairs)
}


identities <- matrix(NA, nrow=length(aliFiles), ncol=length(aliFiles))

for (i in 1:length(alignments)) {
    for (j in i:length(alignments)) {
        if (i == j) {
            cat(paste("Pairwise identity within:", aliFiles[[i]],"\n"))
            identities[i,i] <- get_average_pairwise_identity_within(alignments[[i]])
        } else {
            cat(paste("Pairwise identity between:", aliFiles[[i]],"and",aliFiles[[j]],"\n"))
            identities[i,j] <- get_average_pairwise_identity_between(alignments[[i]],
                                                                     alignments[[j]])
            identities[j,i] <- identities[i,j]
        }
    }
}

name_map <- as_tibble(outbreak_beast_codes, rownames="outbreak_name") %>% rename(beast_name=value)
aliNames <- tibble(file_name=aliFiles) %>%
    mutate(beast_name=str_replace(file_name,"^sequences/", "")) %>%
    mutate(beast_name=str_replace(beast_name,".masked$","")) %>%
    left_join(name_map, by=join_by(beast_name))

colnames(identities) <- aliNames$outbreak_name
rownames(identities) <- aliNames$outbreak_name

p <- as_tibble(identities, rownames="outbreak1") %>%
    pivot_longer(cols=-1, names_to="outbreak2", values_to="identity") %>%
    ggplot(aes(x=outbreak1, y=outbreak2, fill=identity)) + geom_raster() +
    geom_text(aes(label=round(identity, digits=3)), col="white") +
    ## scale_fill_distiller(palette="Purples") +
    theme(axis.title=element_blank(),
          axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5)) +
    guides(fill="none", color="none")
ggsave("figures/sequence_identity.png", p, width=23, height=18, units="cm")


### Shuffled sequence analysis result plots

## Load 

results_indep_shuffled <- NULL
for (i in 1:10) {
    results_indep_shuffled <-
        bind_rows(results_indep_shuffled,
                  process_results(paste0("BD_indep_shuffled.clock_8e-4.bu_36.5.",i,".log")) %>%
                  mutate(shuffle=i))
}

combined_data <- results_indep_shuffled %>%
    mutate(analysis=paste0("Shuffled ",shuffle), type="Shuffled") %>%
    select(-shuffle) %>%
    bind_rows(results_indep %>% mutate(analysis="Unshuffled", type="Unshuffled")) %>%
    bind_rows(results_indep_noseq %>% mutate(analysis="Dates only", type="Dates only")) %>%
    mutate(analysis = fct_relevel(analysis, c("Dates only", "Unshuffled",
                                              paste("Shuffled", 1:10)))) %>%
    mutate(type = fct_relevel(type, c("Dates only", "Unshuffled", "Shuffled")))

## Indep
p <- ggplot(combined_data %>% filter(Population != "Prior")) +
    geom_violin(aes(analysis, Re, fill=type),
                draw_quantiles = c(0.025, 0.5, 0.975)) +
    geom_hline(yintercept=1, linetype="dashed") +
    facet_wrap(facets=vars(Outbreak), ncol=3, scales="free_y") +
    ylab(expression(R[0])) +
    theme_bw() +
    guides(fill=guide_legend(title="Analysis type:")) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "bottom")
ggsave(paste0("figures/Re_comparison_indep_shuffled.pdf"), p,
       width=20, height=20, units="cm")

