## Script used for estimating R0 from sequence collection dates only using EpiEstim.

library("ggplot2")
library("lubridate")
library("readr")
library("EpiEstim")
library("gridExtra")
library("reshape2")
library(RColorBrewer)
library(reshape)
library(ggpubr)
library(wesanderson)
library(plyr)


### Help function for plotting
gg_color <- function(n) {
  
  hues = seq(15, 375, length = n + 1)
  
  hcl(h = hues, l = 65, c = 100)[1:n]
  
}

## Get incidence data from cumulative counts data and restructure dataframe
meltCumulativeData <- function(rawData, dataType) {
  
  cumulData <- rawData
  cumulData$Date <- ymd(cumulData$Date)
 
  incidenceData <- as.data.frame(apply(cumulData[,-1, drop=F], 2, function(x) {  incidence <- x - c(0, x[1:(length(x)-1)]) ; return(incidence)}))
  incidenceData <- cbind(Date=cumulData$Date, incidenceData)
  incidenceData <- melt(incidenceData, id.vars="Date")
  colnames(incidenceData) <- c("date", "region", "value")
  incidenceData$data_type <- dataType
  incidenceData$variable <- "incidence"
  incidenceData$estimate_type <- NA
  
  cumulData <- melt(cumulData, id.vars="Date")
  colnames(cumulData) <- c("date", "region", "value")
  cumulData$data_type <- dataType
  cumulData$variable <- "cumul"
  cumulData$estimate_type <- NA
  
  return(rbind(cumulData, incidenceData))
}


getClusterData <- function(pathToCSV="../sequences/cumulative_incidence_from_clusters_to_last_sample.txt") {
  cumData <- read.table(pathToCSV, header=T)
  colnames(cumData)[1] <- "Date"
  return(meltCumulativeData(cumData, "seq_samples"))
}

### Apply EpiEstim R estimation method to 'incidenceData' timeseries with 'dates' the dates associated
##
## 'estimateOffsetting' is the number of days the estimates are to be shifted towards the past (to account for delay between infection and testing/hospitalization/death..)
## 'method' takes value either 'Cori' or  'WallingaTeunis'. 'Cori' is the classic EpiEstim R(t) method, 'WallingaTeunis' is the method by Wallinga and Teunis (also implemented in EpiEstim)
## 'minimumCumul' is the minimum cumulative count the incidence data needs to reach before the first Re estimate is attempted (if too low, EpiEstim can crash)
## 'mean_si' and 'std_si' are the mean and SD of the serial interval distribution used by EpiEstim
estimate_constant_Re <- function(dates, incidenceData, estimateOffsetting = 0,  method="Cori", minimumCumul = 1,  mean_si = 4.8, std_si  =2.3) {


  ## First, remove missing data at beginning of series
  while(length(incidenceData) > 0 & is.na(incidenceData[1])) {
    incidenceData <- incidenceData[-1]
    dates <- dates[-1]
    if(length(incidenceData) == 0) {
      return(data.frame(date=c(), variable=c(), value=c(), estimate_type=c()))
    }
  }
  
  ## Then, remove missing data at the end of the series
  while(length(incidenceData) > 0 & is.na(incidenceData[length(incidenceData)])) {
    incidenceData <- incidenceData[-length(incidenceData)]
    dates <- dates[-length(dates)]
    if(length(incidenceData) == 0) {
      return(data.frame(date=c(), variable=c(), value=c(), estimate_type=c()))
    }
  }
  
  ## Replace missing data in rest of series by zeroes (required for using EpiEstim)
  incidenceData[is.na(incidenceData)] <- 0
  
  offset <- 1
  cumulativeIncidence <- 0
  while(cumulativeIncidence < minimumCumul) {
    if(offset > length(incidenceData)) {
      return(data.frame(date=c(), variable=c(), value=c(), estimate_type=c()))
    }
    cumulativeIncidence <- cumulativeIncidence + incidenceData[offset]
    offset <- offset + 1
  }
  
  ## offset needs to be at least two for EpiEstim
  offset <- max(2, offset)
  
  ## generate start and end bounds for Re estimates
  t_start <- offset
  t_end <- length(incidenceData)
  
  if(length(incidenceData) < offset) { ## no valid data point, return empty estimate
    return(data.frame(date=c(), variable=c(), value=c(), estimate_type=c()))
  }
  
  if(method == "Cori") {
    
    R_instantaneous <- estimate_R(incidenceData, 
                                  method="parametric_si", 
                                  config = make_config(list(
                                    mean_si = mean_si, std_si = std_si,
                                    t_start = t_start,
                                    t_end = t_end)))
    
  } else if(method == "WallingaTeunis") {
    
    R_instantaneous <- wallinga_teunis(incidenceData,
                                       method="parametric_si",
                                       config = list(
                                         mean_si = mean_si, std_si = std_si,
                                         t_start = t_start,
                                         t_end = t_end,
                                         n_sim = 100))
  } else {
    print("Unknown estimation method")
    return(data.frame(date=c(), variable=c(), value=c(), estimate_type=c()))
  }
  
  outputDates <- dates[t_start:t_end]
  
  ## estimateOffsetting accounts for delay between infection and recorded event (testing, hospitalization, death...)
  outputDates <- outputDates - estimateOffsetting
  
  R_mean <- R_instantaneous$R$`Mean(R)`
  R_highHPD <- R_instantaneous$R$`Quantile.0.975(R)`
  R_lowHPD <- R_instantaneous$R$`Quantile.0.025(R)`
  
  result <- data.frame(date=outputDates,
                       R_mean=R_mean, 
                       R_highHPD=R_highHPD,
                       R_lowHPD=R_lowHPD)
  
  result <- melt(result, id.vars="date")
  colnames(result) <- c("date", "variable", "value")
  result$estimate_type <- method
  
  return(result)
}

## Perform R(t) estimations with EpiEstim on each 'region' of the data, with each 'method' and on each 'data_type'
## 'region' is the geographical region
## 'data_type' is "seq_samples" in the case of the cluster data
doAllClusterReEstimations <- function(data,methods=c("Cori")) {
  
  fullResults <- data.frame(date=c(), region=c(), value=c(),data_type=c(), variable=c(), estimate_type=c())
  
  # do the Re estimations for each region, each type of data, and each Re estimation method
  for(region_i in unique(data$region)) {
    subset_data_region <- subset(data, region == region_i)
    for(data_type_i in unique(subset_data_region$data_type)) {
      subset_data <- subset(subset_data_region, data_type == data_type_i)
      for(method_i in methods) {
        incidence_data <- subset_data$value[subset_data$variable == "incidence"]
        dates <- subset_data$date[subset_data$variable == "incidence"]
        if(method_i == "Cori") {
          result <- estimate_constant_Re(dates=dates, incidenceData=incidence_data, method=method_i)
        } else if(method_i == "WallingaTeunis") {
          result <- estimate_constant_Re(dates=dates, incidenceData=incidence_data, method=method_i)
        }
        
        if(nrow(result) > 0) {
          result$region <- region_i
          result$data_type <- data_type_i
          ## need to reorder columns in 'results' dataframe to do the same as in data
          result <- result[,c(1,5,3,6,2,4)]
          fullResults <- rbind(fullResults, result)
        }
      }
    }
  }
  return(rbind(data, fullResults))
}

saveAndPlotEstimates <- function(meltData, dataType=c("seq_samples"), estimateType="Cori", nameFile="R_estimates_clusters", pathToCSV=".", pathToPlot=".") {
 
   # meltData <- subset(meltData, meltData$data_type == dataType)
  csvData <- subset(meltData, data_type == dataType & !is.na(estimate_type))
  csvData <- subset(csvData, estimate_type ==estimateType )
  
  csvData <-  csvData %>%
    select(region, value, variable) %>%
    group_by(region, variable) %>%
    summarise(R = mean(value))
  
  saveData <- csvData %>% spread(key=variable, value=R)
  saveData <- saveData[,c(1,4,2,3)]
  colnames(saveData) <- c("cluster_id", "R_mean", "R_lowQt", "R_highQt")
  
  
  write.csv(saveData,file= paste0(pathToCSV, nameFile, ".csv"), quote=F, row.names=F)
  
  ggplot(data=saveData, aes(x=cluster_id)) +
    geom_hline(yintercept=1,linetype="dashed", size=1.1) +
    geom_pointrange(aes(y = R_mean,
                          ymin=R_lowQt,
                          ymax=R_highQt),
      stat = "identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text.y= element_text(size=12),
          axis.title.y =  element_text(size=17)) +
    xlab("") + 
    ylab("Reproductive number") + 
    coord_cartesian(ylim=c(0,9))
  ggsave(paste0(pathToPlot, nameFile, ".png"), width = 20, height = 15, units = "cm")
  
}

#################################
#### Start of the script ########
#################################

cluster_data <- getClusterData()
estimates <-  doAllClusterReEstimations(cluster_data)
saveAndPlotEstimates(estimates)
