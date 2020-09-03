#!/usr/bin/R
#author = "Tilman Mehl"
#copyright = "Copyright (C) 2019-2020 Tilman Mehl, Michael Gruenstaeudl"
#contributors = c("Tilman Mehl", "Michael Gruenstaeudl")
#email = "m.gruenstaeudl@fu-berlin.de"
#version = "2020.09.03.1230"

########################################################################

library(ggplot2)
library(tcltk) # For dialog boxes
library(tools) # For function 'file_path_sans_ext'
library(dplyr) # For function '%>%'
library(sqldf)

########################################################################

# GETTING SCRIPT NAME
args = commandArgs(TRUE)
this_script = sub(".*=", "", commandArgs()[4])
script_name = file_path_sans_ext(basename(this_script))

########################################################################

## Load Plastome Availability Table (.csv-format)
inFn1 = tk_choose.files(caption = "Select the plastome availability table (.tsv-format)")
plastData = read.csv(inFn1, sep = "\t")
## Load Reported IR Information Table (.csv-format)
inFn2 = tk_choose.files(caption = "Select the reported IR information table (.tsv-format)")
reportedData = read.csv(inFn2, sep = "\t")

# Only account for accessions that appear in the reported data
plastData <- sqldf("SELECT plastData.* FROM plastData, reportedData WHERE plastData.ACCESSION == reportedData.ACCESSION")

## Replace "n.a." Strings in table with proper <NA> values and cast the columns to numeric
reportedData <- mutate_all(reportedData, funs(replace(.,. == "n.a.", NA)))
reportedData$IRa_REPORTED_START <- as.integer(as.character(reportedData$IRa_REPORTED_START))
reportedData$IRb_REPORTED_START <- as.integer(as.character(reportedData$IRb_REPORTED_START))
reportedData$IRa_REPORTED_END <- as.integer(as.character(reportedData$IRa_REPORTED_END))
reportedData$IRb_REPORTED_END <- as.integer(as.character(reportedData$IRb_REPORTED_END))
reportedData$IRa_REPORTED_LENGTH <- as.integer(as.character(reportedData$IRa_REPORTED_LENGTH))
reportedData$IRb_REPORTED_LENGTH <- as.integer(as.character(reportedData$IRb_REPORTED_LENGTH))

qualityData = data.frame(reportedData$ACCESSION)
qualityData$IS_CG <- reportedData$IRa_REPORTED == "yes" & reportedData$IRb_REPORTED == "yes" & (reportedData$IRa_REPORTED_LENGTH == reportedData$IRb_REPORTED_LENGTH)
qualityData$VERSION <- as.numeric(plastData$VERSION)
colnames(qualityData) <- c("ACCESSION","IS_CG", "VERSION")

qualityData$count <- 1
plotData <- aggregate(qualityData$count, by=list(qualityData$VERSION, qualityData$IS_CG), FUN=sum)
colnames(plotData) <- c("VERSION", "IS_CG", "TOTAL")
acc_per_version <- aggregate(qualityData$count, by=list(qualityData$VERSION), FUN=sum)
colnames(acc_per_version) <- c("VERSION", "TOTAL")

plotData$PERCENTAGE <- as.numeric(sqldf("SELECT (plotData.TOTAL/acc_per_version.TOTAL) FROM plotData, acc_per_version WHERE plotData.VERSION = acc_per_version.VERSION")[,1])
ggplot(data=plotData, mapping = aes(x=VERSION, y=PERCENTAGE)) + geom_col(aes(fill=IS_CG))
