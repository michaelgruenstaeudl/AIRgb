#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2019 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl")
#email = "m.gruenstaeudl@fu-berlin.de"
#version = "2019.11.20.1900"

########################################################################

library(svglite)
library(tcltk) # For dialog boxes
library(grid)  # For 'textGrob'
library(gridExtra) # For 'grid.arrange'

options(scipen=999) # Avoid E-notation in numbers

########################################################################

## Load Plot
PlastomeNumbersPerYear = readRDS(tk_choose.files(caption = "Select Rds file/plot '01_PlastomeNumbersPerYear'"))
PlastomesWithWithoutIRAnnos = readRDS(tk_choose.files(caption = "Select Rds file/plot '02a_PlastomesWithWithoutIRAnnos'"))
PlastomesWithUnequalReportedIRLengths = readRDS(tk_choose.files(caption = "Select Rds file/plot '02b_PlastomesWithUnequalReportedIRLengths'"))
IRLengthDiffDistribution = readRDS(tk_choose.files(caption = "Select Rds file/plot '02c_IRLengthDiffDistribution'"))

########################################################################

svglite(file="Figure1_IRsurvey.svg", width=21, height=22.275) #height=29.7)
grid.arrange(
    PlastomeNumbersPerYear,
    PlastomesWithWithoutIRAnnos,
    PlastomesWithUnequalReportedIRLengths,
    IRLengthDiffDistribution,
    nrow=2, ncol=2, 
    top=textGrob("Comparsion of IR in plastid genomes on NCBI Genbank", gp=gpar(fontsize=16, font=2))
)
dev.off()

########################################################################

### ONE PLOT PER LINE - LINE IS 1/4 OF PAGE

#svglite(file="VIZ_PlastomeNumbersPerYear.svg", width=21, height=7.425)
#PlastomeNumbersPerYear
#dev.off()

### TWO PLOTS PER LINE - LINE IS 1/4 OF PAGE

#svglite(file="VIZ_ComparsionOfReportedIRLengths.svg", width=21, height=7.425)
#grid.arrange(notEqual_finalPlot, LenDiff_finalPlot, ncol=2, 
#    top=textGrob("Comparsion of reported IR lengths", gp=gpar(fontsize=16, font=2))
#)
#dev.off()

########################################################################
