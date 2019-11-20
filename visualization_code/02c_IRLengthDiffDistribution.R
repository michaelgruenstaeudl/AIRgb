#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2019 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl")
#email = "m.gruenstaeudl@fu-berlin.de"
#version = "2019.11.20.1700"

########################################################################

library(ggplot2)
library(tcltk) # For dialog boxes
library(tools) # For function 'file_path_sans_ext'
library(dplyr) # For function '%>%'

########################################################################

# GETTING SCRIPT NAME
args = commandArgs(TRUE)
this_script = sub(".*=", "", commandArgs()[4])
script_name = file_path_sans_ext(basename(this_script))

########################################################################

#GLOBAL VARIABLES
options(scipen=999) # Avoid E-notation in numbers
start_year = 2010

########################################################################

## Load Plastome Availability Table (.csv-format)
inFn1 = tk_choose.files(caption = "Select the plastome availability table (.tsv-format)")
inData1 = read.csv(inFn1, sep = "\t")

## Load Plastome Availability Table (.csv-format)
inFn2 = tk_choose.files(caption = "Select the reported IR stats table (.tsv-format)")
#out_fn = paste(file_path_sans_ext(inFn2), "_", sep='')
out_fn = dirname(inFn2)
inData2 = read.csv(inFn2, sep = "\t")

## Combine the input files
combinedDF = merge(inData1, inData2, by="ACCESSION")

########################################################################

## Ensure that numeric values are recognized as numeric (i.e., transform from factor to numeric via character)
## See: https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-integer-numeric-without-loss-of-information

combinedDF = transform(combinedDF, 
IRa_REPORTED_START = as.integer(as.character(IRa_REPORTED_START)),
IRb_REPORTED_START = as.integer(as.character(IRb_REPORTED_START)),
IRa_REPORTED_END = as.integer(as.character(IRa_REPORTED_END)),
IRb_REPORTED_END = as.integer(as.character(IRb_REPORTED_END)),
IRa_REPORTED_LENGTH = as.integer(as.character(IRa_REPORTED_LENGTH)),
IRb_REPORTED_LENGTH = as.integer(as.character(IRb_REPORTED_LENGTH))
)

########################################################################

## Getting only rows that have reported both IRa and IRb
relevantDF = combinedDF[which(combinedDF[,"IRa_REPORTED"]=="yes" & combinedDF[,"IRb_REPORTED"]=="yes"),]

## Seperate rows by criterion
notEqualDF = relevantDF[!(relevantDF[,"IRa_REPORTED_LENGTH"]==relevantDF[,"IRb_REPORTED_LENGTH"]),]

########################################################################

## PLOTTING IR LENGTH DIFFERENCES AS HISTOGRAM ##

## Add a column that displays the number of nucleotide differences between the IR lengths
notEqualDF = notEqualDF %>% mutate(IR_LEN_DIFF=abs(notEqualDF$IRa_REPORTED_LENGTH-notEqualDF$IRb_REPORTED_LENGTH))

####################################

base_plot = ggplot(data=notEqualDF, aes(x=IR_LEN_DIFF)) + 
    geom_histogram(fill="grey50", alpha=0.5) #+ 
    #geom_density(aes(y=..density..*n), color="grey0")

myPlot = base_plot + 
    xlab("\nLength Difference (in bp)") + 
    ylab("Total Number of Occurrences\n") + 
    ggtitle("Distribution of the differences in length\nbetween the IRa and the IRb",
            subtitle="Note: x-axis is set to logarithmic scale"
    ) + 
    scale_x_log10(
        #limits=c(1,200000),
        breaks=c(1,2,3,4,5,10,100,1000,10000,100000),
        minor_breaks=NULL,
        expand=expand_scale(0),
        labels=c(1,2,3,4,5,10,100,1000,10000,100000)
    ) +
    #scale_fill_manual(values=c("grey50", "grey0"), name="Plastid genome\navailable", labels=c("Yes", "No")) +
    #theme_bw() + 
    theme_minimal() + 
    theme(plot.title = element_text(size=20),
          plot.subtitle = element_text(size=16, face="italic"),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          plot.margin=unit(c(0.5,1.0,0.1,0.1),"cm"),  # Note: margin(t=0, r=0, b=0, l=0, unit="pt")
          legend.key.width=unit(1,"cm"))

####################################

assign(script_name, myPlot)
saveRDS(eval(as.name(script_name)), file=paste(out_fn, '/', script_name, ".Rds", sep=''))

########################################################################
