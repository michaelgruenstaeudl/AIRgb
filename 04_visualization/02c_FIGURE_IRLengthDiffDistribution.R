#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2019-2020 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl")
#email = "m.gruenstaeudl@fu-berlin.de"
#version = "2020.09.03.1230"

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
#start_year = 2000
start_year = 2010

########################################################################

## Load Plastome Availability Table (.csv-format)
AvailTableFn = tk_choose.files(caption = "Select the plastome availability table (.tsv-format)")
AvailTableData = read.csv(AvailTableFn, sep = "\t")

## Load Plastome Availability Table (.csv-format)
IRTableFn = tk_choose.files(caption = "Select the reported IR stats table (.tsv-format)")
#out_fn = paste(file_path_sans_ext(IRTableFn), "_", sep='')
out_fn = dirname(IRTableFn)
IRTableData = read.csv(IRTableFn, sep = "\t")

## Combine the data tables such that only accessions which exist in BOTH tables are maintained
# Note: Does "merge" delete rows in which ACCESSION is only in one of the two infiles? I believe yes.
combinedDF = merge(AvailTableData, IRTableData, by="ACCESSION")

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
    ylab("Cumulative Number of Occurrences\n") + 
    ggtitle("Distribution of the differences in length\nbetween the IRa and the IRb",
            subtitle="Note: only data after 2010 is displayed;\nx-axis is set to logarithmic scale"
    ) + 
    scale_x_log10(
        #limits=c(1,200000),
        breaks=c(1,2,3,4,5,10,100,1000,10000,100000),
        minor_breaks=NULL,
        expand=expansion(0),
        labels=c(1,2,3,4,5,10,100,1000,10000,100000)
    ) +
    #scale_fill_manual(values=c("grey0", "grey50"), name="Plastid genome\navailable", labels=c("Yes", "No")) +
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
