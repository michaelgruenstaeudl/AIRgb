#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2019 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl")
#email = "m.gruenstaeudl@fu-berlin.de"
#version = "2019.11.20.1700"

########################################################################

library(ggplot2)
library(svglite)
library(tcltk) # For dialog boxes
library(tools) # For function 'file_path_sans_ext'
library(dplyr) # For function '%>%'
library(grid)  # For 'textGrob'
library(gridExtra) # For 'grid.arrange'

options(scipen=999) # Avoid E-notation in numbers

########################################################################

## Load Plastome Availability Table (.csv-format)
inFn1 = tk_choose.files(caption = "Select the plastome availability table (.tsv-format)")
inData1 = read.csv(inFn1, sep = "\t")

## Load Plastome Availability Table (.csv-format)
inFn2 = tk_choose.files(caption = "Select the reported IR stats table (.tsv-format)")
out_fn = paste(file_path_sans_ext(inFn2), "_", sep='')
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

## PLOTTING TOTAL NUMBER OF PLASTOMES WITH UNEQUAL IRs AS BARCHART ##

## Convert the column "CREATE_DATE" into a frequency table
# Convert to date
notEqualDatesData = as.Date(notEqualDF$CREATE_DATE, format="%Y-%m-%d")
# Tabulate
notEqualTab = table(cut(notEqualDatesData, 'month'))
# Format as data.frame
notEqualPlotData = data.frame(DATE=as.Date(names(notEqualTab)), FREQ_RECORDS=as.vector(notEqualTab))
## Add a column that displays the growing cumulative sum
notEqualPlotData = notEqualPlotData %>% mutate(CUMFREQ=cumsum(FREQ_RECORDS))

####################################

notEqual_basePlot = ggplot(data=notEqualPlotData, aes(x=DATE, y=CUMFREQ), width=1) + 
    geom_bar(stat="identity", position="identity", fill="grey50", alpha=0.5) #+ 
    #geom_line(color="grey", alpha=0.5)

notEqual_finalPlot = notEqual_basePlot + 
    xlab("\nYear") + 
    ylab("Total Number of Records\n") + 
    ggtitle("Total number of complete plastid genome sequences available on NCBI GenBank in each year,\nwhose reported annotations of the IR do not have identical lengths",
            subtitle="Note: Only plastid genomes of angiosperms are counted") + 
    scale_x_date(limits=c(as.Date("2000-01-01"), as.Date("2020-01-01")),
                 date_breaks="1 year", minor_breaks=NULL, expand=expand_scale(0),
                 date_labels="%Y") + 
    #scale_y_continuous(breaks=seq(0, 6000, 1000), minor_breaks=seq(500, 5500, 1000)) +
    #scale_colour_grey(aesthetics = "fill") + 
    #scale_fill_brewer(palette="Dark2", name="Criterion positive/negative") + 
    scale_fill_manual(values=c("grey50", "grey0"), name="Plastid genome\navailable", labels=c("Yes", "No")) +
    #theme_bw() + 
    theme_minimal() + 
    theme(plot.title = element_text(size=20),
          plot.subtitle = element_text(size=16, face="italic"),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          plot.margin=unit(c(0.5,1.0,0.1,0.1),"cm"),  # Note: margin(t=0, r=0, b=0, l=0, unit="pt")
          legend.key.width=unit(1,"cm"))

####################################

PlastomesWithUnequalReportedIRLengths = notEqual_finalPlot
save(PlastomesWithUnequalReportedIRLengths, file="VIZ_PlastomesWithUnequalReportedIRLengths.Rdata")


####################################

#svglite(file=paste(out_fn, "VIZ_PlastomesWithUnequalReportedIRLengths.svg", sep=''), width=21, height=7.425)
#notEqual_finalPlot
#dev.off()

########################################################################

## PLOTTING IR LENGTH DIFFERENCES AS HISTOGRAM ##

## Add a column that displays the number of nucleotide differences between the IR lengths
notEqualDF = notEqualDF %>% mutate(IR_LEN_DIFF=abs(notEqualDF$IRa_REPORTED_LENGTH-notEqualDF$IRb_REPORTED_LENGTH))

####################################

LenDiff_basePlot = ggplot(data=notEqualDF, aes(x=IR_LEN_DIFF)) + 
    geom_histogram(fill="grey50", alpha=0.5) #+ 
    #geom_density(aes(y=..density..*n), color="grey0")

LenDiff_finalPlot = LenDiff_basePlot + 
    xlab("\nLength Difference") + 
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

IRLengthDiffDistribution = LenDiff_finalPlot
save(IRLengthDiffDistribution, file="VIZ_IRLengthDiffDistribution.Rdata")

####################################

#svglite(file=paste(out_fn, "VIZ_IRLengthDiffDistribution.svg", sep=''), width=21, height=7.425)
#notEqual_finalPlot
#dev.off()

########################################################################

## COMBINING SUBPLOTS INTO LARGE PLOT
svglite(file="VIZ_ComparsionOfReportedIRLengths.svg", width=21, height=7.425)
grid.arrange(notEqual_finalPlot, LenDiff_finalPlot, ncol=2, widths=c(2,1),
             top=textGrob("Comparsion of reported IR lengths", gp=gpar(fontsize=16, font=2))
            )
dev.off()
