#!/usr/bin/R
#author = "Tilman Mehl, Michael Gruenstaeudl"
#copyright = "Copyright (C) 2019-2020 Tilman Mehl, Michael Gruenstaeudl"
#contributors = c("Tilman Mehl", "Michael Gruenstaeudl")
#email = "m.gruenstaeudl@fu-berlin.de"
#version = "2020.09.08.1930"

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

#SMALL HELPER FUNCTION
equalityCheck <- function(lenA, lenB, tolerance){
        if(missing(tolerance)){
                ifelse((lenA == lenB), TRUE, FALSE)
        }else{
                if ((tolerance <= 0) |  (tolerance > 1)){
                        stop("Argument tolerance expects value between 0 and 1")
                }else{
                        ifelse((lenA == lenB) | ((lenA < lenB) & ((lenA + lenA * tolerance) >= lenB)) | ((lenA > lenB) & ((lenB + lenB * tolerance) >= lenA)), TRUE, FALSE)
                }
        }
}

########################################################################

## Load Plastome Availability Table (.csv-format)
AvailTableFn = tk_choose.files(caption = "Select the plastome availability table (.tsv-format)")
AvailTableData = read.csv(AvailTableFn, sep = "\t")

## Load Plastome Availability Table (.csv-format)
IRTableFn = tk_choose.files(caption = "Select the reported IR stats table (.tsv-format)")
#out_fn = paste(file_path_sans_ext(IRTableFn), "_", sep='')
out_fn = dirname(IRTableFn)
IRTableData = read.csv(IRTableFn, sep = "\t")

# Note: Does "merge" delete rows in which ACCESSION is only in one of the two infiles? I believe yes.
combinedDF = merge(AvailTableData, IRTableData, by="ACCESSION")
# Alternative code via sqldf:
#library(sqldf)
#combinedDF = sqldf("SELECT AvailTableData.* FROM AvailTableData, IRTableData WHERE AvailTableData.ACCESSION == IRTableData.ACCESSION")

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

## Extracting information on PUBLICATION STATUS
DataOnPublStatus <- data.frame(combinedDF$ACCESSION)
DataOnPublStatus$IS_PUBLISHED <- !(combinedDF$REFERENCE == "Unpublished")
DataOnPublStatus$IS_CONGRUENT <- combinedDF$IRa_REPORTED == "yes" &
                                 combinedDF$IRb_REPORTED == "yes" &
                                equalityCheck(combinedDF$IRa_REPORTED_LENGTH, combinedDF$IRb_REPORTED_LENGTH)
DataOnPublStatus$count <- 1

# Obtain total numbers
plotData <- aggregate(DataOnPublStatus$count,
                      by=list(DataOnPublStatus$IS_PUBLISHED,
                              DataOnPublStatus$IS_CONGRUENT),
                      FUN=sum)
colnames(plotData) <- c("IS_PUBLISHED", "IS_CONGRUENT", "TOTAL")

# Obtain percentages
plotData <- plotData[order(plotData$IS_PUBLISHED),]
plotData$PERCENTAGE <- c(filter(plotData, IS_PUBLISHED == FALSE)$TOTAL /
                          sum(filter(plotData, IS_PUBLISHED == FALSE)$TOTAL),
                         filter(plotData, IS_PUBLISHED == TRUE)$TOTAL /
                          sum(filter(plotData, IS_PUBLISHED == TRUE)$TOTAL))
# Round the precentages to three comma positions
plotData = plotData %>% mutate_at(vars(PERCENTAGE), list(~ round(., 3)))

########################################################################

## PLOTTING OPERATIONS ##

base_plot = ggplot(data=plotData, aes(x=IS_PUBLISHED, y=PERCENTAGE), width=1) +
            geom_col(aes(fill=IS_CONGRUENT), alpha=0.5, width=0.5) +
            geom_text(data=plotData[which(plotData$IS_CONGRUENT==FALSE),], aes(label=paste("n=", TOTAL, sep="")), y=Inf, size=4.5, vjust=4) + 
            geom_text(data=plotData[which(plotData$IS_CONGRUENT==TRUE),], aes(label=paste("n=", TOTAL, sep="")), y=-Inf, size=4.5, vjust=-3)

myPlot =  base_plot +
        xlab("\nPublication Status") +
        ylab("Percentage of Records\n") +
        #ggtitle("Percentage of complete plastid\ngenomes on NCBI GenBank by publication status\nthat contain annotations for, and\nhave equality in length between, both IRs",
        #    subtitle="Note: Foo bar baz"
        #) +
        scale_x_discrete(labels=c("Unpubl.", "Publ.")) +
        #scale_y_continuous(breaks=seq(0, 6000, 1000), minor_breaks=seq(500, 5500, 1000)) +
        scale_fill_manual(values=c("grey0", "grey50"),
                          name="Presence of annotations for, and equality \nin length between, both IRs",
                          labels=c("No", "Yes")) +
        theme_minimal() +
        theme(plot.title = element_text(size=20),
              plot.subtitle = element_text(size=16, face="italic"),
              axis.text=element_text(size=14),
              axis.title=element_text(size=16, face="bold"),
              plot.margin=unit(c(0.5,1.0,0.1,0.1),"cm"),  # Note: margin(t=0, r=0, b=0, l=0, unit="pt")
              legend.key.width=unit(1,"cm"),
              legend.position = "bottom")


ggsave(file = "04a_FIGURE_IRequalityByPublStatus.svg", plot=myPlot)
ggsave(file = "04a_FIGURE_IRequalityByPublStatus.png", plot=myPlot)
########################################################################

assign(script_name, myPlot)
saveRDS(eval(as.name(script_name)), file=paste(out_fn, '/', script_name, ".Rds", sep=''))

########################################################################
