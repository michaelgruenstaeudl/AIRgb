#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2019 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl")
#email = "m.gruenstaeudl@fu-berlin.de"
#version = "2019.11.14.1300"

########################################################################

library(ggplot2)
library(svglite)
library(tcltk) # For dialog boxes
library(tools) # For function 'file_path_sans_ext'

########################################################################

## Load Plastome Availability Table (.csv-format)
inFn = tk_choose.files(caption = "Select the plastome availability table (.tsv-format)")
out_fn = paste(file_path_sans_ext(inFn), "_", sep='')
inData = read.csv(inFn, sep = "\t")

########################################################################

## Convert the column "CREATE_DATE" into a frequency table

# Convert to date
datesData = as.Date(inData$CREATE_DATE, format="%Y-%m-%d")
# Tabulate
tab = table(cut(datesData, 'month'))
# Format as data.frame
plotData = data.frame(DATE=as.Date(names(tab)), FREQUENCY=as.vector(tab))

########################################################################

base_plot = ggplot(data=plotData, aes(x=DATE, y=cumsum(plotData[,"FREQUENCY"])), width=1) + 
    geom_bar(stat="identity", position="identity", color="grey", fill="grey", alpha=0.5) + 
    #geom_line(color="grey", alpha=0.5)

myPlot = base_plot + 
    xlab("\nYear") + 
    ylab("Total Number of Records\n") + 
    ggtitle("Number of complete plastid genome sequences available on NCBI GenBank per submission year",
            subtitle="Note: Only plastid genomes of angiosperms are counted") + 
    scale_x_date(limits=c(as.Date("2000-01-01"), as.Date("2025-01-01")),
                 date_breaks="1 year", minor_breaks=NULL, expand=expand_scale(0),
                 date_labels="%Y") + 
    #theme_bw() + 
    theme_minimal() + 
    theme(plot.title = element_text(size=20),
          plot.subtitle = element_text(size=16, face="italic"),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          plot.margin=unit(c(0.5,1.0,0.1,0.1),"cm"))  # Note: margin(t=0, r=0, b=0, l=0, unit="pt")

########################################################################

svglite(file=paste(out_fn, "VIZ_PlastomeNumbersPerYear.svg", sep=''), width=17, height=11.5)
myPlot
dev.off()

########################################################################
