#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2019 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl")
#email = "m.gruenstaeudl@fu-berlin.de"
#version = "2019.11.18.1100"

########################################################################

library(ggplot2)
library(svglite)
library(tcltk) # For dialog boxes
library(tools) # For function 'file_path_sans_ext'
library(dplyr) # For function '%>%'

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
## Add a column that displays the growing cumulative sum
plotData = plotData %>% mutate(CUMFREQ=cumsum(FREQUENCY))

########################################################################

base_plot = ggplot(data=plotData, aes(x=DATE, y=CUMFREQ), width=1) + 
    geom_bar(stat="identity", position="identity", fill="grey50", alpha=0.5) #+ 
    #geom_line(color="grey", alpha=0.5)

myPlot = base_plot + 
    xlab("\nYear") + 
    ylab("Total Number of Records\n") + 
    ggtitle("Total number of complete plastid genome sequences available on NCBI GenBank in each year",
            subtitle="Note: Only plastid genomes of angiosperms are counted") + 
    scale_x_date(limits=c(as.Date("2000-01-01"), as.Date("2020-01-01")),
                 date_breaks="1 year", minor_breaks=NULL, expand=expand_scale(0),
                 date_labels="%Y") + 
    scale_y_continuous(breaks=seq(0, 6000, 1000), minor_breaks=seq(500, 5500, 1000)) +
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

########################################################################

PlastomeNumbersPerYear = myPlot
save(PlastomeNumbersPerYear, file="VIZ_PlastomeNumbersPerYear.Rdata")

########################################################################

svglite(file="VIZ_PlastomeNumbersPerYear.svg", width=21, height=7.425)
myPlot
dev.off()

########################################################################
