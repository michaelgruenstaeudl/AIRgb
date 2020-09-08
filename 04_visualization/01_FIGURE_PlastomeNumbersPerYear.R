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
#start_year = 2000
start_year = 2010

########################################################################

## Load Plastome Availability Table (.csv-format)
AvailTableFn = tk_choose.files(caption = "Select the plastome availability table (.tsv-format)")
#out_fn = paste(file_path_sans_ext(AvailTableFn), "_", sep='')
out_fn = dirname(AvailTableFn)
AvailTableData = read.csv(AvailTableFn, sep = "\t")

########################################################################

## Convert the column "CREATE_DATE" into a frequency table

# Convert to date
datesData = as.Date(AvailTableData$CREATE_DATE, format="%Y-%m-%d")
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
    ylab("Cumulative Number of Records\n") +
    ggtitle("Cumulative number of complete plastid\ngenomes on NCBI GenBank per year",
            subtitle="Note: Only data after 2009 is displayed."
    ) +
    scale_x_date(
        limits=c(as.Date(paste(start_year, "-01-01", sep='')), as.Date("2020-01-01")),
        date_breaks="1 year",
        minor_breaks=NULL,
        expand=expansion(0),
        date_labels="%Y"
    ) +
    scale_y_continuous(breaks=seq(0, 6000, 1000), minor_breaks=seq(500, 5500, 1000)) +
    #scale_colour_grey(aesthetics = "fill") +
    #scale_fill_brewer(palette="Dark2", name="Criterion positive/negative") +
    scale_fill_manual(values=c("grey0", "grey50"), name="Plastid genome\navailable", labels=c("Yes", "No")) +
    #theme_bw() +
    theme_minimal() +
    theme(plot.title = element_text(size=20),
          plot.subtitle = element_text(size=16, face="italic"),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          plot.margin=unit(c(0.5,1.0,0.1,0.1),"cm"),  # Note: margin(t=0, r=0, b=0, l=0, unit="pt")
          legend.key.width=unit(1,"cm"))

########################################################################

assign(script_name, myPlot)
saveRDS(eval(as.name(script_name)), file=paste(out_fn, '/', script_name, ".Rds", sep=''))

########################################################################
