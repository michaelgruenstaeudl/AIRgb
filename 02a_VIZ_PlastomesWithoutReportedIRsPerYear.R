#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2019 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl")
#email = "m.gruenstaeudl@fu-berlin.de"
#version = "2019.11.15.2000"

########################################################################

library(ggplot2)
library(svglite)
library(tcltk) # For dialog boxes
library(tools) # For function 'file_path_sans_ext'

########################################################################

## Load Plastome Availability Table (.csv-format)
inFn1 = tk_choose.files(caption = "Select the plastome availability table (.tsv-format)")
inData1 = read.csv(inFn1, sep = "\t")

## Load Plastome Availability Table (.csv-format)
inFn2 = tk_choose.files(caption = "Select the reported IR stats table (.tsv-format)")
out_fn = paste(file_path_sans_ext(inFn2), "_", sep='')
inData2 = read.csv(inFn2, sep = "\t")

## Combine the input files
combData = merge(inData1, inData2, by="ACCESSION")

########################################################################

## Seperate rows by criterion
posMatch = combData[which(combData[,"IRa_REPORTED"]=="yes"),]
negMatch = combData[which(combData[,"IRa_REPORTED"]=="no"),]

## Convert the column "CREATE_DATE" into a frequency table
# Convert to date
posDatesData = as.Date(posMatch$CREATE_DATE, format="%Y-%m-%d")
negDatesData = as.Date(negMatch$CREATE_DATE, format="%Y-%m-%d")
# Tabulate
posTab = table(cut(posDatesData, 'month'))
negTab = table(cut(negDatesData, 'month'))
# Format as data.frame
posPlotData = data.frame(DATE=as.Date(names(posTab)), FREQ_RECORDS=as.vector(posTab), CRITERION="positive")
negPlotData = data.frame(DATE=as.Date(names(negTab)), FREQ_RECORDS=as.vector(negTab), CRITERION="negative")
## Concatenate plot data parts
plotData = rbind(posPlotData, negPlotData)

########################################################################

base_plot = ggplot(data=plotData, aes(x=DATE, y=cumsum(plotData[,"FREQ_RECORDS"])), width=1, fill=CRITERION) + 
    geom_bar(stat="identity", position="dodge", alpha=0.5) ##side-by-side barcharts
    #geom_bar(stat="identity", position="identity", alpha=0.25) + ##overlaid barcharts
    #geom_bar(stat="identity", alpha=0.5) + ##stacked barcharts

myPlot = base_plot + 
    xlab("\nYear") + 
    ylab("Total Number of Records\n") + 
    ggtitle("Number of complete plastid genome sequences available on NCBI GenBank per submission year, \nseparated by presence of annotation for the IR",
            subtitle="Note: Only plastid genomes of angiosperms are counted") + 
    scale_x_date(limits=c(as.Date("2000-01-01"), as.Date("2020-01-01")),
                 date_breaks="1 year", minor_breaks=NULL, expand=expand_scale(0),
                 date_labels="%Y") + 
    #scale_fill_brewer(palette="Dark2", name="Criterion positive/negative") + 
    #theme_bw() + 
    theme_minimal() + 
    theme(plot.title = element_text(size=20),
          plot.subtitle = element_text(size=16, face="italic"),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          plot.margin=unit(c(0.5,1.0,0.1,0.1),"cm"))  # Note: margin(t=0, r=0, b=0, l=0, unit="pt")

################################

#    myPlot_transformed = ggplot(data=inData, aes(x=factor(year), y=accessions, fill=voucher)) +
#        geom_bar(stat="identity", position="dodge", alpha=0.5) + ##side-by-side barcharts
#        #geom_bar(stat="identity", position="identity", alpha=0.25) + ##overlaid barcharts
#        #geom_bar(stat="identity", alpha=0.5) + ##stacked barcharts
#        theme_bw() +
#        scale_fill_brewer(palette="Dark2", name="Specimen voucher") + 
#        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                      labels = trans_format("log10", math_format(10^.x))) +
#        xlab("\nYear") + 
#        ylab("Number of accessions\n") +
#        labs(#tag="Log-transformed y-axis",
#             title="Accession numbers of new submissions of \ngenomic DNA to ENA per submission year",
#             subtitle="Log-transformed y-axis") #+
#        #ggtitle("Accession numbers of new submissions of genomic DNA to ENA per submission year") +
#        #ggsubtitle("Only submissions to the plant database have been counted.\nNote: Y-axis has been square-root transformed for better visibility")
#    myPlot_transformed = myPlot_transformed + annotation_logticks(base=10, sides="l")


########################################################################

svglite(file=paste(out_fn, "VIZ_PlastomesWithWithoutIRAnnos.svg", sep=''), width=17, height=11.5)
myPlot
dev.off()

########################################################################

#    outFn = paste("Figure", taxdiv, sep="_")
#    svglite(file=paste(outFn, ".svg", sep=""), width=20, height=10)
#    grid.arrange(myPlot_regular, myPlot_transformed, ncol=2,
#                 top=textGrob(paste("Only submissions to the", taxdiv, "database have been counted.\n"), gp=gpar(fontsize=16,font=2))
#                 )
#    dev.off()
