#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2019 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl")
#email = "m.gruenstaeudl@fu-berlin.de"
#version = "2019.11.19.1130"

########################################################################

library(ggplot2)
library(svglite)
library(tcltk) # For dialog boxes
library(tools) # For function 'file_path_sans_ext'
library(dplyr) # For function '%>%'

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

## Ensure that numeric values are recognized as numeric (i.e., transform from factor to numeric via character)
## See: https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-integer-numeric-without-loss-of-information

combData = transform(combData, 
IRa_REPORTED_START = as.integer(as.character(IRa_REPORTED_START)),
IRb_REPORTED_START = as.integer(as.character(IRb_REPORTED_START)),
IRa_REPORTED_END = as.integer(as.character(IRa_REPORTED_END)),
IRb_REPORTED_END = as.integer(as.character(IRa_REPORTED_END)),
IRa_REPORTED_LENGTH = as.integer(as.character(IRa_REPORTED_LENGTH)),
IRb_REPORTED_LENGTH = as.integer(as.character(IRa_REPORTED_LENGTH))
)

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
## Resort plot data
plotData = plotData[order(plotData$DATE),]
## Add a column that displays the growing cumulative sum for each criterion
plotData = plotData %>% group_by(CRITERION) %>% mutate(CUMFREQ=cumsum(FREQ_RECORDS))

########################################################################

base_plot = ggplot(data=plotData, aes(x=DATE, y=CUMFREQ, fill=CRITERION), width=1) + 
    geom_bar(stat="identity", position="stack", alpha=0.5) # stacked barcharts
    #geom_bar(stat="identity", position="dodge", alpha=0.5) # side-by-side barcharts

myPlot = base_plot + 
    xlab("\nYear") + 
    ylab("Total Number of Records\n") + 
    ggtitle("Total number of complete plastid genome sequences available on NCBI GenBank in each year, \nseparated by presence of IR annotation",
            subtitle="Note: Only plastid genomes of angiosperms are counted") + 
    scale_x_date(limits=c(as.Date("2000-01-01"), as.Date("2020-01-01")),
                 date_breaks="1 year", minor_breaks=NULL, expand=expand_scale(0),
                 date_labels="%Y") + 
    scale_y_continuous(breaks=seq(0, 6000, 1000), minor_breaks=seq(500, 5500, 1000)) +
    #scale_colour_grey(aesthetics = "fill") + 
    #scale_fill_brewer(palette="Dark2", name="Criterion positive/negative") + 
    scale_fill_manual(values=c("grey50", "grey0"), name="IR annotation\npresent", labels=c("Yes", "No")) + 
    #theme_bw() + 
    theme_minimal() + 
    theme(plot.title = element_text(size=20),
          plot.subtitle = element_text(size=16, face="italic"),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          plot.margin=unit(c(0.5,1.0,0.1,0.1),"cm"),  # Note: margin(t=0, r=0, b=0, l=0, unit="pt")
          legend.key.width=unit(1,"cm"))

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

PlastomesWithWithoutIRAnnos = myPlot
save(PlastomesWithWithoutIRAnnos, file="VIZ_PlastomesWithWithoutIRAnnos.Rdata")

########################################################################

svglite(file="VIZ_PlastomesWithWithoutIRAnnos.svg", width=21, height=7.425)
myPlot
dev.off()

########################################################################

#    outFn = paste("Figure", taxdiv, sep="_")
#    svglite(file=paste(outFn, ".svg", sep=""), width=20, height=10)
#    grid.arrange(myPlot_regular, myPlot_transformed, ncol=2,
#                 top=textGrob(paste("Only submissions to the", taxdiv, "database have been counted.\n"), gp=gpar(fontsize=16,font=2))
#                 )
#    dev.off()
