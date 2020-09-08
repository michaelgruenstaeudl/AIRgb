#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2019-2020 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl")
#email = "m.gruenstaeudl@fu-berlin.de"
#version = "2020.09.08.1930"

########################################################################

library(svglite)
library(tcltk) # For dialog boxes
#library(grid)  # For 'textGrob'
#library(gridExtra) # For 'grid.arrange'
library(ggpubr)

options(scipen=999) # Avoid E-notation in numbers

########################################################################

## Load plots
PlastomeNumbersPerYear__angiosperms = readRDS(tk_choose.files(caption = "Select Rds file/plot '01_FIGURE_PlastomeNumbersPerYear__angiosperms'"))
PlastomeNumbersPerYear__nonangiosperms = readRDS(tk_choose.files(caption = "Select Rds file/plot '01_FIGURE_PlastomeNumbersPerYear__nonangiosperms'"))

## Arrange plots
myFigure = ggarrange(
    PlastomeNumbersPerYear__angiosperms,
    PlastomeNumbersPerYear__nonangiosperms,
    nrow=1, ncol=2, 
    labels=c("(a)","(b)"), font.label=list(size=24)
)

########################################################################

svglite(file="MANUSCRIPT_Figure1.svg", width=21, height=9.8) #height=29.7)
annotate_figure(myFigure,
    top=text_grob("Total number of complete plastid genomes on NCBI GenBank\nper year", size=20, face="bold"),
    fig.lab = "Figure 1", fig.lab.face = "bold", fig.lab.size=24
)
dev.off()

########################################################################
