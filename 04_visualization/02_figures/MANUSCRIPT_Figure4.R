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

## Load Plot
IRequalityByPublStatus__angiosperms = readRDS(tk_choose.files(caption = "Select file '04a_FIGURE_IRequalityByPublStatus__angiosperms.Rds'"))
IRequalityByPublStatus__nonangiosperms = readRDS(tk_choose.files(caption = "Select file '04a_FIGURE_IRequalityByPublStatus__nonangiosperms.Rds'"))

IRequalityBySeqVersion__angiosperms = readRDS(tk_choose.files(caption = "Select file '04b_FIGURE_IRequalityBySeqVersion__angiosperms.Rds'"))
IRequalityBySeqVersion__nonangiosperms = readRDS(tk_choose.files(caption = "Select file '04b_FIGURE_IRequalityBySeqVersion__nonangiosperms.Rds'"))

IRequalityByReleaseYear__angiosperms = readRDS(tk_choose.files(caption = "Select file '04c_FIGURE_IRequalityByReleaseYear__angiosperms.Rds'"))
IRequalityByReleaseYear__nonangiosperms = readRDS(tk_choose.files(caption = "Select file '04c_FIGURE_IRequalityByReleaseYear__nonangiosperms.Rds'"))


## Arrange plots
myFigure = ggarrange(
    ggarrange(
        IRequalityByPublStatus__angiosperms,
        IRequalityByPublStatus__nonangiosperms,
        IRequalityBySeqVersion__angiosperms,
        IRequalityBySeqVersion__nonangiosperms,
        nrow=1, ncol=4, 
        labels=c("(a)","(b)","(c)","(d)"), font.label=list(size=24)
    ),
    ggarrange(
        IRequalityByReleaseYear__angiosperms,
        IRequalityByReleaseYear__nonangiosperms,
        nrow=2, ncol=1, 
        labels=c("(e)","(f)"), font.label=list(size=24)
    ),
    nrow=2
)

########################################################################

# IMPORTANT: Use the "cowplot" R package for complex plot arrangements IN FUTURE #

svglite(file="MANUSCRIPT_Figure4.svg", width=21, height=29.4) #height=29.7)
annotate_figure(myFigure,
    top=text_grob("Foo bar\nbaz", size=20, face="bold"),
    fig.lab = "Figure 4", fig.lab.face = "bold", fig.lab.size=24
)
dev.off()

########################################################################
