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
PlastomesWithWithoutIRAnnos__angiosperms = readRDS(tk_choose.files(caption = "Select file '02a_FIGURE_PlastomesWithWithoutIRAnnos__angiosperms.Rds'"))
PlastomesWithWithoutIRAnnos__nonangiosperms = readRDS(tk_choose.files(caption = "Select file '02a_FIGURE_PlastomesWithWithoutIRAnnos__nonangiosperms.Rds'"))

PlastomesWithUnequalReportedIRLengths__angiosperms = readRDS(tk_choose.files(caption = "Select file '02b_FIGURE_PlastomesWithUnequalReportedIRLengths__angiosperms.Rds'"))
PlastomesWithUnequalReportedIRLengths__nonangiosperms = readRDS(tk_choose.files(caption = "Select file '02b_FIGURE_PlastomesWithUnequalReportedIRLengths__nonangiosperms.Rds'"))

IRLengthDiffDistribution__angiosperms = readRDS(tk_choose.files(caption = "Select file '02c_FIGURE_IRLengthDiffDistribution__angiosperms.Rds'"))
IRLengthDiffDistribution__nonangiosperms = readRDS(tk_choose.files(caption = "Select file '02c_FIGURE_IRLengthDiffDistribution__nonangiosperms.Rds'"))


## Arrange plots
myFigure = ggarrange(
    PlastomesWithWithoutIRAnnos__angiosperms,
    PlastomesWithWithoutIRAnnos__nonangiosperms,
    PlastomesWithUnequalReportedIRLengths__angiosperms,
    PlastomesWithUnequalReportedIRLengths__nonangiosperms,
    IRLengthDiffDistribution__angiosperms,
    IRLengthDiffDistribution__nonangiosperms,
    nrow=3,  ncol=2, 
    labels=c("(a)","(b)","(c)","(d)","(e)","(f)"), font.label=list(size=24)
)

########################################################################

svglite(file="MANUSCRIPT_Figure2.svg", width=21, height=29.4) #height=29.7)
annotate_figure(myFigure,
    top=text_grob("Foo bar\nbaz", size=20, face="bold"),
    fig.lab = "Figure 2", fig.lab.face = "bold", fig.lab.size=24
)
dev.off()

########################################################################
