#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD <mi.gruenstaeudl@gmail.com>"
#copyright = "Copyright (C) 2014-2015 Michael Gruenstaeudl"
#version = "2015.11.17.1600"

## Load libraries
library(ggplot2)

## Load files
setwd("/home/michael_science/Desktop/")
GBrecords = read.csv(file.choose())
Nsubmtools = read.csv(file.choose())

## Plot 1
p1 = ggplot(data=Nsubmtools, aes(x=year, y=total)) +
    geom_line(colour = "black") +
    #xlim=10 + 
    #scale_x_discrete(breaks=c(10000,20000), labels=c(10000,20000)) +
    scale_y_continuous(limit = c(0,10), breaks=c(1:10)) +
    theme_bw() +
    xlab("") + 
    ylab("Total Number of Tools\n")

p2 = ggplot(data=GBrecords, aes(x=year, y=count, width=1)) +
    geom_bar(stat="identity", position = "identity", colour = "black", fill = "gray") +
    #scale_x_discrete(breaks=c(10000,20000), labels=c(10000,20000)) +
    #scale_y_discrete(breaks=GBrecords[,1], labels=GBrecords[,1]) +
    theme_bw() +
    xlab("\nYears") + 
    ylab("Number of New Records\n")

#ggplot(data=woMPI, aes(x=nreps, y=time_h, group=sim, linetype=mpiBool)) +
#    geom_point() +
#    geom_line() +
#    geom_text(data=woMPI_labels, aes(label=sim), size=3, vjust=-0.5, hjust=1.0) +
#    geom_point(data=MPI) +
#    geom_line(data=MPI) +
#    geom_text(data=MPI_labels, aes(label=sim), size=3, vjust=-0.5, hjust=1.0) +
#    scale_x_discrete(breaks=c(1,5,10,20,50,100), labels=c(1,5,10,20,50,100)) +
#    scale_y_discrete(breaks=c(1, seq(0,75,5)[2:16]), labels=c(1, seq(0,75,5)[2:16])) +
#    theme_bw() +
#    xlab("\nNumber of Replicates") + 
#    ylab("Computation Time (hours)\n") +
#    ggtitle("Computation Time by Number of Replicates and\nApplication of OpenMPI\n") +
#    scale_linetype_manual(values = c("dotted", "solid"))



## Save plot
svg("/home/michael_science/Desktop/Figure1.svg", width=6, height=6)

## From: https://github.com/hadley/ggplot2/wiki/Align-two-plots-on-a-page
## Convert plots to gtable objects
library(gtable)
library(grid) # low-level grid functions are required
g1 <- ggplotGrob(p1)
#g1 <- gtable_add_cols(g1, unit(0,"mm")) # add a column for missing legend
g2 <- ggplotGrob(p2)
g <- rbind(g1, g2, size="first") # stack the two plots
g$widths <- unit.pmax(g1$widths, g2$widths) # use the largest widths
# center the legend vertically
g$layout[grepl("guide", g$layout$name),c("t","b")] <- c(1,nrow(g))
grid.newpage()
grid.draw(g)

dev.off()


## Save plot
#svg("/home/michael_science/Desktop/Figure1.svg", width=6, height=6)
#plot
#dev.off()
