#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD"
#copyright = "Copyright (C) 2019-2020 Michael Gruenstaeudl"
#contributors = c("Michael Gruenstaeudl")
#email = "m.gruenstaeudl@fu-berlin.de"
#version = "2020.09.08.1930"

########################################################################

library(tcltk) # For dialog boxes
library(tools) # For function 'file_path_sans_ext'
library(dplyr) # For function '%>%'
library(plyr) # For function 'join_all'

########################################################################

# GETTING SCRIPT NAME
args = commandArgs(TRUE)
this_script = sub(".*=", "", commandArgs()[4])
script_name = file_path_sans_ext(basename(this_script))

########################################################################

#GLOBAL VARIABLES
start_year = 2000

########################################################################

## Load Plastome Availability Table (.csv-format)
AvailTableFn = tk_choose.files(caption = "Select the plastome availability table (.tsv-format)")
AvailTableData = read.csv(AvailTableFn, sep = "\t")

## Load Plastome Availability Table (.csv-format)
IRTableFn = tk_choose.files(caption = "Select the reported IR stats table (.tsv-format)")
#out_fn = paste(file_path_sans_ext(IRTableFn), "_", sep='')
out_fn = dirname(IRTableFn)
IRTableData = read.csv(IRTableFn, sep = "\t")

## Combine the data tables such that only accessions which exist in BOTH tables are maintained
# Note: Does "merge" delete rows in which ACCESSION is only in one of the two infiles? I believe yes.
combinedDF = merge(AvailTableData, IRTableData, by="ACCESSION")

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

## SUMMARIZING FREQUENCY OF ALL PLASTID GENOMES PER YEAR
## Convert to date
datesDF = as.Date(combinedDF$CREATE_DATE, format="%Y-%m-%d")
# Tabulate
allPlastomes_Tab = table(cut(datesDF, 'year'))
# Format as data.frame
allPlastomes_DF = data.frame(DATE=as.Date(names(allPlastomes_Tab)), N_NEW_RECS=as.vector(allPlastomes_Tab))
## Add a column that displays the growing cumulative sum
allPlastomes_DF = allPlastomes_DF %>% mutate(CUM_N_RECS=cumsum(N_NEW_RECS))

########################################################################

## SUMMARIZING FREQUENCY OF PLASTID GENOMES *WITH* REPORTED IR PER YEAR
## Seperate rows by criterion
plastomesWithIR = combinedDF[which(combinedDF[,"IRa_REPORTED"]=="yes"),]
## Convert to date
plastomesWithIR_Dates = as.Date(plastomesWithIR$CREATE_DATE, format="%Y-%m-%d")
# Tabulate
plastomesWithIR_Tab = table(cut(plastomesWithIR_Dates, 'year'))
# Format as data.frame
plastomesWithIR_DF = data.frame(DATE=as.Date(names(plastomesWithIR_Tab)), N_NEW_RECS_WITHIR=as.vector(plastomesWithIR_Tab))
## Add a column that displays the growing cumulative sum
plastomesWithIR_DF = plastomesWithIR_DF %>% mutate(CUM_N_RECS_WITHIR=cumsum(N_NEW_RECS_WITHIR))

########################################################################

## SUMMARIZING FREQUENCY OF PLASTID GENOMES *WITHOUT* REPORTED IR PER YEAR
## Seperate rows by criterion
plastomesWithoutIR = combinedDF[which(combinedDF[,"IRa_REPORTED"]=="no"),]
## Convert to date
plastomesWithoutIR_Dates = as.Date(plastomesWithoutIR$CREATE_DATE, format="%Y-%m-%d")
# Tabulate
plastomesWithoutIR_Tab = table(cut(plastomesWithoutIR_Dates, 'year'))
# Format as data.frame
plastomesWithoutIR_DF = data.frame(DATE=as.Date(names(plastomesWithoutIR_Tab)), N_NEW_RECS_WITHOUTIR=as.vector(plastomesWithoutIR_Tab))
## Add a column that displays the growing cumulative sum
plastomesWithoutIR_DF = plastomesWithoutIR_DF %>% mutate(CUM_N_RECS_WITHOUTIR=cumsum(N_NEW_RECS_WITHOUTIR))

########################################################################

## Joining data frames via "join_all" (as merge is not possible, because column 'DATE' is not a logical)
results_DF = join_all(list(allPlastomes_DF, plastomesWithIR_DF, plastomesWithoutIR_DF), by="DATE", type='left')

####################################

## Write results to table
write.csv(results_DF, file=paste(out_fn, '/', script_name, ".csv", sep=''))

########################################################################
