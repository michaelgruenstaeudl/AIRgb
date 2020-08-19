library(ggplot2)
library(tcltk) # For dialog boxes
library(sqldf)
library(dplyr)

########################################################################

## Load Plastome Availability Table (.csv-format)
inFn1 = tk_choose.files(caption = "Select the plastome availability table (.tsv-format)")
plastData = read.csv(inFn1, sep = "\t")
## Load Reported IR Information Table (.csv-format)
inFn2 = tk_choose.files(caption = "Select the reported IR information table (.tsv-format)")
reportedData = read.csv(inFn2, sep = "\t")

# Only account for accessions that appear in the reported data
plastData <- sqldf("SELECT plastData.* FROM plastData, reportedData WHERE plastData.ACCESSION == reportedData.ACCESSION")

## Replace "n.a." Strings in table with proper <NA> values and cast the columns to numeric
reportedData <- mutate_all(reportedData, funs(replace(.,. == "n.a.", NA)))
reportedData$IRa_REPORTED_START <- as.integer(as.character(reportedData$IRa_REPORTED_START))
reportedData$IRb_REPORTED_START <- as.integer(as.character(reportedData$IRb_REPORTED_START))
reportedData$IRa_REPORTED_END <- as.integer(as.character(reportedData$IRa_REPORTED_END))
reportedData$IRb_REPORTED_END <- as.integer(as.character(reportedData$IRb_REPORTED_END))
reportedData$IRa_REPORTED_LENGTH <- as.integer(as.character(reportedData$IRa_REPORTED_LENGTH))
reportedData$IRb_REPORTED_LENGTH <- as.integer(as.character(reportedData$IRb_REPORTED_LENGTH))

qualityData <- data.frame(reportedData$ACCESSION)
qualityData$is_published <- !(plastData$REFERENCE == "Unpublished")
qualityData$is_cg <- reportedData$IRa_REPORTED == "yes" & reportedData$IRb_REPORTED == "yes" & (reportedData$IRa_REPORTED_LENGTH == reportedData$IRb_REPORTED_LENGTH)
qualityData$count <- 1
plotData <- aggregate(qualityData$count, by=list(qualityData$is_published, qualityData$is_cg), FUN=sum)
colnames(plotData) <- c("IS_PUBLISHED", "IS_CG", "TOTAL")

ggplot(data=plotData, mapping = aes(x=IS_PUBLISHED, y=TOTAL)) + geom_col(aes(fill=IS_CG))