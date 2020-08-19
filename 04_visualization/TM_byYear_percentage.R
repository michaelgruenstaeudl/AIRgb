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

plastData$CREATE_DATE <- as.Date(plastData$CREATE_DATE)

qualityData <- data.frame(reportedData$ACCESSION)
qualityData$IS_CG <- reportedData$IRa_REPORTED == "yes" & reportedData$IRb_REPORTED == "yes" & (reportedData$IRa_REPORTED_LENGTH == reportedData$IRb_REPORTED_LENGTH)
qualityData$YEAR <- format(as.Date(plastData$CREATE_DATE, format="%d/%m/%Y"),"%Y")
colnames(qualityData) <- c("ACCESSION","IS_CG", "YEAR")

qualityData$count <- 1
plotData <- aggregate(qualityData$count, by=list(qualityData$YEAR, qualityData$IS_CG), FUN=sum)
colnames(plotData) <- c("YEAR", "IS_CG", "TOTAL")
acc_per_year <- aggregate(qualityData$count, by=list(as.numeric(qualityData$YEAR)), FUN=sum)
colnames(acc_per_year) <- c("YEAR", "TOTAL")

plotData$PERCENTAGE <- as.numeric(sqldf("SELECT (plotData.TOTAL/acc_per_year.TOTAL) FROM plotData, acc_per_year WHERE plotData.YEAR = acc_per_year.YEAR")[,1])
ggplot(data=plotData, mapping = aes(x=YEAR, y=PERCENTAGE)) + geom_col(aes(fill=IS_CG)) + coord_flip() # coord_flip because otherwise the year labels overlap (can be fixed with wider bars?)