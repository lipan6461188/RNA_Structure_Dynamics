
library("limma")

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#           human
input_dir <- '~/lipan/DYNAMIC/icSHAPE/10-11/HEK293/R_input/'
output_dir <- '~/lipan/DYNAMIC/icSHAPE/10-11/HEK293/R_output/'
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#           mouse
input_dir <- '~/lipan/DYNAMIC/icSHAPE/10-11/mES/R_input/'
output_dir <- '~/lipan/DYNAMIC/icSHAPE/10-11/mES/R_output/'
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=

model <- function(inputFile)
{
    fileHead <- strsplit(inputFile, split='\\.')[[1]][1]
    outputFile <- paste0(fileHead, '.result')

    icTable <- read.csv(paste0(input_dir, inputFile),sep="\t",comment.char="#",head=FALSE, col.names=c('trans_id', 'site_1_based','c1f1', 'c1f2', 'c1f3', 'c1f4','c2f1', 'c2f2', 'c2f3', 'c2f4'))
    X <- matrix(c(rep(1,8), c(0,0,0,0,1,1,1,1)),ncol=2)
    icData <- icTable[,3:10]
    limmaFit <- lmFit(icData, design=X)
    limma.res <- eBayes(limmaFit)
    countSet <- topTable(limma.res, coef=2, number=9999999, sort.by='none')

    bigTable <- cbind(icTable, countSet)
    difference <- rowMeans(bigTable[,c('c1f1','c1f2','c1f3','c1f4')]) - rowMeans(bigTable[,c('c2f1','c2f2','c2f3','c2f4')])
    bigTable <- cbind(bigTable, difference)
    label <- rep('-', nrow(bigTable))
    label[bigTable[,'adj.P.Val'] < 0.05 & abs(difference) > 0.2 ] <- '*'
    bigTable <- cbind(bigTable, label)

    write.table(format(bigTable, digits=3), file=paste0(output_dir, outputFile), sep="\t", row.names=FALSE, quote=FALSE)
}

Files <- list.files(input_dir)
for(file in Files)
{
    cat('process ', file, '...\n')
    model(file)
}





