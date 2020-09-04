library(oligo)
library(affycoretools)
library(limma)
library(pd.hugene.2.0.st)
setwd("/Users/johncsantiago/Documents/Elizabeth Files/")

gfs <- read.celfiles(
  "JT_071318_HUGENE_4A_ElDeiry_(HuGene-2_0-st).CEL", 
  "JT_071318_HUGENE_4B_ElDeiry_(HuGene-2_0-st).CEL", 
  "JT_071318_HUGENE_DMSOA_ElDeiry_(HuGene-2_0-st).CEL",
  "JT_071318_HUGENE_DMSOB_ElDeiry_(HuGene-2_0-st).CEL")

## Do RMA normalization & get ExpressionSet object:
ES <- rma(gfs)

lowExprCutoff <- quantile(exprs(ES), 0.1) ## see Below (*)

## Add sample annotations:
pData(ES)$Cond <- factor(c("Treat", "Treat", "Ctrl", "Ctrl"))

## Annotate probes with gene symbol, etc:
ES <- annotateEset(ES, "pd.hugene.2.0.st")

## Remove probes which don't map to a gene:
omit <- is.na(fData(ES)$SYMBOL)
ES <- ES[!omit,]

## Remove genes with negligible variance (0.01 is an arbitrary cutoff here,
## nothing special):
ES <- ES[apply(exprs(ES), 1, var) > 0.01, ]

## Note that not all gene symbols are unique in ES. If we only want to
## retain the probe for a gene that has the highest mean value, this is
## how to do it:
Means <- apply(exprs(ES), 1, mean)
ES <- ES[order(Means, decreasing=TRUE), ]
ES <- ES[!duplicated(fData(ES)$SYMBOL), ]
ES <- ES[order(fData(ES)$PROBEID), ]

## (*) As an example (not a recommendation), this would remove all probes
## with maximum value (across 4 samples) less than the 10th percentile of
## all expression values:
maxvals <- apply(exprs(ES), 1, max)
ES <- ES[maxvals > lowExprCutoff, ]

## These are the pieces to doo defferential expression analysis with limma:

design <- model.matrix(~ 0 + pData(ES)$Cond)
colnames(design) <- levels(pData(ES)$Cond)
fit <- lmFit(ES, design=design)
fit2 <- contrasts.fit(fit, c(-1, 1))
fit2 <- eBayes(fit2)

## Get results:
Result <- topTable(fit2, number=nrow(ES), sort.by="none")
Result <- cbind(Result, exprs(ES))
names(Result) <- sub("adj.P.Val", "FDR", names(Result))
names(Result) <- sub(".CEL$", "", names(Result))

## Sort by significance. Up- and down-regulated genes will be
## interspersed.
Result <- Result[order(Result$P.Value, decreasing=TRUE), ]

## Write full table to tab-delimited file:
write.table(Result, "limma.TreatVsCtrl.txt", sep="\t",
            row.names=FALSE)

## Note: this tab-delimited file can be opened in Excel, but it is a
## well-known problem that some gene symbols (e.g., SEP1) will get
## automatically (and irreversibly, if you save after opening) converted
## to dates, unless you are careful about how you import the file. Human
## gene symbols are currently being renamed to avoid this problem
## (too hard to fix Excel's behavior, apparently), but many annotation
## tools/files haven't caught up yet.
## See for more info:
## https://www.engadget.com/scientists-rename-genes-due-to-excel-151748790.html
## https://www.biostars.org/p/211861/