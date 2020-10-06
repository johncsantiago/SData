library(oligo)
library(affycoretools)
library(limma)
library(pd.ragene.1.0.st.v1)
setwd("/Users/johncsantiago/Documents/GitHub/SData/Joan_CEL_GCT_Files/ZIP050517EH11_(RaGene-1_0-st-v1)1hr/")


gfs <- read.celfiles(list.files())
ES <- rma(gfs)
lowExprCutoff=log(50,base=2)
pData(ES)$Cond <- factor(c("CtrlF", "CtrlM", "TreatF", "TreatM","CtrlF", "CtrlM", "TreatF", "TreatM","CtrlF", "CtrlM", "TreatF", "TreatM"))
ES <- annotateEset(ES, "pd.ragene.1.0.st.v1")
omit <- is.na(fData(ES)$SYMBOL)
ES <- ES[!omit,]
ES <- ES[apply(exprs(ES), 1, var) > 0.01, ]
maxvals <- apply(exprs(ES), 1, max)
ES <- ES[maxvals > lowExprCutoff, ]


design <- model.matrix(~ 0 + pData(ES)$Cond)
colnames(design) <- levels(pData(ES)$Cond)
fit <- lmFit(ES, design=design)

fit2 <- contrasts.fit(fit, c(0,-1,0,1))
fit2 <- eBayes(fit2)
Result <- topTable(fit2, number=nrow(ES), sort.by="none")
Result <- cbind(Result, exprs(ES))
names(Result) <- sub("adj.P.Val", "FDR", names(Result))
names(Result) <- sub(".CEL$", "", names(Result))
Result <- Result[order(Result$P.Value, decreasing=F), ]
Result$GENENAME[Result$FDR<.05]

MaleEffect=Result

library(pca3d)
library(plot3D)
groups <- colnames(Result)[11:22]
pca <- prcomp(t(Result[,11:22]), scale.=TRUE) 
gr=factor(groups)

eigs <- pca$sdev^2
ve=(eigs / sum(eigs))*100
names(ve)=c("PC1","PC2","PC3")
barplot(ve,ylim=c(0,20),cex.names =.8)
dev.off()
xcoords=pca$x[,1]
ycoords=pca$x[,2]
zcoords=pca$x[,3]

text2D(x=xcoords,y=ycoords, group=gr, labels = gr, col=c("green3","red3","blue3","gold","green3","red3","blue3","gold","green3","red3","blue3","gold"),xlab="PC1",ylab="PC2")

lines2D(x=xcoords[c(1,5,9)],y=ycoords[c(1,5,9)],col="green3",add=T,lty=1)
lines2D(x=xcoords[c(2,6,10)],y=ycoords[c(2,6,10)],col="red3",add=T,lty=1)
lines2D(x=xcoords[c(3,7,11)],y=ycoords[c(3,7,11)],col="blue3",add=T,lty=1)
lines2D(x=xcoords[c(4,8,12)],y=ycoords[c(4,8,12)],col="gold",add=T,lty=1)

abline(h=0,lty=2)
abline(v=0,lty=2)

