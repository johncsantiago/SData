########################
# Load Libraries
library(tidyverse)
library(edgeR)
library(pca3d) # Requires XQuartz on macOS X
library(ggfortify)
library(directlabels)
library(plot3D)
library(rgl)
library(plot3Drgl)
library(org.Hs.eg.db)
library(gplots)
library(RColorBrewer)
x1=as.list(org.Hs.egENSEMBL2EG)
x2=as.list(org.Hs.egSYMBOL)
##using Anders PCA.R script

# User Settings

# Graph All/Lean/Fatty/Viable/Nonviable samples?
samplesToGraph = "Lean"

# Use pca2d (John's method) or ggplot2 (Anders' method)?
# Use 3D for pca3d
typeOfPCA = "ggplot2"

# IF ggplot2, what style of plot? Shape or Label points?
styleOfPCA = "Shape"

# IF ggplot2, What Principal Components? (Default is X=1, Y=2)
PCX = 1
PCY = 2



##Run the script
########################
# Automatically set working directory
currentwd <- dirname(rstudioapi::getSourceEditorContext()$path)
githubwd <- str_replace(currentwd, "R Scripts/PCA", "")
setwd("/Users/johncsantiago/Documents/GitHub/SandersLab/")

########################
# Read input files
countdata=read.csv("RNA-Seq Data Files/mgh_raw_countdata.csv",row.names=1)
groups=read.csv("RNA-Seq Data Files/mgh_metadata.csv",row.names=1)

########################
# Normalize and remove low-count reads

#Saves normalized CPM file and a count table both with low count reads removed
countdata <- countdata[,row.names(groups)]
x <- DGEList(countdata)
#class(x)
geneid <- rownames(x)
x$samples <- merge(x$samples,groups,by="row.names")
row.names(x$samples) <- x$samples[,1]
x$samples <- x$samples[,-1]

##at least is the cpm of a gene with 10 reads
##10/smallest library size=x/1 000 000
##10 000 000/smallest library size = x = smallest acceptable cpm for a gene to be kept
atleast <- (10000000/min(x$samples$lib.size))
keep <- rowSums(cpm(x)>=atleast)>=3
y <- x[keep, ,keep.lib.sizes=FALSE]

#Normalize via TMM method
z <- calcNormFactors(y, method = "TMM") 

##using z gives the normalized cpm, using y gives the not normalized cpm
cpmsc <- cpm(z)
pcadata <- cpmsc

########################
# Set up individual groups

# Relabel for graph display
groups$fat <- factor(groups$fat, labels=c("Fatty", "Lean"))
groups$viability <- factor(groups$viability, levels=c("V","N"), labels=c("Viable","Nonviable"))
groups$timepoint <- factor(groups$timepoint, levels=c("pre","three","six"), labels=c("0h","3h","6h"))

#Set up blank plot variable
plot <- NA

# For all samples
if (samplesToGraph=="All"){
  # Nothing needs doing
}
# For Lean or Fatty
if (samplesToGraph=="Lean"||samplesToGraph=="Fatty"){
  groups <- subset(groups, fat == samplesToGraph)
  pcadata <- pcadata[,row.names(groups)[groups[,1]==samplesToGraph]]
  
  #Need to exclude 1 gene that has 0 expression in Fatty
  pcadata=pcadata[row.names(pcadata)!="ENSG00000126890",]
  pcadata=pcadata[row.names(pcadata)!="CTAG2_ENSG00000126890",]
}
# For Viable or Nonviable
if (samplesToGraph=="Viable"||samplesToGraph=="Nonviable"){
  groups <- subset(groups, viability == samplesToGraph)
  pcadata <- pcadata[,row.names(groups)[groups[,2]==samplesToGraph]]
}

# Run the actual PCA calculation
pca <- prcomp(t(pcadata), scale.=TRUE) 

########################
# ALL
if (samplesToGraph=="All"){
  if (typeOfPCA=="pca2d"){
    gr <- factor(groups[,7])
    pcacols=c(1,1,1,2,2,2,3,3,3,4,4,4)
    par(xpd=T,mfrow=c(1,1))
    #PC1 and PC2
    pca2d(pca, group=gr, show.group.labels = T, palette=pcacols,bg=NA,components = c(1,2),show.ellipses = F,title = metatype[i],radius = .1,legend=NULL)
    #PC1 and PC3
    pca2d(pca, group=gr, show.group.labels = T, palette=pcacols,bg=NA,components = c(1,3),show.ellipses = F,title = metatype[i],radius = .1,legend=NULL)
    #Bar plot of PC distribution
    eigs <- pca$sdev^2
    ve=(eigs / sum(eigs))*100
    names(ve)=c("PC1","PC2","PC3")
    barplot(ve,ylim=c(0,100))
  }
  if (typeOfPCA=="3D") {
    gr <- factor(groups[,7])
    pcacols=c(1,1,1,2,2,2,3,3,3,4,4,4)
    par(xpd=T,mfrow=c(1,1))
    pca3d(pca, group=gr,show.group.labels=T)
  }
  if (typeOfPCA=="ggplot2"){
    if (styleOfPCA=="Shape"){
      #Create plot with this criteria:
      #Colored based on replicate group (i.e. LV1)
      #Group labeled at center point of the 3 replicates
      #Shaped based on fatty (triagle is lean, round is fatty)
      plot <- autoplot(pca, x=PCX, y=PCY, data=groups, colour = 'Group', shape = 'fat') + 
        labs(shape="Liver Type") +
        ggtitle(paste("PCA - ", samplesToGraph, " Samples", sep=""))
      # Use Directlabels to put the group names on the plot  
      plot <- direct.label.ggplot(plot, method=smart.grid)
    }
    # Label each point for individual sample tracking
    if (styleOfPCA=="Label"){
      plot <- autoplot(pca, x=PCX, y=PCY, data = groups, colour = 'Group', shape = FALSE, label.size = 7)+
        labs(color = "Group")+
        ggtitle(paste("PCA - ", samplesToGraph, " Samples (Labeled)", sep=""))
    }
  }
}

########################
# LEAN OR FATTY PCAs
if (samplesToGraph=="Lean"||samplesToGraph=="Fatty"){
  if (typeOfPCA=="pca2d"){
    gr <- factor(groups[groups[,1]==samplesToGraph,7])
    pcacols=c(1,1,1,2,2,2)
    par(xpd=T,mfrow=c(1,1))
    pca2d(pca, group=gr, show.group.labels = T, palette=pcacols,bg=NA,components = c(1,2),show.ellipses = F,title = metatype[i],radius = .1,legend=NULL)
    pca2d(pca, group=gr, show.group.labels = T, palette=pcacols,bg=NA,components = c(1,3),show.ellipses = F,title = metatype[i],radius = .1,legend=NULL)
    eigs <- pca$sdev^2
    ve=(eigs / sum(eigs))*100
    names(ve)=c("PC1","PC2","PC3")
    barplot(ve,ylim=c(0,100))
  }
  if (typeOfPCA=="3D") {
    gr <- factor(groups[groups[,1]==samplesToGraph,7])
    pcacols=c(1,1,1,2,2,2)
    par(xpd=T,mfrow=c(1,1))
    pca3d(pca, group=gr,show.group.labels=T)
  }
  if (typeOfPCA=="ggplot2"){
    #Shapes (circle, triangle, square) for each point
    if (styleOfPCA=="Shape"){
      plot <- autoplot(pca, x=PCX, y=PCY, data=groups)+ # What data to use for the chart
        geom_point(aes(colour=groups$viability, shape=groups$timepoint, size=1))+ # What style of graph to make
        geom_dl(aes(label=NA,colour=groups$viability),method="smart.grid")+ # How to color the chart
        labs(color = "Viability")+ # Label for viability legend
        labs(shape = "Timepoint")+ # Label for timepoint legend
        ggtitle(paste("PCA - All ", samplesToGraph, " Samples", sep=""))+ # Title of chart
        scale_color_manual(values=c("green4","red"))+ # Manually set the colors of the nodes, sequentially
        theme_bw()+ # Change background theme
        scale_size(guide = 'none') # Hide the size from the legend
    }
    # Label each point for individual sample tracking
    if (styleOfPCA=="Label"){
      plot <- autoplot(pca, x=PCX, y=PCY, data = groups, colour = 'viability', shape = FALSE, label.size = 7)+
        labs(color = "Viability")+
        ggtitle(paste("PCA - All ", samplesToGraph, " Samples (Labeled)", sep=""))
    }
  }
}

########################
# VIABLE OR NONVIABLE PCAs
if (samplesToGraph=="Viable"||samplesToGraph=="Nonviable"){
  if (typeOfPCA=="pca2d"){
    gr <- factor(groups[groups[,2]==samplesToGraph,7])
    pcacols=c(1,1,1,2,2,2)
    par(xpd=T,mfrow=c(1,1))
    pca2d(pca, group=gr, show.group.labels = T, palette=pcacols,bg=NA,components = c(1,2),show.ellipses = F,title = metatype[i],radius = .1,legend=NULL)
    pca2d(pca, group=gr, show.group.labels = T, palette=pcacols,bg=NA,components = c(1,3),show.ellipses = F,title = metatype[i],radius = .1,legend=NULL)
    eigs <- pca$sdev^2
    ve=(eigs / sum(eigs))*100
    names(ve)=c("PC1","PC2","PC3")
    barplot(ve,ylim=c(0,100))
  }
  if (typeOfPCA=="3D") {
    gr <- factor(groups[groups[,2]==samplesToGraph,7])
    pcacols=c(1,1,1,2,2,2)
    par(xpd=T,mfrow=c(1,1))
    pca3d(pca, group=gr,show.group.labels=T)
  }
  if (typeOfPCA=="ggplot2"){
    if (styleOfPCA=="Shape"){
      plot <- autoplot(pca, x=PCX, y=PCY, data=groups)+
        geom_point(aes(colour=groups$Group, shape=groups$fat, size=1))+
        geom_dl(aes(label=groups$Group,colour=groups$Group),method="smart.grid")+
        labs(color = "Group")+
        labs(shape = "Fattiness")+
        ggtitle(paste("PCA - All ", samplesToGraph, " Samples", sep=""))+
        theme_bw()+
        scale_size(guide = 'none')
    }
    # Label each point for individual sample tracking
    if (styleOfPCA=="Label"){
      plot <- autoplot(pca, x=PCX, y=PCY, data = groups, colour = 'Group', shape = FALSE, label.size = 7)+
        labs(color = "Group")+
        ggtitle(paste("PCA - All ", samplesToGraph, " Samples (Labeled)", sep=""))
    }
  }
}

#Print out the selected ggplot2
plot


eigs <- pca$sdev^2
ve=(eigs / sum(eigs))*100
names(ve)=c("PC1","PC2","PC3")
barplot(ve,ylim=c(0,30),cex.names =.8)
lines(x=c(0,1.5),y=c((ve[1]/10),(ve[1]/10)),col="blue",lty=2,lwd=2)
lines(x=c(0,1.5),y=c((ve[1]/10),(ve[1]/10)),col="blue",lty=2,lwd=2)

##gene contribution to each PC

contribution=pca$rotation
PC1contribute=contribution[,1]
PC1contribute=PC1contribute[order(PC1contribute,decreasing=T)]
PC2contribute=contribution[,2]
PC3contribute=contribution[,3]

##Because the sum of the squares of all loadings for an individual principal component must sum to one, we can calculate what the loadings would be if all variables contributed equally to that principal component.

##Any variable that has a larger loading than this value (positive or negative) contributes more than one variable’s worth of information and would be regarded as an important contributor to that principal component.



equalloading.cutoff=sqrt(1/nrow(contribution))
loading=abs(PC1contribute)
hist(loading, breaks = 100)
abline(v=equalloading.cutoff,col="red",lty=2)
abline(v=(2*equalloading.cutoff),col="red",lty=2)

hist(PC1contribute, breaks = 100)
abline(v=equalloading.cutoff,col="red",lty=2)
abline(v=-equalloading.cutoff,col="red",lty=2)

hist(PC3contribute, breaks = 100)
abline(v=equalloading.cutoff,col="red",lty=2)
abline(v=-equalloading.cutoff,col="red",lty=2)

PC=1

PCcontribute=contribution[,PC]
loading=abs(PCcontribute)
eigs <- pca$sdev^2
ve=(eigs / sum(eigs))*100
percent.variation=signif(ve,4)

percents=(PCcontribute*PCcontribute)*percent.variation[PC]

loading=loading[order(loading,decreasing=T)]
length(loading[loading>=equalloading.cutoff])
percents=percents[order(percents,decreasing=T)]
percent2=percents
temppercent=0
i=1
while(i<=length(percents)){
  percents[i]=percents[i]+temppercent
  temppercent=percents[i]
  i=i+1
}
plot(x=c(1:length(percents)),y=percents,main=paste("PC",PC,sep=""),ylab="Percent Variation",xlab="Cumulative Gene Contribution  to Variance")
lines(x=c(length(loading[loading>=equalloading.cutoff]),length(loading[loading>=equalloading.cutoff])),y=c(0,percents[length(loading[loading>=equalloading.cutoff])]),col="red",lty=2,lwd=3)
lines(x=c(0,length(loading[loading>=equalloading.cutoff])),y=c(percents[length(loading[loading>=equalloading.cutoff])],percents[length(loading[loading>=equalloading.cutoff])]),col="red",lty=2,lwd=3)

plot(x=c(1:length(percents)),y=percents,main=paste("PC",PC,sep=""),ylab="Percent Variation",xlab="Cumulative Gene Contribution  to Variance",cex=.1)

lines(x=c(length(loading[percents<=(percent.variation[PC]/10)]),length(loading[percents<=(percent.variation[PC]/10)])),y=c(0,(percent.variation[PC]/10)),col="red",lty=2,lwd=1.5)

lines(x=c(0,length(loading[percents<=(percent.variation[PC]/10)])),y=c((percent.variation[PC]/10),(percent.variation[PC]/10)),col="red",lty=2,lwd=1.5)
##points(x=c(1:length(percents)),y=percents,main=paste("PC",PC,sep=""),cex=.1)
abline(h=0)
abline(v=0)


percent2=as.matrix(percent2)
percent2=percent2[,c(1,1,1)]
percent2[,1]=(as.character(x2[as.character(x1[row.names(percent2)])]))
percent2[,3]=PCcontribute[row.names(percent2)]
colnames(percent2)=c("symbol",paste("PC",PC," Percent Contribution",sep=""),"Loading")

most.var=percent2[1:length(loading[loading>=equalloading.cutoff]),]

top25=most.var[1:length(percents[percents<(percent.variation[PC]/4)]),]
top10=most.var[1:length(percents[percents<(percent.variation[PC]/10)]),]

write.csv(percent2,"/Users/johncsantiago/Documents/R outputs/AllPC1contributions.csv")
write.csv(top25,"/Users/johncsantiago/Documents/R outputs/Top25Percent_PC1contributions.csv")
write.csv(top10,"/Users/johncsantiago/Documents/R outputs/Top10Percent_PC1contributions.csv")

hm10up=cpmsc[row.names(top10[as.numeric(top10[,3])>0,]),18:35]
row.names(hm10up)=top10[as.numeric(top10[,3])>0,1]
hm10down=cpmsc[row.names(top10[as.numeric(top10[,3])<0,]),18:35]
row.names(hm10down)=top10[as.numeric(top10[,3])<0,1]
hmcolors=colorpanel(100,"royalblue4","white","firebrick4")
heatmap.2(hm10up,scale="row",trace="none",Colv = F,col=hmcolors,cexRow = .1,main="Genes in top 10% contribution to PC1 \nwith positive correlation")
heatmap.2(hm10down,scale="row",trace="none",Colv = F,col=hmcolors,main="Genes in top 10% contribution to PC1 \nwith negative correlation",cexRow=.4)

          

PC=3

PCcontribute=contribution[,PC]
loading=abs(PCcontribute)
eigs <- pca$sdev^2
ve=(eigs / sum(eigs))*100
percent.variation=signif(ve,4)

percents=(PCcontribute*PCcontribute)*percent.variation[PC]

loading=loading[order(loading,decreasing=T)]
length(loading[loading>=equalloading.cutoff])
percents=percents[order(percents,decreasing=T)]
percent2=percents
temppercent=0
i=1
while(i<=length(percents)){
  percents[i]=percents[i]+temppercent
  temppercent=percents[i]
  i=i+1
}
plot(x=c(1:length(percents)),y=percents,main=paste("PC",PC,sep=""),ylab="Percent Variation",xlab="Cumulative Gene Contribution  to Variance")
lines(x=c(length(loading[loading>=equalloading.cutoff]),length(loading[loading>=equalloading.cutoff])),y=c(0,percents[length(loading[loading>=equalloading.cutoff])]),col="red",lty=2,lwd=3)
lines(x=c(0,length(loading[loading>=equalloading.cutoff])),y=c(percents[length(loading[loading>=equalloading.cutoff])],percents[length(loading[loading>=equalloading.cutoff])]),col="red",lty=2,lwd=3)


percent2=as.matrix(percent2)
percent2=percent2[,c(1,1,1)]
percent2[,1]=(as.character(x2[as.character(x1[row.names(percent2)])]))
percent2[,3]=PCcontribute[row.names(percent2)]
colnames(percent2)=c("symbol",paste("PC",PC," Percent Contribution",sep=""),"Loading")

most.var=percent2[1:length(loading[loading>=equalloading.cutoff]),]

top25=most.var[1:length(percents[percents<(percent.variation[PC]/4)]),]
top10=most.var[1:length(percents[percents<(percent.variation[PC]/10)]),]

write.csv(percent2,"/Users/johncsantiago/Documents/R outputs/AllPC3contributions.csv")
write.csv(top25,"/Users/johncsantiago/Documents/R outputs/Top25Percent_PC3contributions.csv")
write.csv(top10,"/Users/johncsantiago/Documents/R outputs/Top10Percent_PC3contributions.csv")

hm10up=cpmsc[row.names(top10[as.numeric(top10[,3])>0,]),18:35]
row.names(hm10up)=top10[as.numeric(top10[,3])>0,1]
hm10down=cpmsc[row.names(top10[as.numeric(top10[,3])<0,]),18:35]
row.names(hm10down)=top10[as.numeric(top10[,3])<0,1]
hmcolors=colorpanel(100,"royalblue4","white","firebrick4")
heatmap.2(hm10up,scale="row",trace="none",Colv = F,col=hmcolors,cexRow = .1,main="Genes in top 10% contribution to PC3 \nwith positive correlation")
heatmap.2(hm10down,scale="row",trace="none",Colv = F,col=hmcolors,main="Genes in top 10% contribution to PC3 \nwith negative correlation",cexRow=.4)


PC=3

PCcontribute=contribution[,PC]
loading=abs(PCcontribute)
eigs <- pca$sdev^2
ve=(eigs / sum(eigs))*100
percent.variation=signif(ve,4)

percents=(PCcontribute*PCcontribute)*percent.variation[PC]

loading=loading[order(loading,decreasing=T)]
length(loading[loading>=equalloading.cutoff])
percents=percents[order(percents,decreasing=T)]
percent2=percents
temppercent=0
i=1
while(i<=length(percents)){
  percents[i]=percents[i]+temppercent
  temppercent=percents[i]
  i=i+1
}
plot(x=c(1:length(percents)),y=percents,main=paste("PC",PC,sep=""),ylab="Percent Variation",xlab="Cumulative Gene Contribution  to Variance")
lines(x=c(length(loading[loading>=equalloading.cutoff]),length(loading[loading>=equalloading.cutoff])),y=c(0,percents[length(loading[loading>=equalloading.cutoff])]),col="red",lty=2,lwd=3)
lines(x=c(0,length(loading[loading>=equalloading.cutoff])),y=c(percents[length(loading[loading>=equalloading.cutoff])],percents[length(loading[loading>=equalloading.cutoff])]),col="red",lty=2,lwd=3)


percent2=as.matrix(percent2)
percent2=percent2[,c(1,1,1)]
percent2[,1]=(as.character(x2[as.character(x1[row.names(percent2)])]))
percent2[,3]=PCcontribute[row.names(percent2)]
colnames(percent2)=c("symbol",paste("PC",PC," Percent Contribution",sep=""),"Loading")

most.var=percent2[1:length(loading[loading>=equalloading.cutoff]),]

top25=most.var[1:length(percents[percents<(percent.variation[PC]/4)]),]
top10=most.var[1:length(percents[percents<(percent.variation[PC]/10)]),]

##write.csv(percent2,"/Users/johncsantiago/Documents/R outputs/AllPC3contributions.csv")
##write.csv(top25,"/Users/johncsantiago/Documents/R outputs/Top25Percent_PC3contributions.csv")
##write.csv(top10,"/Users/johncsantiago/Documents/R outputs/Top10Percent_PC3contributions.csv")

hm10up=cpmsc[row.names(top10[as.numeric(top10[,3])>0,]),18:35]
row.names(hm10up)=top10[as.numeric(top10[,3])>0,1]
hm10down=cpmsc[row.names(top10[as.numeric(top10[,3])<0,]),18:35]
row.names(hm10down)=top10[as.numeric(top10[,3])<0,1]
hmcolors=colorpanel(100,"royalblue4","white","firebrick4")
heatmap.2(hm10up,scale="row",trace="none",Colv = F,col=hmcolors,cexRow = .1,main="Genes in top 10% contribution to PC3 \nwith positive correlation")
heatmap.2(hm10down,scale="row",trace="none",Colv = F,col=hmcolors,main="Genes in top 10% contribution to PC3 \nwith negative correlation",cexRow=.4)

```{r}

polygon3D(y=xcoords[1:3],x=ycoords[1:3],z=zcoords[1:3],colkey=F,col="royalblue3",xlim=c(min(ycoords),max(ycoords)),zlim=c(min(zcoords),max(zcoords)),ylim=c(min(xcoords),max(xcoords)),zlab="PC3 (%)",xlab="PC2 (%)",ylab="PC1 (%)",axes=T,ticktype="detailed",box=T,bty="b2", theta=50, phi=30,r=4)

polygon3D(y=xcoords[4:6],x=ycoords[4:6],z=zcoords[4:6],colkey=F,col="royalblue3",add=T,border="black")
polygon3D(y=xcoords[7:9],x=ycoords[7:9],z=zcoords[7:9],colkey=F,col="royalblue3",add=T,border="black")

polygon3D(y=xcoords[c(10:12)],x=ycoords[c(10:12)],z=zcoords[c(10:12)],colkey=F,col="gold3",add=T,border="black")
polygon3D(y=xcoords[c(13:15)],x=ycoords[c(13:15)],z=zcoords[c(13:15)],colkey=F,col="gold3",add=T,border="black")
lines3D(y=xcoords[c(16:17)],x=ycoords[c(16:17)],z=zcoords[c(16:17)],colkey=F,col="gold3",add=T,border="black")

polygon3D(y=xcoords[c(18:20)],x=ycoords[c(18:20)],z=zcoords[c(18:20)],colkey=F,col="darkgreen",add=T,border="black")
polygon3D(y=xcoords[c(21:23)],x=ycoords[c(21:23)],z=zcoords[c(21:23)],colkey=F,col="darkgreen",add=T,border="black")
polygon3D(y=xcoords[c(24:26)],x=ycoords[c(24:26)],z=zcoords[c(24:26)],colkey=F,col="darkgreen",add=T,border="black")

polygon3D(y=xcoords[c(27:29)],x=ycoords[c(27:29)],z=zcoords[c(27:29)],colkey=F,col="firebrick3",add=T,border="black")
polygon3D(y=xcoords[c(30:32)],x=ycoords[c(30:32)],z=zcoords[c(30:32)],colkey=F,col="firebrick3",add=T,border="black")
polygon3D(y=xcoords[c(33:35)],x=ycoords[c(33:35)],z=zcoords[c(33:35)],colkey=F,col="firebrick3",add=T,border="black")
##plotrgl()


legend("topleft", legend=c("Lean Viable","Fatty Viable","Lean Non-Viable","Fatty Non-Viable"),col=pcacolors[c(3,1,4,2)],lty=1,bty="n",lwd=2)

```
