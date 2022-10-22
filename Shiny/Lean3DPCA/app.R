list.of.packages <- c("shiny", "scales","gplots","colourpicker", "edgeR","plot3D","plot3Drgl")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0){
  if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")}
  BiocManager::install("scales")
  BiocManager::install("edgeR")
  install.packages('shiny')
  install.packages("gplots")
  install.packages("colourpicker")
  install.packages('plot3D')
  install.packages('plot3Drgl')
  
}

library(shiny)
library(scales)
library(gplots)
library(colourpicker)
library(edgeR)
library(plot3D)
library(plot3Drgl)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  titlePanel("Customize Plot"),
  sidebarLayout(position="left",
                sidebarPanel(
                  selectInput("colselect", h3("Select Functional Group for Color"), 
                              choices = list("Functional" = 1, 
                                             "Nonfunctional" = 2), 
                              selected = 1),
                  
                  conditionalPanel(
                    condition = "input.colselect == 1",
                    colourInput("F.line.color", "Line Color",
                                "darkgreen",showColour = 'background'),
                    colourInput("F.fill.color", "Line Color",
                                "limegreen",showColour = 'background'),
                    ),
                  conditionalPanel(
                    condition = "input.colselect == 2",
                    colourInput("N.line.color", "Line Color",
                                "red3",showColour = 'background'),
                    colourInput("N.fill.color", "Line Color",
                                "red",showColour = 'background'),
                    ),
                  selectInput("shapeselect", h3("Select Time for Shape"), 
                              choices = list("0 Hours" = 1, 
                                             "3 Hours" = 2,
                                             "6 Hours" = 3), 
                              selected = 1),
                  conditionalPanel(
                    condition = "input.shapeselect == 1",
                    sliderInput("shape0", h3("Shape"), min = 1, max = 5, value = 1)),
                  conditionalPanel(
                    condition = "input.shapeselect == 2",
                    sliderInput("shape3", h3("Shape"), min = 1, max = 5, value = 4)),
                  conditionalPanel(
                    condition = "input.shapeselect == 3",
                    sliderInput("shape6", h3("Shape"), min = 1, max = 5, value = 2)),
                  
                  sliderInput("size", h3("Size"), min = 0, max = 5, value = 2, step=.1),
                  sliderInput("theta", "Left-Right",
                              min = -180, max = 180,
                              value = 40, step = 1,
                              animate =
                                animationOptions(interval = 100, loop = TRUE)),
                  ##sliderInput("theta", h3("Left-Right"), min = -180, max = 180, value = 40, step=1),
                  sliderInput("phi", "Up-Down",
                              min = -180, max = 180,
                              value = 20, step = 1,
                              animate =
                                animationOptions(interval = 100, loop = TRUE)),
                  ##sliderInput("phi", h3("Up-Down"), min = -180, max = 180, value = 20, step=1),

                  width = 3),
                
                mainPanel(
                  
                  plotOutput(outputId = "distPlot", height=900, width=1000)
                  
                )
  )
)


# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  countdata=read.csv("https://raw.githubusercontent.com/johncsantiago/Ohman-2021/master/Data/mgh_raw_countdata.csv",row.names=1)
  groups=read.csv("https://raw.githubusercontent.com/johncsantiago/Ohman-2021/master/Data/mgh_metadata.csv",row.names=1)
  
  ##EdgeR comparisons
  countdata=countdata[,row.names(groups)]
  countdata=countdata[,row.names(groups)]
  x <- countdata
  group <- factor(groups$Group)
  y <- DGEList(counts=x,group=group)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE] 
  z <- calcNormFactors(y, method = "TMM") 
  cpmdata=cpm(z)
  
  groups <- factor(group[18:35])
  pca <- prcomp(t(cpmdata[,18:35]), scale.=TRUE) 
  
  eigs <- pca$sdev^2
  ve=signif(((eigs / sum(eigs))*100)[1:3],4)
  names(ve)=c("PC1","PC2","PC3")
  
  zcoords=pca$x[,1]
  xcoords=pca$x[,2]
  ycoords=pca$x[,3]
    
    output$distPlot <- renderPlot({
      ##pcacolors=c("darkgreen","red3")
      lines3D(x=xcoords[c(1,4,7)],y=ycoords[c(1,4,7)],z=zcoords[c(1,4,7)],
              type="l",colkey=F,col=input$F.line.color,
              xlim=c(min(xcoords),max(xcoords)),
              zlim=c(min(zcoords),max(zcoords)),
              ylim=c(min(ycoords),max(ycoords)),
              zlab=paste0("PC1 (",ve[1],"%)"),
              xlab=paste0("PC2 (",ve[2],"%)"),
              ylab=paste0("PC3 (",ve[3],"%)"),
              lty=1,axes=T,ticktype="detailed",box=T,bty="b2", 
              theta=input$theta, phi=input$phi,r=4,lwd=input$size)
      
      lines3D(x=xcoords[c(2,5,8)],y=ycoords[c(2,5,8)],z=zcoords[c(2,5,8)],
              type="l",colkey=F,col=input$F.line.color,add=T,lty=2,lwd=input$size)
      lines3D(x=xcoords[c(3,6,9)],y=ycoords[c(3,6,9)],z=zcoords[c(3,6,9)],
              type="l",colkey=F,col=input$F.line.color,add=T,lty=3,lwd=input$size)
      
      lines3D(x=xcoords[c(10,13,16)],y=ycoords[c(10,13,16)],z=zcoords[c(10,13,16)],
              type="l",colkey=F,col=input$N.line.color,add=T,lty=1,lwd=input$size)
      lines3D(x=xcoords[c(11,14,17)],y=ycoords[c(11,14,17)],z=zcoords[c(11,14,17)],
              type="l",colkey=F,col=input$N.line.color,add=T,lty=2,lwd=input$size)
      lines3D(x=xcoords[c(12,15,18)],y=ycoords[c(12,15,18)],z=zcoords[c(12,15,18)],
              type="l",colkey=F,col=input$N.line.color,add=T,lty=3,lwd=input$size)
      
      points3D(x=xcoords[1:3],y=ycoords[1:3],z=zcoords[1:3],type="p",colkey=F,
               col=input$F.line.color,add=T,pch=input$shape0+20,
               bg=input$F.fill.color, cex=input$size)
      points3D(x=xcoords[4:6],y=ycoords[4:6],z=zcoords[4:6],type="p",colkey=F,
               col=input$F.line.color,add=T,pch=input$shape3+20,
               bg=input$F.fill.color, cex=input$size)
      points3D(x=xcoords[7:9],y=ycoords[7:9],z=zcoords[7:9],type="p",colkey=F,
               col=input$F.line.color,add=T,pch=input$shape6+20,
               bg=input$F.fill.color, cex=input$size)
      
      points3D(x=xcoords[10:12],y=ycoords[10:12],z=zcoords[10:12],type="p",
               colkey=F,col=input$N.line.color,add=T,pch=input$shape0+20,
               bg=input$N.fill.color,cex=input$size)
      points3D(x=xcoords[13:15],y=ycoords[13:15],z=zcoords[13:15],type="p",
               colkey=F,col=input$N.line.color,add=T,pch=input$shape3+20,
               bg=input$N.fill.color,cex=input$size)
      points3D(x=xcoords[16:18],y=ycoords[16:18],z=zcoords[16:18],type="p",
               colkey=F,col=input$N.line.color,add=T,pch=input$shape6+20,
               bg=input$N.fill.color,cex=input$size)

      
      legend("right", legend=c("F1","F2","F3","N1","N2","N3",NA,NA,NA,NA),col=c(input$F.line.color,input$F.line.color,input$F.line.color, input$N.line.color,input$N.line.color,input$N.line.color,NA,NA,NA,NA),lty=c(1,2,3,1,2,3,NA,NA,NA,NA),bty="n",lwd=input$size,cex=1.5)
      legend("right", legend=c(rep(NA,6),"0h","3H","6H"),
             pch=c(rep(NA,6),input$shape0+20,input$shape3+20,input$shape6+20),
             bty="n",pt.bg="black",cex=1.5)
    
    }
    ##, height = 700, width = 900
    )
}


shinyApp(ui = ui, server = server, options)