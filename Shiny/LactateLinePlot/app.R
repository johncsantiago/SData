list.of.packages <- c("shiny", "scales","gplots","colourpicker")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0){
  if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")}
  BiocManager::install("scales")
  install.packages('shiny')
  install.packages("gplots")
  install.packages("colourpicker")
  
}

library(shiny)
library(scales)
library(gplots)
library(colourpicker)
# Define UI for app that draws a histogram ----
ui <- fillPage(
  titlePanel("Customize Plot"),
  sidebarLayout(position="left",
                sidebarPanel(
                  selectInput("select", h3("Select Functional Group for Color"), 
                              choices = list("High Performance" = 1, 
                                             "Intermediate Performance" = 2,
                                             "Low Performance" = 3), 
                              selected = 1),
                  
                  conditionalPanel(
                    condition = "input.select == 1",
                    colourInput("H.line.color", "Line Color",
                                "royalblue3",showColour = 'background'),
                    ),
                  conditionalPanel(
                    condition = "input.select == 2",
                    colourInput("I.line.color", "Line Color",
                                "gold",showColour = 'background'),
                    ),
                  conditionalPanel(
                    condition = "input.select == 3",
                    colourInput("L.line.color", "Line Color",
                                "firebrick",showColour = 'background'),
                    ),
                  sliderInput("shape", h3("Shape"), min = 1, max = 5, value = 2),
                  sliderInput("size", h3("Size"), min = 0, max = 5, value = 2, step=.1),
                  sliderInput("back", h4("Background"), min = 0, max = 1, value = 1,step=.01),
                , width = 3),
                mainPanel(
                  
                  plotOutput(outputId = "distPlot", height=600, width=900)
                  
                )
  )
)


# Define server logic required to draw a histogram ----
server <- function(input, output) {

    ##Liver Trait Data 
    LTD= read.csv("https://raw.githubusercontent.com/johncsantiago/Lactate-Correlation-Manuscript/master/Data/LiverTraitData.csv")
    
    output$distPlot <- renderPlot({
      ##FunctionalGroups
      FG=c(rep(input$H.line.color,5),
           rep(input$I.line.color,5,2),
           rep(input$L.line.color,5,3))
      names(FG)=c("FV12","LV12","LV13","FV13","LV11","FN11","FN12","LN11","FN13","LN13")
      
      plot.new()
      col.pan <- colorpanel(101, "black","white")
      # Change the plot region color
      rect(par("usr")[1], par("usr")[3],
           par("usr")[2], par("usr")[4],
           col = col.pan[(input$back*100)+1]) # Color
      par(new = TRUE)
      
      x=LTD$time[LTD$libraryID==paste0(substring(names(FG)[1],1,2),substring(names(FG)[1],4,4))&LTD$time!=0&LTD$time<=6]
      y=LTD$alactate[LTD$libraryID==paste0(substring(names(FG)[1],1,2),substring(names(FG)[1],4,4))&LTD$time!=0&LTD$time<=6]
      
      plot(x,y, pch=input$shape+20, bg=FG[1],col="black",ylim=c(min(na.omit(LTD$alactate)),max(na.omit(LTD$alactate))), ylab= "Lactate (mmol/L)", xlab = "Time (Hours)", cex=input$size, cex.lab=1.5, cex.axis=1.5)
      lines(spline(x,y, n = 201,ties = "mean",method="natural"), col = FG[1], lwd=input$size)
      i=2
      while(i<=length(FG)){
        x=LTD$time[LTD$libraryID==paste0(substring(names(FG)[i],1,2),substring(names(FG)[i],4,4))&LTD$time!=0&LTD$time<=6]
        y=LTD$alactate[LTD$libraryID==paste0(substring(names(FG)[i],1,2),substring(names(FG)[i],4,4))&LTD$time!=0&LTD$time<=6]
        points(x,y, pch=input$shape+20, bg=FG[i],col="black", cex=input$size)
        lines(spline(x,y, n = 201,ties = "mean",method="natural"), col = FG[i],lwd=input$size)
        i=i+1
      }
      legend("topright",
             c("High Performance","Intermediate Performance","Low Performance"),
             fill=unique(FG),bty="n",cex=1.5, 
             text.col=col.pan[((((.6-input$back)/abs(.6-input$back))+1)*50)+1])
    
    }
    ##, height = 700, width = 900
    )
}


shinyApp(ui = ui, server = server, options)