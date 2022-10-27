list.of.packages <- c("shiny", "ggfortify")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0){
  install.packages('shiny')
  install.packages("ggfortify")
}

library(shiny)
library(ggfortify)

git.dir = "https://raw.githubusercontent.com/johncsantiago/SData/master/FlySectionsProject/"
cpmdata = read.csv(paste0(git.dir,"cpmdata.csv"),row.names = 1)
groups  = read.csv(paste0(git.dir,"groups.csv"),row.names = 1)
clickdata=data.frame(Sample=1,
                    PC1=1,
                    PC2=1,
                    sdev=1,
                    Species=1, 
                    Sex=1, 
                    Treatment=1, 
                    Tissue=1,
                    Percent1=1,
                    Percent2=1)
clickdata=clickdata[0,]


# Define UI for app that draws a histogram ----
ui <- fluidPage(
  titlePanel("Customize PCA"),
  sidebarLayout(position="left", sidebarPanel(
    checkboxGroupInput("species",
                       h5("Species"),
                       choices = list("OreR;OreR" = "OO",
                                      "sm21;OreR" = "SO"),
                       selected = c("OO","SO")),
    
    checkboxGroupInput("sex",
                       h5("Sex"),
                       choices = list("Male"   = "Male",
                                      "Female" = "Female"),
                       selected = c("Male","Female")),
    
    checkboxGroupInput("treatment",
                       h5("Treatment"),
                       choices = list("Control" = "Control",
                                      "Rapa"    = "Rapa"),
                       selected = c("Control","Rapa")),
    
    checkboxGroupInput("tissue",
                       h5("Tissue"),
                       choices = list("Head"    = "Head",
                                      "Thorax"  = "Thorax",
                                      "Abdomen" = "Abdomen"),
                       selected = c("Head","Thorax","Abdomen")),
    
    selectInput("shape", h5("Shape Variable"),
                choices = list("Species"   = "Species",
                               "Sex"       = "Sex",
                               "Treatment" = "Treatment",
                               "Tissue"    = "Tissue")),
    
    selectInput("color", h5("Color Variable"),
                choices = list("Species"   = "Species",
                               "Sex"       = "Sex",
                               "Treatment" = "Treatment",
                               "Tissue"    = "Tissue")),
    
    sliderInput("size", h5("Size"), min = 0, max = 5, value = 2, step=.1),
                
    width = 3),
    mainPanel(
      plotOutput(outputId = "distPlot",
                 click = "clicked",
                 height=600, width=600),
      tableOutput("data"))
  )
)


# Define server logic required to draw a histogram ----
server <- function(input, output) {



  ##pc1.click=reactiveVal(rep(0,72))
  ##pc2.click=reactiveVal(rep(0,72))
  ##observeEvent(input$plot_click,
    ##           pc1.click(rep(input$clicker[1,1],71)))  
  ##observeEvent(input$plot_click,
    ##           pc2.click(rep(input$clicker[1,2],72)))  
  ##use.groups=groups
  ##makeReactiveBinding("use.groups")
  
  newdata=reactive({
    use.groups =  groups
    
    if(length(input$species) == 1){
      use.groups = use.groups[use.groups$Species == input$species,]
    }
    
    if(length(input$sex) == 1){
      use.groups = use.groups[use.groups$Sex == input$sex,]
    }
    
    if(length(input$treatment) == 1){
      use.groups = use.groups[use.groups$Treatment == input$treatment,]
    }
    
    if(length(input$tissue) == 1){
      use.groups = use.groups[use.groups$Tissue == input$tissue,]
    }
    
    if(length(input$tissue) == 2){
      use.groups = use.groups[use.groups$Tissue == input$tissue[1]|
                                use.groups$Tissue == input$tissue[2],]
    }
    use.data   = cpmdata[,row.names(use.groups)]
    use.data   = use.data[apply(use.data,1,sum)>0,]
    
    pca <- prcomp(t(use.data), scale.=TRUE) 
    gr <- as.data.frame(use.groups)
    PC1=scale(pca$x[,1])
    PC2=scale(pca$x[,2])
    eigs <- pca$sdev^2
    ve=signif(((eigs / sum(eigs))*100)[1:2],4)
    ##distance=sqrt(((highlight()-PC1[1])^2)+((highlight()-PC2[1])^2))
    pca.data=data.frame(Sample=row.names(gr),
                        PC1=scale(pca$x[,1]),
                        PC2=scale(pca$x[,2]),
                        sdev=pca$sdev,
                        Species=gr$Species, 
                        Sex=gr$Sex, 
                        Treatment=gr$Treatment, 
                        Tissue=gr$Tissue,
                        Percent1=ve[1],
                        Percent2=ve[2])
    
    return(pca.data)
  })


  highlight=reactive({
    temp=newdata()
    clickdata<<-rbind(clickdata,nearPoints(temp,input$clicked, xvar="PC1", yvar="PC2", maxpoints=1))
    clickdata<<-clickdata[!(duplicated(clickdata) |
                            duplicated(clickdata, fromLast=TRUE)),]
    return(clickdata)
  })

  output$data <- renderTable({
    highlight()
  })

output$distPlot <- renderPlot({
  pca.data=newdata()
      ggplot(pca.data, aes(x = PC1, y = PC2)) +
        geom_point(aes_string(shape = input$shape, 
                       colour       = input$color),
                   size  = input$size,
                   alpha = 1) +
        labs(x=paste0("PC1 (",pca.data$Percent1[1],"%)"), 
             y=paste0("PC2 (",pca.data$Percent2[1],"%)")) +
        geom_point(data=highlight(), colour = "black", size = 6) +
        geom_point(data=highlight(), colour = "yellow", size = 5)
      
    }, res=96)
    
}


shinyApp(ui = ui, server = server, options)