list.of.packages <- c("shiny", "ggfortify")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0){
  install.packages('shiny')
  install.packages("ggfortify")
}

library(shiny)
library(ggfortify)
# Define UI for app that draws a histogram ----
ui <- fluidPage(
  titlePanel("Customize PCA"),
  sidebarLayout(position="left",
                sidebarPanel(
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
                                     choices = list("Head"      = "Head", 
                                                    "Thorax"    = "Thorax",
                                                    "Abdomen"   = "Abdomen"),
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
                  
                  plotOutput(outputId = "distPlot", height=600, width=600)
                  
                )
  )
)


# Define server logic required to draw a histogram ----
server <- function(input, output) {
    git.dir = "/Users/johncsantiago/Documents/GitHub/SData/FlySectionsProject/"
    cpmdata = read.csv(paste0(git.dir,"cpmdata.csv"),row.names = 1)
    groups  = read.csv(paste0(git.dir,"groups.csv"),row.names = 1)

    output$distPlot <- renderPlot({
      
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
      
      use.groups = na.omit(use.groups)
      use.data   = cpmdata[,row.names(use.groups)]
      use.data   = use.data[apply(use.data,1,sum)>0,]
    
      pca <- prcomp(t(use.data), scale.=TRUE) 
      gr <- as.data.frame(use.groups)
    
      autoplot(pca, data=gr, shape=input$shape, colour=input$color, size=input$size)
    })
}


shinyApp(ui = ui, server = server, options)