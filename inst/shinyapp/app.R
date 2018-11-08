library(shiny)
library(shinydashboard)
library(simDeNet)
library(DT)
library(GeneNet)
library(abind)

# function for generating beta distributions
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

# load in data files
data("celltype")
data("microRNAome_avg_exprs")

ui <- dashboardPage(
  dashboardHeader(title="One-Step Simulation"),
  dashboardSidebar(disable=TRUE), 
  dashboardBody(
    fluidRow(
      column(width=6, 
             box(
               title="Dataset", status="info", width = NULL,
              actionButton("generate", "Load data to your workspace")
               ),
             box(
               title="Parameters", status="info", width=NULL,
               wellPanel(style="overflow-y:scroll; max-height:800px",
                         p("One Step Simulation"),
                         p("This function simulates pure and mixed cell type expression data for two cell types, namely 1 and 2. The simulation is based on multivariate normal (MVN) distribution. Enter parameters for the variables below!"),
                         textInput(inputId="nSamp", label="Enter # of samples to simulate for each celltype: ", value="20"),
                         selectInput(inputId="muT", label="Enter cell type 1 (Target): ", c("Choose your dataset first!"), selected="Choose your dataset first!"),
                         selectInput(inputId="muN", label="Enter cell type 2 (Null): ", c("Choose your dataset first!"), selected="Choose your dataset first!"),
                         selectInput(inputId="propT", label="Choose vector of mixing proportions",
                                     c("seq", "unif", "customize", "beta distribution"), selected="seq"),
                         conditionalPanel(
                           condition="input.propT=='customize'",
                           textInput("customProp", "Enter customized prop: ")
                           ),
                         conditionalPanel(
                           condition="input.propT=='beta distribution'",
                           sliderInput("beta_mu", "Mean", min = 0, max = 1, value = 0),
                           sliderInput("beta_var", "Variance", min = 0, max = 1, value = 0)
                         ),
                         textInput(inputId="blocksizeT", label="Enter # of genes for each block (separate by comma): ", value="30,30,30"),
                         textInput(inputId="rhoT", label="Enter corr. coefficients for each block (separate by comma): ", value="0.7,0.7,0.7"),
                         selectInput(inputId="selectedGenes", label="Gene selection",
                                     c("random", 
                                       "most differentially expressed",
                                       "least differentially expressed"), selected="random"))
               )
             ), 
      column(width=6, 
             actionButton(inputId="OSSButton", label= "Run OneStepSim"),
             textInput(inputId="niter", label="Num. of Iterations: ", value="5"),
             plotOutput(outputId = "ROC.Plot", width="100%", height="400px"),
             downloadButton("dld.ROC", label="Download ROC plot")
             )

  )
  )
  )

server <- function(input, output, session) {
  
  # modal dialog asks user to choose a dataset before proceeding with application
  dataModal <- function(failed=FALSE){
    modalDialog(
      selectInput(inputId='ds', label='Select dataset to work with:', c('celltype', 'microRNAome_avg_exprs'), selected='celltype'),
      
      # if the dataset can't be found
      if(failed)
        div(tags$b("Invalid name of data object", style="color: red;")),
      
      footer = tagList(
        modalButton("Cancel"),
        actionButton("ok", "OK")
      )
    )
  }
  
  # brings up the modal for choosing datasets
  observeEvent(input$generate, {
    showModal(dataModal())
  })
  
  data <- reactive({paste(getwd(),"/data/",input$ds,'.rda', sep='')})
  vars <- reactive({load(data())}) # loads in correct data file
 
  ######################################
  ######## Choosing the dataset ########
  ######################################
  observeEvent(input$ok, {
  
    if(input$ds == "celltype") {
      ctype_choices <- as.list(ctab$Full_name) # retrieve the 'Full_Name' column in ctab table for clear names
    }
    else if(input$ds == "microRNAome_avg_exprs") {
      ctype_choices <- as.list(colnames(get(vars()))) # retrieve column names in microRNAome_avg_exprs
    }
    
    # functions updates the celltype options according to what dataset is entered
    updateSelectInput(session, "muT", label = "Enter cell type 1 (Target): ", choices = ctype_choices) 
    updateSelectInput(session, "muN", label = "Enter cell type 2 (Null): ", choices = ctype_choices)
  
    removeModal()
    }
  )

  N <-reactive(as.numeric(input$nSamp))
  num_rows <- reactive(nrow(get(vars())))
  niter <- reactive(as.numeric(input$niter))
  
  MUT <- reactive({
    if(input$ds == 'celltype'){
      get(vars())[,ctab$Fastq_file_name[which(ctab$Full_name==input$muT)]]
    }
    else if(input$ds == 'microRNAome_avg_exprs'){
      get(vars())[, input$muT]
    }
  })
  
  MUN <- reactive({
    if(input$ds == 'celltype'){
     get(vars())[,ctab$Fastq_file_name[which(ctab$Full_name==input$muN)]]
    }
    else if(input$ds == 'microRNAome_avg_exprs'){
      get(vars())[, input$muN]
    }
  })
 

  TProp <- reactive({
    if(input$propT=="seq"){
      seq(0,1, length=N())
    } 
    else if(input$propT=="unif"){
      runif(N(), 0, 1)
    } 
    else if(input$propT=="customize"){
      as.numeric(unlist(strsplit(input$customProp, ",")))
    }
    else if(input$propT == "beta distribution"){
      params <- estBetaParams(mu = input$beta_mu, var = input$beta_var)
      rbeta(N(), params$alpha, params$beta)
    }
  
  })
  
  BS <- reactive({as.numeric(unlist(strsplit(input$blocksizeT, ",")))})
  RHO <- reactive({as.numeric(unlist(strsplit(input$rhoT, ",")))})
  SG <- reactive({unlist(input$selectedGenes)})
  gene_selection_options <- c("random", "DEG", "non-DEG")
  names(gene_selection_options) <- c("random", 
                                     "most differentially expressed", 
                                     "least differentially expressed")
  SG <- reactive({gene_selection_options[unlist(input$selectedGenes)]})
  
  # our own functions:
  acor <- function(m){abs(cor(t(m)))}
  apcor <- function(x){abs(ggm.estimate.pcor(t(x), verbose=FALSE))}
  # The palette with grey:
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  observeEvent(input$OSSButton, {
    reactives <- reactiveValues(
      true.str.pooled = array(numeric(), c(num_rows(), num_rows(), 0)),
      acor.pure.pooled = array(numeric(), c(num_rows(), num_rows(), 0)),
      acor.mixed.pooled=array(numeric(), c(num_rows(), num_rows(), 0))
      )
    
    for(i in 1:niter()){
      oss <- reactive({oneStepSim(N(), MUT(), MUN(), Sigma.T=NULL, Sigma.N=NULL, prop.T=TProp(), dd=NULL, 
                                  rho=RHO(), block.size=BS(), str.type='interchangeable', select.gene=SG())})
      
      true.iter <- reactive({oss()$true.str.T})
      pure.iter <- reactive({acor(oss()$expr.pure.T)})
      mixed.iter <- reactive({acor(oss()$expr.mixed)})
      
      reactives$true.str.pooled <- abind(isolate(reactives$true.str.pooled), true.iter(), along=3)
      reactives$acor.pure.pooled <- abind(isolate(reactives$acor.pure.pooled), isolate(pure.iter()), along=3)
      reactives$acor.mixed.pooled <- abind(isolate(reactives$acor.mixed.pooled), isolate(mixed.iter()), along=3)
    }
    
    output$ROC.Plot <- renderPlot({
      input$OSSButton
      pure = eval.ROC(est.str=reactives$acor.pure.pooled, true.str=reactives$true.str.pooled, plot.ROC=TRUE, show.AUC=FALSE, lightPDF=TRUE, lwd=2, col=cbPalette[1])
      mixed = eval.ROC(est.str=reactives$acor.mixed.pooled, true.str=reactives$true.str.pooled, plot.ROC=TRUE, show.AUC=FALSE, lightPDF=TRUE, lwd=2, col=cbPalette[2], add=TRUE)
      abline(a=0,b=1,lty=3,lwd=3)
      legend("bottomright", title='AUC', fill=c(cbPalette[1], cbPalette[2]), legend=c(paste('Pure cell type 1:', format(round(pure[['AUC']], 3))),paste('Mixed:', format(round(mixed[['AUC']],3)))))
      
    })
    
    Rplot <- function(){
      pure = eval.ROC(est.str=reactives$acor.pure.pooled, true.str=reactives$true.str.pooled, plot.ROC=TRUE, show.AUC=FALSE, lightPDF=TRUE, lwd=2, col=cbPalette[1])
      mixed = eval.ROC(est.str=reactives$acor.mixed.pooled, true.str=reactives$true.str.pooled, plot.ROC=TRUE, show.AUC=FALSE, lightPDF=TRUE, lwd=2, col=cbPalette[2], add=TRUE)
      abline(a=0,b=1,lty=3,lwd=3)
      legend("bottomright", title='AUC', fill=c(cbPalette[1], cbPalette[2]), legend=c(paste('Pure cell type 1:', format(round(pure[['AUC']], 3))),paste('Mixed:', format(round(mixed[['AUC']],3)))))
      
    }
    
    output$dld.ROC <- downloadHandler(
      filename = function(){
        paste(Sys.Date(),'.png', sep='')
      },
      content = function(file){
        png(file)
        print(Rplot())
        dev.off()
      }
    )
    
    })
}


shinyApp(ui, server)