## app.R ##

# flexdashboard implementation of GOvisualiser
# Alexandra Garnham and Goknur Giner, The Walter and Eliza Hall Institute of Medical Research
# Created: 18th June 2018 Modified: June 26th, 2018
#---------------------------------------------------------------------------
library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(RColorBrewer)
library(limma)
library(edgeR)
#library(matrixStats)
#library(Glimma)
library(igraph)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(topGO)
#library(network)
#library(GGally)
#library(threejs)
library(htmlwidgets)
library(visNetwork)
library(AnnotationDbi)
library(GO.db)
#library(shiny)
library(DT)
library(plyr)

options(shiny.maxRequestSize = 100*1024^2) # increase the size of the file can be uploaded

header <- dashboardHeader(
  #title = "GOvisualiser",
  title = tags$a(tags$img(src="logoname-02.png", height=36, width=206)),
  titleWidth = 300
)

sidebar <- ## Sidebar content
  dashboardSidebar(
    sidebarMenu(id = "tabs",
                menuItem("UPLOAD DATA", icon = icon("upload"), tabName = "upload"),
                menuItem("GO Network", icon = icon("database"), tabName = "network")
    ),
    width = 300
  )

body <- dashboardBody(
  tabItems(
    tabItem(
      tabName = "upload",
      box(title = "Welcome to Pathway-Visualiser", width = 12, solidHeader = TRUE, status = "primary",
          p("Pathway-Visualiser is an interactive website that allows users to visualise network structure within ", a(href='http://www.geneontology.org/', 'Gene Ontology (GO)'),  "database based on their gene input.")),
      box(
        title = "UPLOAD DATA", width = 12, solidHeader = TRUE, status = "primary",
        selectInput("typeofinput", 
                    h3("Select input type"), choices = list("List of genes", "Matrix of counts", "Fit object"), selected = "List of genes"),
        conditionalPanel(condition = "input.typeofinput == 'List of genes'",
                         fileInput(inputId = "genes", label = h3("Gene list upload"), buttonLabel = "Browse", placeholder = "Enter file name"))
        # conditionalPanel(condition = "input.typeofinput == 'Matrix of counts'",
        #                  fileInput(inputId = "counts", label = h3("Upload a matrix of counts"))),
        # conditionalPanel(condition = "input.typeofinput == 'Fit object'",
        #                  fileInput("fit", label = h3("Upload a fit object")))
      ),
      box(
        title = "DEFAULT SETTINGS", width = 12, solidHeader = TRUE, status = "primary",
        fluidRow(
          column(width = 2,
                 radioButtons(inputId = "inputType", label = h3("Gene list type"), choices = list("Gene Symbol" = "Symbol", "Entrez ID" = "Entrez"))),
          column(width = 2, radioButtons(inputId = "species", label = h3("Species"), choices = list("Human" = "Hs", "Mouse" = "Mm"))),
          column(width = 2, radioButtons(inputId = "ont", label = h3("Ontology"), choices = list("Biological Process" = "BP", "Cellular Component" = "CC", "Molecular Function" = "MF")))
        ),
        fluidRow(
          column(width = 12,
                 actionButton(inputId = "Run", label = h4("Create network")),
                 conditionalPanel(condition = "input.Run",
                                  h4("Please proceed to the GO Network tab")
                                  )
          )
        )
      )
    ),
    tabItem(
      tabName = "network",
      box(title = paste("Pathway-Visualiser results"), width = 12, solidHeader = TRUE, status = "primary",
          fluidRow(
            conditionalPanel(condition = "input.Run",
                             column(width=3, uiOutput("nodes")),
                             column(width=3, uiOutput("genes"))
            )
            ),
          fluidRow(
            conditionalPanel(condition = "input.Run",
                             column(width=3, actionButton(inputId = "TermSearch", label = "Search")),
                             column(width=3, actionButton(inputId = "GeneSearch", label= "Search"))
            )
          ),
          fluidRow(
            conditionalPanel(condition = "input.Run",
                             column(width=12,
                                    visNetworkOutput("network") %>% withSpinner()
                             ))
          )
      ),
      box(title="GO terms", width=8, solidHeader = TRUE, status = "primary",
          conditionalPanel(condition = "input.Run",
                           h3(textOutput("Ontology")),
                           h3("Selected Term"),
                           DT::dataTableOutput("table"),
                           br(),
                           h3("Periphery Terms"),
                           DTOutput("table2")
                           )
      ),
      box(title="Genes", width=4, solidHeader = TRUE, status = "primary",
          conditionalPanel(condition="input.Run",
                           h3("GO term"),
                           textOutput("View_selectedTerm"),
                           h4("Genes for this term in your data"),
                           textOutput("View_genesInput"),
                           h4("All other genes for this term"),
                           textOutput("View_genesNotInput")
          )
      )
    )
  )
)


server <- function(input, output, session) {
  
  shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }
  
  output$Ontology <- renderText(switch(isolate(input$ont),
                                       BP="Selected ontology: Biological Process",
                                       CC="Selected ontology: Cellular Component",
                                       MF="Selected ontology: Molecular Function"))
  
  #get gene info
  GeneInfo <- reactive({
    
    input$Run
    
    gene_type <- isolate(input$inputType)
    species <- isolate(input$species)
    #ont <- isolate(input$ont)
    
    genes <- read.table(input$genes$datapath, header=FALSE, stringsAsFactors = FALSE)[,1]
    
    #identify both gene symbol and entrez id
    if(gene_type=="Symbol"){
      genes <- alias2Symbol(alias = genes, species = species)
      switch(species,
             Hs=assign("GeneInfo", select(org.Hs.eg.db, keys=genes, keytype="SYMBOL", columns = "ENTREZID")),
             Mm=assign("GeneInfo", select(org.Mm.eg.db, keys=genes, keytype = "SYMBOL", columns="ENTREZID"))
      )
    }else{
      switch(species,
             Hs=assign("GeneInfo", select(org.Hs.eg.db, keys=genes, keytype = "ENTREZID", columns = "SYMBOL")),
             Mm=assign("GeneInfo", select(org.Mm.eg.db, keys=genes, keytype = "ENTREZID", columns = "SYMBOL"))
      )
    }
  })
  
  GO_terms <- reactive({
    
    species <- isolate(input$species)
    
    #get GO terms
    switch(species,
           Hs=assign("GO_terms", as.list(org.Hs.egGO2ALLEGS)),
           Mm=assign("GO_terms", as.list(org.Mm.egGO2ALLEGS))
    )
    GO_terms <- lapply(GO_terms, function(x) unique(x))
  })
  
  # GO_terms_all <- reactive({
  #   species <- isolate(input$species)
  #   switch(species,
  #          Hs=assign("GO_terms_all", as.list(org.Hs.egGO2ALLEGS)),
  #          Mm=assign("GO_terms_all", as.list(org.Mm.egGO2ALLEGS))
  #   )
  #   GO_terms_all <- lapply(GO_terms_all, function(x) unique(x))
  # })
  
  #get GO results
  GOres <- reactive({
    
    species <- isolate(input$species)
    
    #for gene list only, run goana to identify significant GO terms
    GOres <- goana(de = GeneInfo()$ENTREZID, species = species)
  })
  
  GO_terms_sig <- reactive({
   GO_terms_sig <- GO_terms()[names(GO_terms()) %in% row.names(GOres()[GOres()$P.DE<0.05,])]
   #terms_10genes <- sapply(GO_terms_sig, function(x) length(x)>=10)
   #GO_terms_sig <- GO_terms_sig[terms_10genes]
  })
  
  #create factor designating genes of interest - required by topGO
  GOI <- reactive({
    temp_GOI <- factor(GeneInfo()$ENTREZID %in% GeneInfo()$ENTREZID, levels = c(FALSE,TRUE))
    levels(temp_GOI) <- c(0,1)
    names(temp_GOI) <- GeneInfo()$ENTREZID
    temp_GOI <- temp_GOI
  })
  
  #create network
  network_input <- reactive({  
    
    ont <- isolate(input$ont)
    
    #all terms need to have a min. of 10 genes to be included - may have branchs of non-significant terms due to node size!
    #may be some significant GO terms that have small node size and are therefore removed, however the non-significant conecting ones may stay
    topGO_net <- new("topGOdata", ontology=ont, allGenes=GOI(), annot=annFUN.GO2genes, GO2genes=GO_terms_sig(), nodeSize=10)
  })
  
  #get network edges
  net_edges <- reactive({
    temp_edges <- edges(graph(network_input()))
    temp_edges <- setNames(unlist(temp_edges, use.names = FALSE), rep(names(temp_edges), lengths(temp_edges)))
    temp_edges <- data.frame(from=names(temp_edges), to=temp_edges)
  })
  
  #get network nodes
  net_nodes <- reactive({
    
    #node size
    temp_nodes <- nodes(graph(network_input()))
    value <- rep(10, length(temp_nodes))
    names(value) <- temp_nodes
    value[names(value) %in% names(GO_terms_sig())] <- 25
    
    #hyperlinks
    GO_link <- paste0("<a href='http://amigo.geneontology.org/amigo/term/", temp_nodes, "'> ", temp_nodes, "</a>")
    #GO_gene_link <- sprintf('<input type="radio" name="%s" value="%s"/>',
    #                        temp_nodes, temp_nodes)
    GO_genes <- shinyInput(actionButton, length(temp_nodes), "button_", label = "View", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)')
    
    #genes with each term
    Genes <- genesInTerm(network_input())
    NGenes <- sapply(Genes, function(x) length(x))
    
    temp_nodes <- data.frame(id=temp_nodes, 
                             label=suppressMessages(select(GO.db, keys=temp_nodes, columns = "TERM")$TERM), 
                             title=suppressMessages(select(GO.db, keys=temp_nodes, columns = "TERM")$TERM), 
                             Term=sapply(suppressMessages(select(GO.db, keys=temp_nodes, columns = "TERM")$TERM), function(x) gsub("_", " ", x)), 
                             Ontology=suppressMessages(select(GO.db, keys=temp_nodes, columns="ONTOLOGY")$ONTOLOGY),
                             Definition=suppressMessages(select(GO.db, keys=temp_nodes, columns = "DEFINITION")$DEFINITION),
                             value=value,
                             "GO ID"=GO_link,
                             PValue=round(GOres()$P.DE[match(temp_nodes, row.names(GOres()))], digits = 3),
                             "View Genes"=GO_genes, 
                             "Number Genes"=NGenes, stringsAsFactors = FALSE, check.names = FALSE)
    
    
    #Genes <- lapply(Genes, function(x){switch(species,
    #                                           Hs=suppressMessages(select(org.Hs.eg.db, keys=x, keytype = "ENTREZID", columns = "SYMBOL")[,"SYMBOL"]),
    #                                           Mm=suppressMessages(select(org.Mm.eg.db, keys=x, keytype = "ENTREZID", columns = "SYMBOL")[,"SYMBOL"])
    #                                           )
    #   })
    Genes <- ldply(genesInTerm(network_input()), function(x){paste0(x, collapse = ", ")})
    
    temp_nodes <- data.frame(temp_nodes, Genes=Genes[match(temp_nodes$id, Genes[,1]), 2], row.names = NULL, check.names = FALSE)
  })
  
  #network visualization
  output$network <- renderVisNetwork({
    visNetwork(nodes=net_nodes(), edges = net_edges()) %>%
      visIgraphLayout(randomSeed = 100, layout = "layout_nicely") %>%
      visEdges(smooth =FALSE, color = list(hover="red", highlight="purple"), physics = FALSE) %>%
      visInteraction(hover=TRUE, dragNodes = FALSE,  dragView = TRUE, hideEdgesOnDrag = FALSE, hoverConnectedEdges = TRUE, selectable = TRUE, selectConnectedEdges = TRUE, multiselect = TRUE) %>%
      visNodes(color=list(hover="red", highlight="purple"), chosen=TRUE, physics = FALSE) %>%
      visOptions(highlightNearest = list(enabled=TRUE, hover=FALSE), nodesIdSelection = list(enabled=FALSE)) %>%
      visEvents(select = "function(nodes){Shiny.onInputChange('current_node', nodes.nodes);;}") %>%
      visLegend(addNodes = list(list(label="Significant \nterm", size=5, shape="dot", text.align="left", font.size=10),
                                list(label="Non-significant \nterm", size=1.5, shape="dot", font.align="top"),
                                list(label="Selected \nterm", color="purple", size=5, shape="dot", font.align="top")),
                useGroups = FALSE, zoom=FALSE) %>%
      visExport(type="pdf")
  })
  
  #search list for GO terms(nodes)
  output$nodes <- renderUI({
    selectizeInput(inputId = "term_text", label=h4("Term search"), choices=as.list(as.character(net_nodes()$label)), selected=character(0), options = list(placeholder="Enter term", maxItems=1))
  })
  
  #search list for genes
  output$genes <- renderUI({
    selectizeInput(inputId = "gene_search", label=h4("Gene search"), choices=as.list(GeneInfo()$SYMBOL), options=list(placeholder="Enter gene"))
  })
  
  #event for clicking on term(node) in network
  observeEvent(input$current_node, {
    current_node <- input$current_node
    #visNetworkProxy("network") %>% visSelectNodes(id = current_node)
    output$table <- renderDataTable({net_nodes()[net_nodes()$id==current_node, c("GO ID", "Term","Definition", "PValue", "Number Genes", "View Genes")]}, 
                                    escape=FALSE,  selection="none")
    output$table2 <- renderDataTable({net_nodes()[net_nodes()$id %in% unique(c(as.character(net_edges()$to[net_edges()$from==current_node]), as.character(net_edges()$from[net_edges()$to==current_node]))), 
                                                  c("GO ID", "Term","Definition", "PValue", "Number Genes", "View Genes")]}, escape=FALSE, selection="none")
  })
  
  #event for serching for term
  observeEvent(input$TermSearch,  {
    isolate({
      current_node <- net_nodes()[grepl(input$term_text, net_nodes()$label), "id"]
      visNetworkProxy("network") %>% visSelectNodes(id = current_node)
      output$table <- renderDataTable({net_nodes()[net_nodes()$id %in% current_node, c("GO ID", "Term","Definition", "PValue", "Number Genes", "View Genes")]}, 
                                      escape=FALSE, selection="none")
      output$table2 <- renderDataTable({net_nodes()[net_nodes()$id %in% unique(c(as.character(net_edges()$to[net_edges()$from %in% current_node]), as.character(net_edges()$from[net_edges()$to %in% current_node]))), 
                                                    c("GO ID", "Term","Definition", "PValue", "Number Genes", "View Genes")]}, escape=FALSE, selection="none")
    })
  })
  
  #event for searching for gene
  observeEvent(input$GeneSearch, {
    isolate({
      selected_gene <- GeneInfo()$ENTREZID[GeneInfo()$SYMBOL==input$gene_search]
      current_node <- net_nodes()[grepl(selected_gene, net_nodes()$Genes), "id"]
      visNetworkProxy("network") %>% visSelectNodes(id = current_node)
      output$table <- renderDataTable({net_nodes()[net_nodes()$id %in% current_node, c("GO ID", "Term","Definition", "PValue", "Number Genes", "View Genes")]}, 
                                      escape=FALSE, selection="none")
      output$table2 <- renderDataTable({net_nodes()[net_nodes()$id %in% unique(c(as.character(net_edges()$to[net_edges()$from %in% current_node]), as.character(net_edges()$from[net_edges()$to %in% current_node]))), 
                                                    c("GO ID", "Term","Definition", "PValue", "Number Genes", "View Genes")]}, escape=FALSE, selection="none")
    })
  })
  
  observeEvent(input$select_button, {

    selectedRow <- as.numeric(strsplit2(input$select_button, "_")[1,2])
    term_selected <- net_nodes()$id[selectedRow]
    species <- isolate(input$species)
    gene_list <- GO_terms()[[term_selected]]
    gene_list <- switch(species,
                       Hs=suppressMessages(select(org.Hs.eg.db, keys=gene_list, keytype = "ENTREZID", columns = "SYMBOL")[,"SYMBOL"]),
                       Mm=suppressMessages(select(org.Mm.eg.db, keys=gene_list, keytype = "ENTREZID", columns = "SYMBOL")[,"SYMBOL"])
                       )
    genes_Input <- gene_list[gene_list %in% GeneInfo()$SYMBOL]
    genes_Input <- paste0(genes_Input[order(genes_Input)], collapse = ", ")
    genes_NotInput <- gene_list[!(gene_list %in% GeneInfo()$SYMBOL)]
    genes_NotInput <- paste0(genes_NotInput[order(genes_NotInput)], collapse = ", ")
    output$View_selectedTerm <- renderText(term_selected)
    output$View_genesInput <- renderText(genes_Input)
    output$View_genesNotInput <- renderText(genes_NotInput)
  })
  
  
  
  ####when have data
  #create index
  #go_index <- ids2indices(gene.sets = GO_terms, identifiers = GeneInfo$ENTREZID)
  
  #run fry and identify significant GO terms
  #fry_res <- fry(y=GeneInfo$ENTREZID, index = go_index)
  #GO_terms_sig <- GO_terms[names(GO_terms) %in% row.names(fry_res)[fry_res$FDR.Mixed<0.05]]
  
  
}
ui <- dashboardPage(header, sidebar, body, skin = "blue")
shinyApp(ui, server)
