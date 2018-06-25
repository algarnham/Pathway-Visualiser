



library(shiny)
library(RColorBrewer)
library(limma)
library(edgeR)
library(matrixStats)
library(Glimma)
library(igraph)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(topGO)
library(network)
library(GGally)
library(threejs)
library(htmlwidgets)
library(visNetwork)
library(AnnotationDbi)
library(GO.db)
library(shiny)
library(DT)
library(plyr)

ui <- fluidPage(
   
   titlePanel("Pathway VisualiseR"),
   
   navlistPanel(
     tabPanel("Data upload", 
              p("Welcome blah blah blah"),
              p("Please upload your gene list below and specify the appropriate species and ontology. Currently Pathway VisualiseR supports human and mouse data."),
              fluidRow(
                column(width=5, fileInput(inputId = "genes", label = h3("Gene list upload"), buttonLabel = "Browse", placeholder = "Enter file name"),
                       radioButtons(inputId = "inputType", label=h3("Gene list type"), choices = list("Gene Symbol"="Symbol", "Entrez ID"="Entrez"))),
                column(width=2, radioButtons(inputId = "species", label = h3("Species"), choices = list("Human"="Hs", "Mouse"="Mm"))),
                column(width=3, radioButtons(inputId = "ont", label = h3("Ontology"), choices = list("Biological Process"= "BP", "Cellular Component"="CC", "Molecular Function"="MF")))
                ),
              fluidRow(
                actionButton(inputId = "Run", label = h4("Create network"))
              )
              ),
     tabPanel("GO network",
              sidebarLayout(
                sidebarPanel(
                  #textInput(inputId = "term_text", label=h4("Term search"), value="Enter term"),
                  #selectizeInput(inputId = "term_text", label=h4("Term search"), choices=as.list(as.character("nodes$label")), selected=character(0), options = list(placeholder="Enter term", maxItems=1)),
                  conditionalPanel(condition = "input.Run",
                                   uiOutput("nodes"),
                                   actionButton(inputId = "TermSearch", label = "Search"),
                                   uiOutput("genes"),
                                   actionButton(inputId = "GeneSearch", label= "Search")
                                   )
                ),
                
                mainPanel(
                  conditionalPanel(condition = "input.Run",
                                   visNetworkOutput("network"), 
                                   h2("Selected Term"),
                                   DTOutput("table"),
                                   br(),
                                   h2("Periphery Terms"),
                                   DTOutput("table2")
                  )
                  )
                )
              )
     , widths = c(2,10)
     )
)

server <- function(input, output, session) {
  
  #get gene info
  GeneInfo <- reactive({
    
    input$Run
    
    gene_type <- isolate(input$inputType)
    species <- isolate(input$species)
    #ont <- isolate(input$ont)
    
    genes <- read.table(input$genes$datapath, header=FALSE)[,1]
    
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
  
  #get GO results
  GO_terms_sig <- reactive({
    
    species <- isolate(input$species)
    ont <- isolate(input$ont)
    
    #get GO terms
    switch(species,
           Hs=assign("GO_terms", as.list(org.Hs.egGO2EG)),
           Mm=assign("GO_terms", as.list(org.Mm.egGO2EG))
    )
    GO_terms <- lapply(GO_terms, function(x) unique(x))
    
    #for gene list only run goana to identify significant GO terms
    GOres <- goana(de = GeneInfo()$ENTREZID, species = species)
    GOres <- GOres[GOres$P.DE<0.05,]
    GO_terms_sig <- GO_terms[names(GO_terms) %in% row.names(GOres)]
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
    
    #all terms need to have a min. of 10 genes to be included
    topGO_net <- new("topGOdata", ontology=ont, allGenes=GOI(), annot=annFUN.GO2genes, GO2genes=GO_terms_sig(), nodeSize=10)
    #topGO_graph <- graph(topGO_net)
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
    
    temp_nodes <- data.frame(id=temp_nodes, 
                            label=suppressMessages(select(GO.db, keys=temp_nodes, columns = "TERM")$TERM), 
                            title=suppressMessages(select(GO.db, keys=temp_nodes, columns = "TERM")$TERM), 
                            Term=suppressMessages(select(GO.db, keys=temp_nodes, columns = "TERM")$TERM), 
                            Ontology=suppressMessages(select(GO.db, keys=temp_nodes, columns="ONTOLOGY")$ONTOLOGY),
                            Definition=suppressMessages(select(GO.db, keys=temp_nodes, columns = "DEFINITION")$DEFINITION),
                            value=value,
                            GO_ID=GO_link)
    
    #genes with each term
    Genes <- ldply(genesInTerm(network_input()), function(x){paste0(x, collapse = ",")})
    # Genes <- lapply(Genes, function(x){switch(species,
    #                                           Hs=select(org.Hs.eg.db, keys=genes, keytype = "ENTREZID", columns = "SYMBOL"),
    #                                           Mm=select(org.Mm.eg.db, keys=genes, keytype = "ENTREZID", columns = "SYMBOL")
    #                                           )
    #   })
    
    temp_nodes <- data.frame(temp_nodes, Genes=Genes[match(temp_nodes$id, Genes[,1]), 2], row.names = NULL)
  })
  
  #network visualization
  output$network <- renderVisNetwork({
    visNetwork(nodes=net_nodes(), edges = net_edges()) %>%
      visIgraphLayout(randomSeed = 100, layout = "layout_nicely") %>%
      visEdges(smooth =FALSE, color = list(hover="red", highlight="purple")) %>%
      visInteraction(hover=TRUE, dragNodes = FALSE,  dragView = TRUE, hideEdgesOnDrag = FALSE, hoverConnectedEdges = TRUE, selectable = TRUE, selectConnectedEdges = TRUE, multiselect = TRUE) %>%
      visNodes(color=list(hover="red", highlight="purple"), chosen=TRUE) %>%
      visOptions(highlightNearest = list(enabled=TRUE, hover=FALSE), nodesIdSelection = list(enabled=FALSE)) %>%
      visEvents(select = "function(nodes){Shiny.onInputChange('current_node', nodes.nodes);;}")
  })
  
  #search list for GO terms(nodes)
  output$nodes <- renderUI({
    selectizeInput(inputId = "term_text", label=h4("Term search"), choices=as.list(as.character(net_nodes()$label)), selected=character(0), options = list(placeholder="Enter term", maxItems=1))
  })

  #search list for genes
  output$genes <- renderUI({
    selectizeInput(inputId = "gene_search", label=h4("Gene search"), choices=as.list(GeneInfo()$SYMBOL))
  })

  #event for clicking on term(node) in network
  observeEvent(input$current_node, {
    current_node <- input$current_node
    #visNetworkProxy("network") %>% visSelectNodes(id = current_node)
    output$table <- renderDataTable({net_nodes()[net_nodes()$id==current_node, c("GO_ID","Term","Definition")]}, escape=FALSE)
    output$table2 <- renderDataTable({net_nodes()[net_nodes()$id %in% unique(c(as.character(net_edges()$to[net_edges()$from==current_node]), as.character(net_edges()$from[net_edges()$to==current_node]))), c("GO_ID","Term","Definition")]}, escape=FALSE)
  })

  #event for serching for term
  observeEvent(input$TermSearch,  {
      isolate({
        current_node <- net_nodes()[grepl(input$term_text, net_nodes()$label), "id"]
        visNetworkProxy("network") %>% visSelectNodes(id = current_node)
        output$table <- renderDataTable({net_nodes()[net_nodes()$id %in% current_node, c("GO_ID","Term","Definition")]}, escape=FALSE)
        output$table2 <- renderDataTable({net_nodes()[net_nodes()$id %in% unique(c(as.character(net_edges()$to[net_edges()$from %in% current_node]), as.character(net_edges()$from[net_edges()$to %in% current_node]))), c("GO_ID","Term","Definition")]}, escape=FALSE)
      })
  })

  #event for searching for gene
  observeEvent(input$GeneSearch, {
    isolate({
      selected_gene <- GeneInfo()$ENTREZID[GeneInfo()$SYMBOL==input$gene_search]
      current_node <- net_nodes()[grepl(selected_gene, net_nodes()$Genes), "id"]
      visNetworkProxy("network") %>% visSelectNodes(id = current_node)
      output$table <- renderDataTable({net_nodes()[net_nodes()$id %in% current_node, c("GO_ID","Term","Definition")]}, escape=FALSE)
      #output$table2 <- renderDataTable()
      output$table2 <- renderDataTable({net_nodes()[net_nodes()$id %in% unique(c(as.character(net_edges()$to[net_edges()$from %in% current_node]), as.character(net_edges()$from[net_edges()$to %in% current_node]))), c("GO_ID","Term","Definition")]}, escape=FALSE)
    })
  })
  
  #event for hyperlink - selecting GO id in table
  #observeEvent(input$table_cell_clicked {
  #  cell <- input$table_cell_clicked
    
  #})
  
  


  
  
  
  ####when have data
  #create index
  #go_index <- ids2indices(gene.sets = GO_terms, identifiers = GeneInfo$ENTREZID)
  
  #run fry and identify significant GO terms
  #fry_res <- fry(y=GeneInfo$ENTREZID, index = go_index)
  #GO_terms_sig <- GO_terms[names(GO_terms) %in% row.names(fry_res)[fry_res$FDR.Mixed<0.05]]
  

}

shinyApp(ui = ui, server = server)

