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
library(GO.db)
library(shiny)
library(DT)

ui <- fluidPage(
   
   titlePanel("Pathway VisualiseR"),
   
   navlistPanel(
     tabPanel("Data upload", 
              p("Welcome blah blah blah"),
              p("Please upload your gene list below and specify the appropriate species and ontology. Currently Pathway VisualiseR supports human and mouse data."),
              fluidRow(
                column(width=5, fileInput(inputId = "genes", label = h3("Gene list upload"), buttonLabel = "Browse", placeholder = "Enter file name"),
                       checkboxGroupInput(inputId = "inputType", label=h3("Gene list type"), choices = list("Gene Symbol"="Symbol", "Entrez ID"="Entrez"))),
                column(width=2, checkboxGroupInput(inputId = "species", label = h3("Species"), choices = list("Human"="Hs", "Mouse"="Mm"))),
                column(width=3, checkboxGroupInput(inputId = "ont", label = h3("Ontology"), choices = list("Biological Process"= "BP", "Cellular Component"="CC", "Molecular Function"="MF")))
                ),
              fluidRow(
                actionButton(inputId = "Run", label = h4("Create network"))
              )
              ),
     tabPanel("GO network",
              sidebarLayout(
                sidebarPanel(
                  conditionalPanel(condition = "input.Run",
                    uiOutput("terms")
                    #textInput(uiOutput("goterms")),
                    )
                  # DTOutput("goterms"),
                  # verbatimTextOutput("terms")
                  
                  
                    #textInput(inputId = "", label=h4("Term search"), value="Enter term")
                  #selectizeInput(inputId = "term_text", label=h4("Term search"), choices=as.list(as.character(net_nodes$label)), selected=character(0), options = list(placeholder="Enter term", maxItems=1))
                ),
                
                mainPanel(
                  conditionalPanel(condition = "input.Run",
                  visNetworkOutput("network")), 
                  #verbatimTextOutput("shiny_return")
                  h2("Selected Term"),
                  DTOutput("table"),
                  #dataTableOutput("text")
                  br(),
                  h2("Perefery Terms"),
                  DTOutput("table2")
                  )
                )
              )
     , widths = c(2,10))
)

server <- function(input, output, session) {
  network_input <- reactive({
    
    input$Run
    
    gene_type <- isolate(input$inputType)
    species <- isolate(input$species)
    ont <- isolate(input$ont)
    
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
    
    #create factor designating genes of interest - required by topGO
    GOI <- factor(GeneInfo$ENTREZID %in% GeneInfo$ENTREZID, levels = c(FALSE,TRUE))
    levels(GOI) <- c(0,1)
    names(GOI) <- GeneInfo$ENTREZID
    
    #get GO terms
    switch(species,
           Hs=assign("GO_terms", as.list(org.Hs.egGO2EG)),
           Mm=assign("GO_terms", as.list(org.Mm.egGO2EG))
    )
    GO_terms <- lapply(GO_terms, function(x) unique(x))
    
    #for gene list only run goana to identify significant GO terms
    GOres <- goana(de = GeneInfo$ENTREZID, species = species)
    GOres <- GOres[GOres$P.DE<0.05,]
    GO_terms_sig <- GO_terms[names(GO_terms) %in% row.names(GOres)]
    
    #create network
    #all terms need to have a min. of 10 genes to be included
    topGO_net <- new("topGOdata", ontology=ont, allGenes=GOI, annotationFun=annFUN.GO2genes, GO2genes=GO_terms_sig, nodeSize=10)
    topGO_graph <- graph(topGO_net)
  })
  
  #get network edges
  net_edges <- reactive({
    temp_edges <- edges(network_input())
    temp_edges <- setNames(unlist(temp_edges, use.names = FALSE), rep(names(temp_edges), lengths(temp_edges)))
    temp_edges <- data.frame(from=names(temp_edges), to=temp_edges)
  })
  
  #get network nodes
  net_nodes <- reactive({
    temp_nodes <- nodes(network_input())
    temp_nodes <- data.frame(id=temp_nodes, 
                            label=suppressMessages(select(GO.db, keys=temp_nodes, columns = "TERM")$TERM), 
                            title=suppressMessages(select(GO.db, keys=temp_nodes, columns = "TERM")$TERM), 
                            Term=suppressMessages(select(GO.db, keys=temp_nodes, columns = "TERM")$TERM), 
                            Ontology=suppressMessages(select(GO.db, keys=temp_nodes, columns="ONTOLOGY")$ONTOLOGY),
                            Definition=suppressMessages(select(GO.db, keys=temp_nodes, columns = "DEFINITION")$DEFINITION))
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
  
  
  
  observeEvent(input$current_node, {
    current_node <- input$current_node
    output$table <- renderDataTable({net_nodes()[net_nodes()$id==current_node, c("id","Term","Definition")]})
    output$table2 <- renderDataTable({net_nodes()[net_nodes()$id %in% unique(c(as.character(net_edges()$to[net_edges()$from==current_node]), as.character(net_edges()$from[net_edges()$to==current_node]))), c("id","Term","Definition")]})
    # output$goterms <- renderText({as.list(as.character(net_nodes()$label))})
    # output$goterms <- renderUI({current_node})
    output$goterms <- renderDataTable({net_nodes()[,1:2]})
    # output$terms <- renderPrint({
    #     print(as.list(as.character(net_nodes()$label)))
    #   })
    output$terms <- renderUI({
     selectInput("term_text", label = h4("Term search"),
       choices = as.list(as.character(net_nodes()$label)), multiple = TRUE)
     })
  })
  
  # observe(input$term_text, {
  #   isolate({
  #     curent_node <- net_nodes()[grepl(input$term_text, net_nodes()$label), "id"]
  #     visNetworkProxy("network") %>% visSelectNodes(id = curent_node)
  #     output$table <- renderDataTable({net_nodes()[net_nodes()$id %in% curent_node, c("id","Term","Definition")]})
  #     output$table2 <- renderDataTable({net_nodes()[net_nodes()$id %in% unique(c(as.character(net_edges()$to[net_edges()$from %in% curent_node]), as.character(net_edges()$from[net_edges()$to %in% curent_node]))), c("id","Term","Definition")]}, width="50%")
  #   })
  # })

  
  
  ####when have data
  #create index
  #go_index <- ids2indices(gene.sets = GO_terms, identifiers = GeneInfo$ENTREZID)
  
  #run fry and identify significant GO terms
  #fry_res <- fry(y=GeneInfo$ENTREZID, index = go_index)
  #GO_terms_sig <- GO_terms[names(GO_terms) %in% row.names(fry_res)[fry_res$FDR.Mixed<0.05]]
  

}

shinyApp(ui = ui, server = server)

