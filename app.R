#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library( tidyverse )
library( tidygraph )
library( cowplot )
library( R.utils )
library( data.table )
library( visNetwork )
library( DT )
library( bslib )

dt_in <-
  fread(
    "data/9606.protein.physical.links.detailed.v12.0.txt.gz"
  )[combined_score >= 500, 
  ] %>%
  mutate(
    filt_me = 
      if_else(protein1 < protein2, 
              paste(protein1,"_",protein2, sep = ""),
              paste(protein2,"_",protein1, sep="") ))%>% 
  filter( !duplicated(filt_me))

info_dt <-
  fread(
    "data/9606.protein.info.v12.0.txt.gz"
  )

goi_fct <- fct_inorder(info_dt$preferred_name)

gois <- sort( info_dt$preferred_name )

base_graph <- tbl_graph(
  nodes = info_dt,
  edges = unique(dt_in),
  directed = FALSE
) 

# Define UI for application that draws a histogram
ui <- page_fluid(

  theme = bs_theme(bootswatch = "minty"),
    # Application title
    # titlePanel("STRING DB Multi-gene viewer"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          h2("STRING DB Multi-gene viewer"),
            selectizeInput("Genes",
                           "Select Gene(s) of Interest",
                           choices = NULL,
                           multiple = T),
            sliderInput("graph_size",
                        "Graph size (values above ~100 will load slowly)",
                        min = 1,
                        max = 250,
                        value = 50),
            radioButtons("score_cutoff",
                         "Minimum combined score",
                         choices = c(500,600,700,800,900),
                         selected = 500),
            submitButton( text = "Update Graph")
        ),

        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(
            tabPanel(
              "Network",
          visNetworkOutput("graph", width = "100%", height = "600px")
            ),
          tabPanel(
            "Node Table",
            dataTableOutput("node_DT")
            
          )
          )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

  updateSelectizeInput(session = session,
                       "Genes", choices = gois,
                       server = T)
  
  gene_graph <-
    eventReactive( {
      input$Genes
      input$graph_size
      input$score_cutoff
    },
                   {
      
                     req(input$Genes)

      # %>%
      #   activate(nodes) %>%
      #   mutate( nd1 = case_when(
      #     preferred_name %in% input$Genes ~ 0,
      #     node_is_adjacent( as.integer(
      #     goi_fct[as.character(goi_fct) %in% input$Genes])) ~ 1,
      #     T ~ 10),
      #     group = case_when( 
      #       nd1 == 0 ~ "GOI", 
      #       nd1 == 1 ~ "Other",
      #       T ~ "Exclude")
          # color.background = if_else( nd1 == 0, "red", "lightblue"),
          # color.border = "black",
          # color.hover = "yellow"
        # )
      
      # for( gene in input$Genes  ){
      #   
      #   temp_graph <-  temp_graph %>%
      #     mutate( 
      #       nd_temp =
      #               node_distance_to(
      #                 as.integer(
      #                   goi_fct[as.character(goi_fct) == gene])),
      #             nd1 = pmin( nd1, nd_temp, na.rm = T))
      #   
      # }
                     
      temp_graph <-
        base_graph %>%
        # activate(nodes) %>%
        # filter( nd1 <= 1  ) %>%
        activate(edges) %>% 
        mutate( counter = if_else( (.N()$preferred_name[to] %in% input$Genes | 
                                      .N()$preferred_name[from] %in% input$Genes), 1, 0),
                to_name = .N()$preferred_name[to],
                from_name = .N()$preferred_name[from],
                title = paste0("combined_score:", combined_score),
                value = combined_score/1000,
                color.opacity = value,
                physics = TRUE) %>%
        arrange( desc( counter ), desc( combined_score ) ) %>% 
        filter( !(row_number() > input$graph_size & counter == 1) &
                  !( to == from)) %>%
        activate(nodes) %>%
        mutate( 
          nd1 = case_when(
            preferred_name %in% input$Genes ~ 0,
            node_is_adjacent( as.integer(
              goi_fct[as.character(goi_fct) %in% input$Genes])) ~ 1,
            T ~ 10),
          group = case_when( 
            nd1 == 0 ~ "GOI", 
            nd1 == 1 ~ "Other",
            T ~ "Exclude"),
          connected = node_coreness(),
                id = preferred_name,
                label = preferred_name,
                title = preferred_name,
                physics = TRUE) %>%
        filter( group != "Exclude" &
                  connected >= 1 )
      
    })
  
  
    output$graph <- 
      renderVisNetwork({
        
        print("visnet")
        
        
        # gene_list <- dt_in %>% 
        #   filter(protein1 %in% input$Genes | protein2 %in% input$Genes) %>%
        #   select( protein1, protein2) %>%
        #   unlist() %>%
        #   unname()
        
        # temp_dt <- dt_in %>%
        #   # filter(protein1 %in% gene_list | protein2 %in% gene_list) %>%
        #   mutate( from = protein1, #as.integer(factor( protein1, levels = levels(goi_fct)) ), 
        #           to = protein2,#as.integer(factor( protein2, levels = levels(goi_fct))), 
        #           title = paste0("combined_score:", combined_score), 
        #           value = combined_score/1000, 
        #           physics = FALSE) %>%
        #   select( from, to, title, value, physics,
        #           protein1, protein2) #%>%
        #   
        # 
        # # print( info_dt )
        #   # head(n = 100)
        # 
        # temp_dt2 <- rbind( temp_dt %>% mutate(id = protein1), 
        #                    temp_dt%>% mutate(id = protein2)
        #                    ) %>%
        #   left_join( info_dt, by = c("id" = "#string_protein_id")) %>%
        #   # filter(id %in% gene_list) %>%
        #   mutate( label = preferred_name,
        #           title = preferred_name,
        #           physics = TRUE) %>%
        #   select( id, label, title, physics) %>%
        #   unique() %>%
        #   head(n = 75L )
        # 
        # print( temp_dt )
        # 
        # print( temp_dt2 )
          # mutate( id = as.integer(factor( id, levels = levels(goi_fct) ))) %>%
          # arrange( id )
        
        
        edge_dt <- gene_graph() %E>%
          as_tibble() %>%
          mutate( to = to_name,
                  from = from_name) %>%
          select( from, to, title, value )
        
        print(edge_dt %>% arrange( value ))
        
        visNetwork(as_tibble( gene_graph() ), edges = edge_dt,
                   height = "1000px", width = "100%") %>% 
          visNodes( font = list(size = "30")) %>%
          visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, collapse = TRUE) %>%
          visInteraction(navigationButtons = TRUE) %>%
          visGroups(groupname = "GOI", color = list(background = "red",
                                                            border = "black",
                                                            highlight =
                                                      list( background = "yellow",
                                                            border = "black"))) %>%
          visGroups(groupname = "Other", color = list(background = "lightblue",
                                                      border = "black",
                                                      highlight = list( background = "yellow",
                                                                        border = "black"))) %>%
          visLegend(width = 0.1, position = "right", main = "Group") %>%
          visPhysics(solver = "barnesHut", 
                     barnesHut = list(gravitationalConstant = -1500,
                                      avoidOverlap = 0.5,
                                      springConstant = 0.001,
                                      # springLengh = 200,
                                      # nodeDistance = 250,
                                      damping = 1),
                     stabilization = list(
                       iterations = 300
                     ), maxVelocity = 25) %>%
          visEdges(scaling = list( min = 3, max = 15),
                   color = list(
                     border = "grey50",
                     opacity = 0.50,
                     inherit = F,
                     highlight = "yellow")) %>%
          visExport(type = "png", 
                    name = paste0("network_",paste(input$Genes, sep = "_", collapse = "")), 
                    float = "left", 
                    label = "Save network", background = "white", style= "") 
        
        
      })
    
    output$node_DT <-
      DT::renderDataTable(
        datatable(
          gene_graph() %>%
            activate( nodes) %>%
            as_tibble() %>%
            mutate( Ensembl_ID =
                      paste0(
                        '<a href="https://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=',
                        str_remove(`#string_protein_id`,"9606."),
                        '" target="_blank"class="btn btn-primary">',
                        str_remove(`#string_protein_id`,"9606."),
                        '</a>')) %>%
            select( Ensembl_ID,
                    preferred_name,
                    protein_size,
                    annotation),
          filter = "top",
          escape = FALSE,
          extensions = "Buttons", 
          options = list(paging = TRUE,
                         scrollX=TRUE, 
                         searching = TRUE,
                         ordering = TRUE,
                         dom = 'Blfrtip',
                         buttons = c('copy', 'csv', 'excel'),
                         pageLength=5, 
                         lengthMenu=c(5,10,25) )
        ),
        escape = FALSE,
        server = FALSE
      )
      
}

# Run the application 
shinyApp(ui = ui, server = server)
