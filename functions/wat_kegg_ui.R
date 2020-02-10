wat_kegg_ui <- tabPanel(
  'KEGG Maps',
  fluidRow(
    pickerInput(
      inputId = "wat_kegg_map",
      label = 'KEGG Pathway Search',
      selected = kegg_pathway_names[1],
      choices = kegg_pathway_names,
      options = list(`live-search` = TRUE,
                     title = 'KEGG Pathways')
    )
  ),
  fluidRow(
    plotOutput('wat_kegg_map') %>%
      withSpinner(type = 6)
    
  )
)