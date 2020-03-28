comb_kegg_ui <- tabPanel(
  'KEGG Maps',
  fluidRow(
    pickerInput(
      inputId = "comb_kegg_map",
      label = 'KEGG Pathway Search',
      selected = kegg_pathway_names[1],
      choices = kegg_pathway_names,
      options = list(`live-search` = TRUE,
                     title = 'KEGG Pathways')
    )
  ),
  fluidRow(
    plotOutput('comb_kegg_map') %>%
      withSpinner(type = 6)
    
  )
)