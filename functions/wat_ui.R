wat_ui <- tabPanel(
  'WAT Data',
  tabsetPanel(
    tabPanel(
      'Volcano Plot',
      plotOutput('wat_vol_plot',
                 brush = "wat_vol_plot_brush",
                 dblclick = 'wat_vol_plot_click',
                 height = "600px") %>%
        withSpinner(type = 6)
    ),
    tabPanel(
      'Gene Bar Plot',
      plotOutput('wat_gene_bar_plot',
                 height = "600px") %>%
        withSpinner(type = 6)
    ),
    tabPanel(
      'Transcript Bar Plot',
      plotOutput('wat_transcript_bar_plot',
                 height = "600px") %>%
        withSpinner(type = 6)
    )
  ),
  tabsetPanel(
    id = 'wat_tables',
    tabPanel(
      'Gene Table',
      DTOutput('wat_gene_table') %>% withSpinner(type = 6)
    ),
    tabPanel('GO Term Table',
             DTOutput('wat_go_terms_table'))
  )
)