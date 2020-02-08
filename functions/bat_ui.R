bat_ui <- tabPanel(
  'BAT Data',
  tabsetPanel(
    tabPanel(
      'Volcano Plot',
      plotOutput('bat_vol_plot',
                 brush = "bat_vol_plot_brush",
                 dblclick = 'bat_vol_plot_click',
                 height = "600px") %>%
        withSpinner(type = 6)
    ),
    tabPanel(
      'Transcript Bar Plot',
      plotOutput('bat_transcript_bar_plot',
                 height = "600px") %>%
        withSpinner(type = 6)
    ),
    tabPanel(
      'Gene Bar Plot',
      plotOutput('bat_gene_bar_plot',
                 height = "600px") %>%
        withSpinner(type = 6)
    )
  ),
  tabsetPanel(
    id = 'bat_tables',
    tabPanel(
      'Gene Table',
      DTOutput('bat_gene_table') %>% withSpinner(type = 6)
    ),
    tabPanel('GO Term Table',
             DTOutput('bat_go_terms_table'))
  )
)