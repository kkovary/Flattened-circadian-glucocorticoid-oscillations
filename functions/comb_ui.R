comb_ui <- tabPanel(
  'RNA-seq Plots',
  tabsetPanel(
    tabPanel(
      'Volcano Plot',
      plotOutput('comb_vol_plot',
                 brush = "comb_vol_plot_brush",
                 dblclick = 'comb_vol_plot_click',
                 height = "600px") %>%
        withSpinner(type = 6)
    ),
    tabPanel(
      'Transcript Bar Plot',
      plotOutput('comb_transcript_bar_plot',
                 height = "600px") %>%
        withSpinner(type = 6)
    ),
    tabPanel(
      'Gene Bar Plot',
      plotOutput('comb_gene_bar_plot',
                 height = "600px") %>%
        withSpinner(type = 6)
    )
  ),
  tabsetPanel(
    id = 'comb_tables',
    tabPanel(
      'Gene Table',
      DTOutput('comb_gene_table') %>% withSpinner(type = 6)
    ),
    tabPanel('GO Term Table',
             DTOutput('comb_go_terms_table')
    ),
    tabPanel('KEGG Table',
             DTOutput('comb_kegg_table')
    )
  )
)