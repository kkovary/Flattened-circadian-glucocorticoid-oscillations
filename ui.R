source('functions/wat_ui.R')
source('functions/bat_ui.R')
source('functions/bat_kegg_ui.R')
source('functions/wat_kegg_ui.R')

ui <- navbarPage(
  title = 'In Vivo RNAseq',
  tabPanel(
    'Introduction',
    
    shiny::includeHTML(
      'HTML and Markdown/notebook.html'
    )
    
    
  ),
  navbarMenu('WAT Data',
             wat_ui,
             wat_kegg_ui
  ),
  navbarMenu('BAT Data',
             bat_ui,
             bat_kegg_ui
  )
)
