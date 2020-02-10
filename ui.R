source('functions/wat_ui.R')
source('functions/bat_ui.R')

ui <- navbarPage(
  title = 'In Vivo RNAseq',
  tabPanel(
    'Introduction',
    
    shiny::includeHTML(
      'HTML and Markdown/notebook.html'
    )

    
  ),
  wat_ui,
  bat_ui
)
