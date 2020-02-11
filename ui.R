source('functions/wat_ui.R')
source('functions/bat_ui.R')
source('functions/bat_kegg_ui.R')
source('functions/wat_kegg_ui.R')

ui <- navbarPage(
  title = 'In Vivo RNAseq',
  tabPanel(
    'Introduction',
    
    # shiny::includeHTML(
    #   'HTML and Markdown/notebook.html'
    # )
    
    tags$iframe(src = 'notebook.html', # put myMarkdown.html to /www
                width = '100%', height = '800px', 
                frameborder = 0, scrolling = 'auto'
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
