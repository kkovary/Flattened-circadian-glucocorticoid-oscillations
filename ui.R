source('functions/wat_ui.R')
source('functions/bat_ui.R')

ui <- navbarPage(
  title = 'In Vivo RNAseq',
  tabPanel(
    'Introduction',
    shiny::includeHTML(
      'Flattened-circadian-glucocorticoid-oscillations-cause-obesity-due-to-increased.html'
    )
  ),
  wat_ui,
  bat_ui
)
