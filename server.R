shinyServer(function(input, output, session) {
  
  source('functions/bat_server.R', local = TRUE)
  source('functions/wat_server.R', local = TRUE)
  source('functions/bat_kegg_server.R', local = TRUE)
  source('functions/wat_kegg_server.R', local = TRUE)
  
})