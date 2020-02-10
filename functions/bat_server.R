
#######################
# Input based filters #
#######################


##############################
# Clear gene table selection #
##############################

# When clicking on GO Terms table all rows are unselected
proxy_bat_gene_table = dataTableProxy('bat_gene_table')
observeEvent(input$bat_go_terms_table_rows_selected, {
  selectRows(proxy_bat_gene_table, NULL)
})

# All rows are unselected when KEGG row is selected
observeEvent(input$bat_kegg_table_rows_selected, {
  selectRows(proxy_bat_gene_table, NULL)
})

# When selecting points on plot all selected gene names are selected
observeEvent(input$bat_vol_plot_brush, {
  dims <- input$bat_vol_plot_brush
  genes <- filter(joined_bat_results,
                  -log10(lrt_qval) > dims$ymin, 
                  -log10(lrt_qval) < dims$ymax,
                  log2(fc) > dims$xmin,
                  log2(fc) < dims$xmax,
                  !duplicated(ext_gene))$ext_gene
  
  selectRows(proxy_bat_gene_table, which(joined_bat_results$ext_gene %in% genes))
})

##################################
# Clear GO Terms table selection #
##################################

# When clicking on gene table all rows are unselected
proxy_bat_go_terms_table = dataTableProxy('bat_go_terms_table')
observeEvent(input$bat_gene_table_rows_selected, {
  selectRows(proxy_bat_go_terms_table, NULL)
})

# All rows are unselected when KEGG row is selected
observeEvent(input$bat_kegg_table_rows_selected, {
  selectRows(proxy_bat_go_terms_table, NULL)
})

# When selecting points on plot all rows are unselected
observeEvent(input$bat_vol_plot_brush, {
  selectRows(proxy_bat_go_terms_table, NULL)
})

##############################
# Clear KEGG table selection #
##############################

# When clicking on gene table all rows are unselected
proxy_bat_kegg_table = dataTableProxy('bat_kegg_table')
observeEvent(input$bat_gene_table_rows_selected, {
  selectRows(proxy_bat_go_terms_table, NULL)
})

# When clicking on GO Terms table all rows are unselected
observeEvent(input$bat_go_terms_table_rows_selected, {
  selectRows(proxy_bat_kegg_table, NULL)
})

# When selecting points on plot all rows are unselected
observeEvent(input$bat_vol_plot_brush, {
  selectRows(proxy_bat_kegg_table, NULL)
})

##########################################
# Filter transcript data based on inputs #
##########################################
bat_highlighted_genes <- reactive({
  
  if(length(input$bat_gene_table_rows_selected) > 0){
    selected_genes <- joined_bat_results[input$bat_gene_table_rows_selected,]$ext_gene
    filter(joined_bat_results, ext_gene %in% selected_genes)
    
  } else if(length(input$bat_go_terms_table_rows_selected) > 0){
    selected_go_terms <- bat_go_terms[input$bat_go_terms_table_rows_selected,]$ID
    bat_go_terms %>% filter(ID %in% selected_go_terms) %>% 
      mutate(ENTREZID = strsplit(geneID,'/')) %>%
      dplyr::select(ONTOLOGY, ID, Description, ENTREZID) %>%
      unnest(cols = ENTREZID) %>% left_join(bat_universe, by = 'ENTREZID') %>%
      dplyr::rename(ext_gene = SYMBOL) %>% left_join(joined_bat_results, by = 'ext_gene') %>%
      filter(lrt_qval < 0.05)
  } else if(length(input$bat_kegg_table_rows_selected) > 0){
    
    selected_kegg <- bat_kegg[input$bat_kegg_table_rows_selected,]$Description
    
    left_join(joined_bat_results, kegg_genes, 'ext_gene') %>%
      filter(Description %in% selected_kegg, lrt_qval < 0.05)
  }
})

#######################
# Volcano Plot Output #
#######################

output$bat_vol_plot <- renderPlot({
  if(length(input$bat_gene_table_rows_selected) > 0){
    
    vol_plot_base(joined_bat_results, 'BAT') +
      
      geom_point(alpha = 0.8, 
                 size = 4, 
                 stroke = 0, 
                 data = bat_highlighted_genes(),
                 aes(color = ext_gene)) +
      geom_label_repel(data = bat_highlighted_genes(),
                       aes(label = external_transcript_name, 
                           color = ext_gene), 
                       size = 6, 
                       force = 20) +
      labs(color = "Gene Name")
  } else if(length(input$bat_go_terms_table_rows_selected) > 0){
    vol_plot_base(joined_bat_results, 'BAT') +
      geom_point(alpha = 1, 
                 size = 3, 
                 stroke = 0, 
                 data = bat_highlighted_genes(),
                 aes(color = Description)) +
      labs(color = "Gene Name")
  } else if(length(input$bat_kegg_table_rows_selected) > 0){
    vol_plot_base(joined_bat_results, 'BAT') +
      geom_point(alpha = 1, 
                 size = 3, 
                 stroke = 0, 
                 data = bat_highlighted_genes(),
                 aes(color = Description)) +
      labs(color = "Gene Name")
  } else{
    vol_plot_base(joined_bat_results, 'BAT')
  }
})

########################
# Gene Bar Plot Output #
########################
output$bat_gene_bar_plot <- renderPlot({
  if(length(input$bat_gene_table_rows_selected) > 0){
    
    gene_bp_gn(joined_bat_gene_results, bat_highlighted_genes(), 'BAT')
    
  } else if(length(input$bat_go_terms_table_rows_selected) > 0){
    
    gene_bp_go(bat_go_terms, 
               input$bat_go_terms_table_rows_selected,
               bat_universe,
               joined_bat_gene_results,
               'BAT',
               'go')
  } else if(length(input$bat_kegg_table_rows_selected) > 0){
    
    gene_bp_go(bat_kegg, 
               input$bat_kegg_table_rows_selected,
               kegg_genes,
               joined_bat_gene_results,
               'BAT',
               'kegg')
  } else{
    NULL
  }
})

##############################
# Transcript Bar Plot Output #
##############################
output$bat_transcript_bar_plot <- renderPlot({
  if(length(input$bat_gene_table_rows_selected) > 0){
    
    tran_bp_gn(bat_highlighted_genes(),'BAT')
    
  } else if(length(input$bat_go_terms_table_rows_selected) > 0){
    
    tran_bp_go(bat_highlighted_genes(),'BAT')
    
  } else if(length(input$bat_kegg_table_rows_selected) > 0){
    
    tran_bp_go(bat_highlighted_genes(),'BAT')
    
  } else{
    NULL
  }
})





#####################
# Gene Table Output #
#####################
output$bat_gene_table <- DT::renderDataTable({
  DT::datatable(joined_bat_results,
                selection = 'multiple',
                options = list(buttons = c('csv', 'excel'))
  )
})

########################
# GO Term Table Output #
########################
output$bat_go_terms_table <- DT::renderDataTable({
  DT::datatable(dplyr::select(bat_go_terms,
                              ONTOLOGY,
                              ID,
                              Description,
                              GeneRatio,
                              Count,
                              p.adjust,
                              zscore),
                selection = 'multiple',
                options = list(buttons = c('csv', 'excel'))
  )
})

#####################
# KEGG Table Output #
#####################
output$bat_kegg_table <- DT::renderDataTable({
  
  DT::datatable(bat_kegg,
                selection = 'multiple',
                options = list(buttons = c('csv', 'excel'))
  )
})