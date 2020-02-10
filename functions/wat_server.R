
#######################
# Input based filters #
#######################


##############################
# Clear gene table selection #
##############################

# When clicking on GO Terms table all rows are unselected
proxy_wat_gene_table = dataTableProxy('wat_gene_table')
observeEvent(input$wat_go_terms_table_rows_selected, {
  selectRows(proxy_wat_gene_table, NULL)
})

# All rows are unselected when KEGG row is selected
observeEvent(input$wat_kegg_table_rows_selected, {
  selectRows(proxy_wat_gene_table, NULL)
})

# When selecting points on plot all selected gene names are selected
observeEvent(input$wat_vol_plot_brush, {
  dims <- input$wat_vol_plot_brush
  genes <- filter(joined_wat_results,
                  -log10(lrt_qval) > dims$ymin, 
                  -log10(lrt_qval) < dims$ymax,
                  log2(fc) > dims$xmin,
                  log2(fc) < dims$xmax,
                  !duplicated(ext_gene))$ext_gene
  
  selectRows(proxy_wat_gene_table, which(joined_wat_results$ext_gene %in% genes))
})

##################################
# Clear GO Terms table selection #
##################################

# When clicking on gene table all rows are unselected
proxy_wat_go_terms_table = dataTableProxy('wat_go_terms_table')
observeEvent(input$wat_gene_table_rows_selected, {
  selectRows(proxy_wat_go_terms_table, NULL)
})

# All rows are unselected when KEGG row is selected
observeEvent(input$wat_kegg_table_rows_selected, {
  selectRows(proxy_wat_go_terms_table, NULL)
})

# When selecting points on plot all rows are unselected
observeEvent(input$wat_vol_plot_brush, {
  selectRows(proxy_wat_go_terms_table, NULL)
})

##############################
# Clear KEGG table selection #
##############################

# When clicking on gene table all rows are unselected
proxy_wat_kegg_table = dataTableProxy('wat_kegg_table')
observeEvent(input$wat_gene_table_rows_selected, {
  selectRows(proxy_wat_go_terms_table, NULL)
})

# When clicking on GO Terms table all rows are unselected
observeEvent(input$wat_go_terms_table_rows_selected, {
  selectRows(proxy_wat_kegg_table, NULL)
})

# When selecting points on plot all rows are unselected
observeEvent(input$wat_vol_plot_brush, {
  selectRows(proxy_wat_kegg_table, NULL)
})

##########################################
# Filter transcript data based on inputs #
##########################################
wat_highlighted_genes <- reactive({
  
  if(length(input$wat_gene_table_rows_selected) > 0){
    selected_genes <- joined_wat_results[input$wat_gene_table_rows_selected,]$ext_gene
    filter(joined_wat_results, ext_gene %in% selected_genes)
    
  } else if(length(input$wat_go_terms_table_rows_selected) > 0){
    selected_go_terms <- wat_go_terms[input$wat_go_terms_table_rows_selected,]$ID
    wat_go_terms %>% filter(ID %in% selected_go_terms) %>% 
      mutate(ENTREZID = strsplit(geneID,'/')) %>%
      dplyr::select(ONTOLOGY, ID, Description, ENTREZID) %>%
      unnest(cols = ENTREZID) %>% left_join(wat_universe, by = 'ENTREZID') %>%
      dplyr::rename(ext_gene = SYMBOL) %>% left_join(joined_wat_results, by = 'ext_gene') %>%
      filter(lrt_qval < 0.05)
  } else if(length(input$wat_kegg_table_rows_selected) > 0){
    
    selected_kegg <- wat_kegg[input$wat_kegg_table_rows_selected,]$Description
    
    left_join(joined_wat_results, kegg_genes, 'ext_gene') %>%
      filter(Description %in% selected_kegg, lrt_qval < 0.05)
  }
})

#######################
# Volcano Plot Output #
#######################

output$wat_vol_plot <- renderPlot({
  if(length(input$wat_gene_table_rows_selected) > 0){
    
    vol_plot_base(joined_wat_results, 'WAT') +
      
      geom_point(alpha = 0.8, 
                 size = 4, 
                 stroke = 0, 
                 data = wat_highlighted_genes(),
                 aes(color = ext_gene)) +
      geom_label_repel(data = wat_highlighted_genes(),
                       aes(label = external_transcript_name, 
                           color = ext_gene), 
                       size = 6, 
                       force = 20) +
      labs(color = "Gene Name")
  } else if(length(input$wat_go_terms_table_rows_selected) > 0){
    vol_plot_base(joined_wat_results, 'WAT') +
      geom_point(alpha = 1, 
                 size = 3, 
                 stroke = 0, 
                 data = wat_highlighted_genes(),
                 aes(color = Description)) +
      labs(color = "Gene Name")
  } else if(length(input$wat_kegg_table_rows_selected) > 0){
    vol_plot_base(joined_wat_results, 'WAT') +
      geom_point(alpha = 1, 
                 size = 3, 
                 stroke = 0, 
                 data = wat_highlighted_genes(),
                 aes(color = Description)) +
      labs(color = "Gene Name")
  } else{
    vol_plot_base(joined_wat_results, 'WAT')
  }
})

########################
# Gene Bar Plot Output #
########################
output$wat_gene_bar_plot <- renderPlot({
  if(length(input$wat_gene_table_rows_selected) > 0){
    
    gene_bp_gn(joined_wat_gene_results, wat_highlighted_genes(), 'WAT')
    
  } else if(length(input$wat_go_terms_table_rows_selected) > 0){
    
    gene_bp_go(wat_go_terms, 
               input$wat_go_terms_table_rows_selected,
               wat_universe,
               joined_wat_gene_results,
               'WAT',
               'go')
    
  } else if(length(input$wat_kegg_table_rows_selected) > 0){
    
    gene_bp_go(wat_kegg, 
               input$wat_kegg_table_rows_selected,
               kegg_genes,
               joined_wat_gene_results,
               'WAT',
               'kegg')
  } else{
    NULL
  }
})

##############################
# Transcript Bar Plot Output #
##############################
output$wat_transcript_bar_plot <- renderPlot({
  if(length(input$wat_gene_table_rows_selected) > 0){
    
    tran_bp_gn(wat_highlighted_genes(),'WAT')
    
  } else if(length(input$wat_go_terms_table_rows_selected) > 0){
    
    tran_bp_go(wat_highlighted_genes(),'WAT')
    
  } else if(length(input$wat_kegg_table_rows_selected) > 0){
    
    tran_bp_go(wat_highlighted_genes(),'WAT')
    
  } else{
    NULL
  }
})





#####################
# Gene Table Output #
#####################
output$wat_gene_table <- DT::renderDataTable({
  DT::datatable(joined_wat_results,
                selection = 'multiple',
                options = list(buttons = c('csv', 'excel'))
  )
})

########################
# GO Term Table Output #
########################
output$wat_go_terms_table <- DT::renderDataTable({
  DT::datatable(dplyr::select(wat_go_terms,
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
output$wat_kegg_table <- DT::renderDataTable({
  
    DT::datatable(wat_kegg,
                  selection = 'multiple',
                  options = list(buttons = c('csv', 'excel'))
    )
})
