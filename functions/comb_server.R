
#######################
# Input based filters #
#######################


##############################
# Clear gene table selection #
##############################

# When clicking on GO Terms table all rows are unselected
proxy_comb_gene_table = dataTableProxy('comb_gene_table')
observeEvent(input$comb_go_terms_table_rows_selected, {
  selectRows(proxy_comb_gene_table, NULL)
})

# All rows are unselected when KEGG row is selected
observeEvent(input$comb_kegg_table_rows_selected, {
  selectRows(proxy_comb_gene_table, NULL)
})

# When selecting points on plot all selected gene names are selected
observeEvent(input$comb_vol_plot_brush, {
  dims <- input$comb_vol_plot_brush
  genes <- filter(joined_comb_results,
                  -log10(lrt_qval) > dims$ymin, 
                  -log10(lrt_qval) < dims$ymax,
                  log2(fc) > dims$xmin,
                  log2(fc) < dims$xmax,
                  !duplicated(ext_gene))$ext_gene
  
  selectRows(proxy_comb_gene_table, which(joined_comb_results$ext_gene %in% genes))
})

##################################
# Clear GO Terms table selection #
##################################

# When clicking on gene table all rows are unselected
proxy_comb_go_terms_table = dataTableProxy('comb_go_terms_table')
observeEvent(input$comb_gene_table_rows_selected, {
  selectRows(proxy_comb_go_terms_table, NULL)
})

# All rows are unselected when KEGG row is selected
observeEvent(input$comb_kegg_table_rows_selected, {
  selectRows(proxy_comb_go_terms_table, NULL)
})

# When selecting points on plot all rows are unselected
observeEvent(input$comb_vol_plot_brush, {
  selectRows(proxy_comb_go_terms_table, NULL)
})

##############################
# Clear KEGG table selection #
##############################

# When clicking on gene table all rows are unselected
proxy_comb_kegg_table = dataTableProxy('comb_kegg_table')
observeEvent(input$comb_gene_table_rows_selected, {
  selectRows(proxy_comb_go_terms_table, NULL)
})

# When clicking on GO Terms table all rows are unselected
observeEvent(input$comb_go_terms_table_rows_selected, {
  selectRows(proxy_comb_kegg_table, NULL)
})

# When selecting points on plot all rows are unselected
observeEvent(input$comb_vol_plot_brush, {
  selectRows(proxy_comb_kegg_table, NULL)
})

##########################################
# Filter transcript data based on inputs #
##########################################
comb_highlighted_genes <- reactive({
  
  if(length(input$comb_gene_table_rows_selected) > 0){
    selected_genes <- joined_comb_results[input$comb_gene_table_rows_selected,]$ext_gene
    filter(joined_comb_results, ext_gene %in% selected_genes)
    
  } else if(length(input$comb_go_terms_table_rows_selected) > 0){
    selected_go_terms <- comb_go_terms[input$comb_go_terms_table_rows_selected,]$ID
    comb_go_terms %>% filter(ID %in% selected_go_terms) %>% 
      mutate(ENTREZID = strsplit(geneID,'/')) %>%
      dplyr::select(ONTOLOGY, ID, Description, ENTREZID) %>%
      unnest(cols = ENTREZID) %>% left_join(comb_universe, by = 'ENTREZID') %>%
      dplyr::rename(ext_gene = SYMBOL) %>% left_join(joined_comb_results, by = 'ext_gene') %>%
      filter(lrt_qval < 0.05)
  } else if(length(input$comb_kegg_table_rows_selected) > 0){
    
    selected_kegg <- comb_kegg[input$comb_kegg_table_rows_selected,]$Description
    
    left_join(joined_comb_results, kegg_genes, 'ext_gene') %>%
      filter(Description %in% selected_kegg, lrt_qval < 0.05)
  }
})

#######################
# Volcano Plot Output #
#######################

output$comb_vol_plot <- renderPlot({
  if(length(input$comb_gene_table_rows_selected) > 0){
    
    vol_plot_base(joined_comb_results, 'comb') +
      
      geom_point(alpha = 0.8, 
                 size = 4, 
                 stroke = 0, 
                 data = comb_highlighted_genes(),
                 aes(color = ext_gene)) +
      geom_label_repel(data = comb_highlighted_genes(),
                       aes(label = external_transcript_name, 
                           color = ext_gene), 
                       size = 6, 
                       force = 20) +
      labs(color = "Gene Name")
  } else if(length(input$comb_go_terms_table_rows_selected) > 0){
    vol_plot_base(joined_comb_results, 'comb') +
      geom_point(alpha = 1, 
                 size = 3, 
                 stroke = 0, 
                 data = comb_highlighted_genes(),
                 aes(color = Description)) +
      labs(color = "Gene Name")
  } else if(length(input$comb_kegg_table_rows_selected) > 0){
    vol_plot_base(joined_comb_results, 'comb') +
      geom_point(alpha = 1, 
                 size = 3, 
                 stroke = 0, 
                 data = comb_highlighted_genes(),
                 aes(color = Description)) +
      labs(color = "Gene Name")
  } else{
    vol_plot_base(joined_comb_results, 'comb')
  }
})

########################
# Gene Bar Plot Output #
########################
output$comb_gene_bar_plot <- renderPlot({
  if(length(input$comb_gene_table_rows_selected) > 0){
    
    gene_bp_gn(joined_comb_gene_results, comb_highlighted_genes(), 'comb')
    
  } else if(length(input$comb_go_terms_table_rows_selected) > 0){
    
    gene_bp_go(comb_go_terms, 
               input$comb_go_terms_table_rows_selected,
               comb_universe,
               joined_comb_gene_results,
               'comb',
               'go')
  } else if(length(input$comb_kegg_table_rows_selected) > 0){
    
    gene_bp_go(comb_kegg, 
               input$comb_kegg_table_rows_selected,
               kegg_genes,
               joined_comb_gene_results,
               'comb',
               'kegg')
  } else{
    NULL
  }
})

##############################
# Transcript Bar Plot Output #
##############################
output$comb_transcript_bar_plot <- renderPlot({
  if(length(input$comb_gene_table_rows_selected) > 0){
    
    tran_bp_gn(comb_highlighted_genes(),'comb')
    
  } else if(length(input$comb_go_terms_table_rows_selected) > 0){
    
    tran_bp_go(comb_highlighted_genes(),'comb')
    
  } else if(length(input$comb_kegg_table_rows_selected) > 0){
    
    tran_bp_go(comb_highlighted_genes(),'comb')
    
  } else{
    NULL
  }
})





#####################
# Gene Table Output #
#####################
output$comb_gene_table <- DT::renderDataTable({
  DT::datatable(joined_comb_results,
                selection = 'multiple',
                options = list(buttons = c('csv', 'excel'))
  )
})

########################
# GO Term Table Output #
########################
output$comb_go_terms_table <- DT::renderDataTable({
  DT::datatable(dplyr::select(comb_go_terms,
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
output$comb_kegg_table <- DT::renderDataTable({
  
  DT::datatable(comb_kegg,
                selection = 'multiple',
                options = list(buttons = c('csv', 'excel'))
  )
})


