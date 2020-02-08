
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

# When selecting points on plot all rows are unselected
observeEvent(input$bat_vol_plot_brush, {
  selectRows(proxy_bat_go_terms_table, NULL)
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
  }
})

#######################
# Volcano Plot Output #
#######################
bat_vol_plot <- 
  reactive(
    {
      ggplot(joined_bat_results, 
             aes(x = log2(fc), 
                 y = -log10(lrt_qval))) + 
        geom_point(alpha = 0.25, 
                   stroke = 0,
                   data = filter(joined_bat_results, 
                                 !is.na(Direction))) +
        geom_hline(yintercept = -log10(0.05), 
                   color = 'red', 
                   linetype = 'dashed') + 
        geom_vline(xintercept = 0, 
                   color = 'red', 
                   linetype = 'dashed') + 
        xlim(-9,9) +
        #scale_color_manual(values = c('#4393c3','#878787','#d6604d')) +
        xlab('log2(fold change)') + 
        ylab('-log10(qValue)') +
        ggtitle('BAT Cort Pellet vs Placebo Pellet - Gene level') +
        theme_bw() + 
        theme(text = element_text(size=20),
              legend.position="bottom")
    }
  )

output$bat_vol_plot <- renderPlot({
  if(length(input$bat_gene_table_rows_selected) > 0){
    bat_vol_plot() +
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
    bat_vol_plot() +
      geom_point(alpha = 1, 
                 size = 3, 
                 stroke = 0, 
                 data = bat_highlighted_genes(),
                 aes(color = Description)) +
      labs(color = "Gene Name")
  } else{
    bat_vol_plot()
  }
})

########################
# Gene Bar Plot Output #
########################
output$bat_gene_bar_plot <- renderPlot({
  if(length(input$bat_gene_table_rows_selected) > 0){
    joined_bat_gene_results %>% filter(target_id %in% bat_highlighted_genes()$ext_gene,
                                       lrt_qval < 0.05) %>%
      ggplot(., aes(x = fct_reorder(target_id, fc), y = log2(fc), fill = target_id)) +
      geom_bar(stat = 'identity', aes(alpha = I(((lrt_qval < 0.05) / 1.25) + 0.2))) +
      geom_hline(yintercept = 0) +
      theme_bw() +
      guides(alpha = FALSE) +
      theme(axis.text.x = element_text(angle = 50, hjust = 1),
            panel.spacing = unit(0, "lines"),
            legend.position="bottom") +
      ggtitle('Selected gene abundances - gene level') +
      xlab('') + ylab('log2 bat fold change cort\n pellet / placebo pellet')
    
  } else if(length(input$bat_go_terms_table_rows_selected) > 0){
    
    selected_go_terms <- bat_go_terms[input$bat_go_terms_table_rows_selected,]$ID
    bat_go_terms %>% filter(ID %in% selected_go_terms) %>% 
      mutate(ENTREZID = strsplit(geneID,'/')) %>%
      dplyr::select(ONTOLOGY, ID, Description, ENTREZID) %>%
      unnest(cols = ENTREZID) %>% left_join(bat_universe, by = 'ENTREZID') %>%
      dplyr::rename(target_id = SYMBOL) %>% left_join(joined_bat_gene_results, by = 'target_id') %>%
      filter(lrt_qval < 0.05) %>% 
      
      ggplot(., aes(x = target_id, y = log2(fc), fill = Description)) +
      geom_bar(stat = 'identity') +
      geom_hline(yintercept = 0) +
      facet_wrap(~Description, scales = 'free_x') +
      theme_bw() +
      guides(alpha = FALSE) +
      theme(axis.text.x = element_text(angle = 50, hjust = 1),
            panel.spacing = unit(0, "lines"),
            legend.position="bottom") +
      ggtitle('Selected gene abundances - gene level') +
      xlab('') + ylab('log2 bat fold change cort\n pellet / placebo pellet')
    
  } else{
    NULL
  }
})

##############################
# Transcript Bar Plot Output #
##############################
output$bat_transcript_bar_plot <- renderPlot({
  if(length(input$bat_gene_table_rows_selected) > 0){
    ggplot(bat_highlighted_genes(), aes(x = external_transcript_name, y = log2(fc), fill = ext_gene)) +
      geom_bar(stat = 'identity', aes(alpha = I(((lrt_qval < 0.05) / 1.25) + 0.2))) +
      geom_hline(yintercept = 0) +
      theme_bw() +
      guides(alpha = FALSE) +
      theme(panel.spacing = unit(0, "lines"),
            axis.text.x = element_text(angle = 50, hjust = 1),
            legend.position="bottom") +
      ggtitle('Selected gene abundances - transcript level') +
      xlab('') + ylab('log2 bat fold change cort\n pellet / placebo pellet')
    
  } else if(length(input$bat_go_terms_table_rows_selected) > 0){
    ggplot(bat_highlighted_genes(), aes(x = external_transcript_name, y = log2(fc), fill = Description)) +
      geom_bar(stat = 'identity') +
      geom_hline(yintercept = 0) +
      facet_wrap(~Description, scales = 'free_x') +
      theme_bw() +
      guides(alpha = FALSE) +
      theme(panel.spacing = unit(0, "lines"),
            axis.text.x = element_text(angle = 50, hjust = 1),
            legend.position="bottom") +
      ggtitle('Selected gene abundances - transcript level') +
      xlab('') + ylab('log2 bat fold change cort\n pellet / placebo pellet')
    
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