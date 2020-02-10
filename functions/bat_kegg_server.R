###################
# KEGG MAP Output #
###################
output$bat_kegg_map <- renderImage({
  
  pathway <- filter(kegg_genes, Description == input$bat_kegg_map)$kegg_pathway[1]
  
  data <- kegg_genes %>% filter(kegg_pathway %in% pathway) %>%
    dplyr::rename(target_id = ext_gene) %>%
    left_join(joined_bat_gene_results, by = 'target_id') %>%
    filter(lrt_qval < 0.05)
  
  genes <- as.vector(log2(data$fc))
  names(genes) <- data$entrez_id
  pv.out <- pathview(genes, 
                     pathway.id = pathway, 
                     species = 'mmu', 
                     same.layer = F,
                     low = list(gene = 'blue'),
                     mid = list(gene = 'transparent'),
                     high = list(gene = 'red'))
  
  
  # Delete unused files
  currentFiles <- list.files('.')
  png <- grep('.png',currentFiles)
  if(length(png) > 1){
    file.remove(currentFiles[png][2])
  }
  
  xml <- grep('.xml',currentFiles)
  if(length(xml) > 0){
    file.remove(currentFiles[xml])
  }
  
  # Load in KEGG map PNG file
  imageFile <- list.files('.')[grep('pathview.png',list.files('.'))]
  list(src = imageFile)
  
  
  
}, 
deleteFile = TRUE
)