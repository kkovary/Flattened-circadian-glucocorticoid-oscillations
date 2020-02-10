library(shiny)
library(tidyverse)
library(shinyWidgets)
library(shinycssloaders)
library(ggrepel)
library(DT)
library(shinyjs)

###################
# Load data files #
###################

# WAT Data
joined_wat_results <- readRDS('data/joined_wat_results.RDS') %>% 
  filter(!is.na(lrt_qval), !is.na(ext_gene)) %>%
  dplyr::select(ext_gene, external_transcript_name, everything())

joined_wat_gene_results <- readRDS('data/joined_wat_gene_results.RDS') %>% filter(!is.na(lrt_qval))

wat_go_terms <- readRDS('data/wat_go_terms.RDS')
wat_universe <- readRDS('data/universe_wat.RDS')

# BAT Data
joined_bat_results <- readRDS('data/joined_bat_results.RDS') %>% 
  filter(!is.na(lrt_qval), !is.na(ext_gene)) %>%
  dplyr::select(ext_gene, external_transcript_name, everything())

joined_bat_gene_results <- readRDS('data/joined_bat_gene_results.RDS') %>% filter(!is.na(lrt_qval))

bat_go_terms <- readRDS('data/bat_go_terms.RDS')
bat_universe <- readRDS('data/universe_bat.RDS')

# KEGG data
kegg_genes <- readRDS('data/kegg_genes.RDS')

bat_kegg <- left_join(joined_bat_results, kegg_genes, 'ext_gene') %>%
  group_by(kegg_pathway, Description) %>% 
  summarise(num_transcripts = n(),
            zscore = (sum(log2(fc) > 0 & lrt_qval < 0.05, na.rm = T) - sum(log2(fc) < 0 & lrt_qval < 0.05, na.rm = T)) / sqrt(n()))

wat_kegg <- left_join(joined_wat_results, kegg_genes, 'ext_gene') %>%
  group_by(kegg_pathway, Description) %>% 
  summarise(num_transcripts = n(),
            zscore = (sum(log2(fc) > 0 & lrt_qval < 0.05, na.rm = T) - sum(log2(fc) < 0 & lrt_qval < 0.05, na.rm = T)) / sqrt(n()))


########################
# Additional functions #
########################


#####################
# Volcano Plot Base #
#####################

vol_plot_base <- function(base_data, depot){
  ggplot(base_data, 
         aes(x = log2(fc), 
             y = -log10(lrt_qval))) + 
    geom_point(alpha = 0.25, 
               stroke = 0,
               data = filter(base_data, 
                             !is.na(Direction))) +
    geom_hline(yintercept = -log10(0.05), 
               color = 'red', 
               linetype = 'dashed') + 
    geom_vline(xintercept = 0, 
               color = 'red', 
               linetype = 'dashed') + 
    xlim(-9,9) +
    xlab('log2(fold change)') + 
    ylab('-log10(qValue)') +
    ggtitle(paste0(depot,' Cort Pellet vs Placebo Pellet - Gene level')) +
    theme_bw() + 
    theme(text = element_text(size=20),
          legend.position="bottom")
}

##########################################
# Transcript bar plot: input = gene name #
##########################################
tran_bp_gn <- function(selected_genes, 
                       depot){
  
  ggplot(selected_genes, aes(x = external_transcript_name, y = log2(fc), fill = ext_gene)) +
    geom_bar(stat = 'identity', aes(alpha = I(((lrt_qval < 0.05) / 1.25) + 0.2))) +
    geom_hline(yintercept = 0)  +
    theme_bw() +
    guides(alpha = FALSE) +
    theme(panel.spacing = unit(0, "lines"),
          axis.text.x = element_text(angle = 50, hjust = 1),
          legend.position="bottom") +
    ggtitle('Selected gene abundances - transcript level') +
    xlab('') + ylab(paste0('log2 ',depot,' fold change cort\n pellet / placebo pellet'))
}


########################################
# Transcript bar plot: input = GO/KEGG #
########################################
tran_bp_go <- function(selected_genes, 
                       depot){
  
  ggplot(selected_genes, aes(x = external_transcript_name, y = log2(fc), fill = Description)) +
    geom_bar(stat = 'identity') +
    geom_hline(yintercept = 0) +
    facet_wrap(~Description, scales = 'free_x') +
    theme_bw() +
    theme(panel.spacing = unit(0, "lines"),
          axis.text.x = element_text(angle = 50, hjust = 1),
          legend.position="bottom") +
    guides(alpha = FALSE) +
    ggtitle('Selected gene abundances - transcript level') +
    xlab('') + ylab(paste0('log2 ',depot,' fold change cort\n pellet / placebo pellet'))
}

####################################
# Gene bar plot: input = gene name #
####################################
gene_bp_gn <- function(selected_genes, 
                       gene_names, 
                       depot){
  
  selected_genes %>% filter(target_id %in% gene_names$ext_gene) %>%
    ggplot(., aes(x = fct_reorder(target_id, fc), y = log2(fc), fill = target_id)) +
    geom_bar(stat = 'identity', aes(alpha = I(((lrt_qval < 0.05) / 1.25) + 0.2))) +
    geom_hline(yintercept = 0) +
    theme_bw() +
    guides(alpha = FALSE) +
    theme(axis.text.x = element_text(angle = 50, hjust = 1),
          panel.spacing = unit(0, "lines"),
          legend.position="bottom") +
    ggtitle('Selected gene abundances - gene level') +
    xlab('') + ylab(paste0('log2 ',depot,' fold change cort\n pellet / placebo pellet'))
}


##################################
# Gene bar plot: input = GO/KEGG #
##################################
gene_bp_go <- function(all_cat,
                       selected_cat,
                       universe,
                       gene_df,
                       depot,
                       type){
  if(type == 'go'){
    selected_cat <- all_cat[selected_cat,]$ID
    
    all_cat %>% filter(ID %in% selected_cat) %>% 
      mutate(ENTREZID = strsplit(geneID,'/')) %>%
      dplyr::select(ONTOLOGY, ID, Description, ENTREZID) %>%
      unnest(cols = ENTREZID) %>% left_join(universe, by = 'ENTREZID') %>%
      dplyr::rename(target_id = SYMBOL) %>% left_join(gene_df, by = 'target_id') %>%
      filter(lrt_qval < 0.05) %>% 
      
      ggplot(., aes(x = target_id, y = log2(fc), fill = Description)) +
      geom_bar(stat = 'identity') +
      geom_hline(yintercept = 0) +
      facet_wrap(~Description, scales = 'free_x') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 50, hjust = 1),
            panel.spacing = unit(0, "lines"),
            legend.position="bottom") +
      ggtitle('Selected gene abundances - gene level') +
      xlab('') + ylab(paste0('log2 ',depot,' fold change cort\n pellet / placebo pellet'))
  } else if(type == 'kegg'){
    
    selected_cat <- all_cat[selected_cat,]$Description
    
    universe %>% filter(Description %in% selected_cat) %>%
      dplyr::rename(target_id = ext_gene) %>%
      left_join(gene_df, by = 'target_id') %>%
      filter(lrt_qval < 0.05) %>% 
      
      ggplot(., aes(x = target_id, y = log2(fc), fill = Description)) +
      geom_bar(stat = 'identity') +
      geom_hline(yintercept = 0) +
      facet_wrap(~Description, scales = 'free_x') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 50, hjust = 1),
            panel.spacing = unit(0, "lines"),
            legend.position="bottom") +
      ggtitle('Selected gene abundances - gene level') +
      xlab('') + ylab(paste0('log2 ',depot,' fold change cort\n pellet / placebo pellet'))
  }
  
}