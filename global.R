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

########################
# Additional functions #
########################

formatHTML <- function(origFile, newFile){
  
  # This function reformats an HTLM file produced by Rmarkdown
  # so that it's compatible with a Shiny app that uses navbarMenu.
  
  require(magrittr)
  require(xml2)
  require(rvest)
  
  read_html(origFile) %>% 
    html_node('body') %>% 
    write_html(newFile)
  
}

