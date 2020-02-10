library(limma)
library(tidyverse)
library(BiocManager)
library(pathview)
library(org.Mm.eg.db)
library(xml2)

# Import list of kegg pathway names
paths.mmu <- kegga(joined_wat_results$ens_gene, species = 'Mmu')
paths.mmu <- tibble(pathway_id = rownames(paths.mmu),
                    pathway_name = paths.mmu$Pathway,
                    N = paths.mmu$N)
paths.mmu <- paths.mmu %>% separate(col = pathway_id, into = c('a','pathway_id'), sep = 'mcc') %>%
  dplyr::select(-a)


# Kegg
kegg_to_dataframe <- function(id,
                              name = NULL,
                              species = 'mmu'){
  
  # Download kegg xml file
  download.kegg(pathway.id = id, species = 'mmu',file.type = 'xml')
  
  # Read xml file in and delete from folder
  xml_file_name <- paste0(species, id, '.xml')
  x <- read_xml(xml_file_name)
  file.remove(xml_file_name)
  
  # Reformat xml file into a data frame
  x2 <- xml_children(x) %>% as.character() %>% unlist() %>% as_tibble()
  x2 <- x2[grep('gene',x2$value),]
  x3 <- x2 %>% separate(value, into = c('a','b'), sep = 'name=\"') %>%
    separate(b, into = c('b','c'), sep = '\" type') %>%
    mutate(list = as.list(strsplit(b, split = ' '))) %>% dplyr::select(list) %>%
    unnest(list) %>% separate(list, into = c('spec','entrez_id'), sep = ':') %>%
    filter(spec == 'mmu') %>% mutate(kegg_pathway = id, pathway_name = name) %>%
    dplyr::select(kegg_pathway, pathway_name, kegg_gene)
  return(x3)
}

# Remove existing kegg data frame if already exists
if(exists('kegg_genes')){
  rm(kegg_genes)
}

# Loop over paths.mmu data frame and list correspoding kegg genes
for(i in 1:nrow(paths.mmu)){
  if(!exists('kegg_genes')){
    kegg_genes <- kegg_to_dataframe(id = paths.mmu$pathway_id[i],
                                    name = paths.mmu$pathway_name[i])
  } else{
    kegg_genes <- rbind(kegg_genes, 
                        kegg_to_dataframe(id = paths.mmu$pathway_id[i],
                                          name = paths.mmu$pathway_name[i]))
  }
}

# Import data to convert entrez IDs to gene names
library(biomaRt)
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "mmusculus_gene_ensembl",
                         host = "dec2015.archive.ensembl.org")
# host = "ensembl.org")
t2g <- biomaRt::getBM(
  attributes = c(
    "external_gene_name",
    'entrezgene'
  ),
  mart = mart
) %>% as_tibble() %>% dplyr::rename(entrez_id = entrezgene,
                                    ext_gene = external_gene_name) %>%
  mutate(entrez_id = as.character(entrez_id))

# Save files
kegg_genes <- kegg_genes %>% dplyr::rename(Description = pathway_name)
kegg_genes <- kegg_genes %>% left_join(t2g, by = 'entrez_id')

saveRDS(paths.mmu,'data/paths.mmu.RDS')
saveRDS(kegg_genes,'data/kegg_genes.RDS')