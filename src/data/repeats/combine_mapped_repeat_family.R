library(tidyverse)
library(glue)

families <- system('ls data/external/repeats/"$family"*fa | xargs -n1 basename | sed "s/\\..*$//"', intern = TRUE)

# Exclude minor satellites since they are mostly not present
# and LINE_L2 as well
families <- families[!(families %in% c("minor_satellite", "LINE_L2"))]

meta <- read_csv('data/raw/main/tables/samples_barcodes.csv', col_types = cols())

output_file <- "data/processed/main/features/cell_to_repeat_families_umis.tsv"

repeats <- lapply(families, function(family) {

	files <- Sys.glob(glue('data/raw/main/scrnaseq/transcriptome/*/outs/aligned_{family}_mapped.sam'))

	cat(glue(">> Looking into family {family}...\n"))

	cell_to_repeats_umis <- tibble(file = files) %>% 
  	tidyr::extract(file, "barcode", "(SIGA\\w\\d{2})", remove = FALSE) %>%
 	mutate(data = map(file, ~ read_tsv(.x, col_names = F, col_types = 'cccccccccccccccccc') %>%
		select(X1) %>%
		mutate(cell = str_match(X1, "CB:Z:(\\w+-\\d+)")[,2],
  		       molecule = str_match(X1, "UB:Z:(\\w+)")[,2], 
  		       barcode = .x %>% str_extract("SIGA\\w\\d{2}")) %>%
  		left_join(meta %>% select(transfection_lane, barcode)) %>%
  		mutate(transfection_lane = paste0('-', str_replace(transfection_lane, '_', '-'))) %>%
  		mutate(cell = ifelse(is.na(cell), cell, str_replace(cell, '-1', transfection_lane))) %>%
  		select(-transfection_lane, -barcode) %>%
		filter(!is.na(cell)) %>%
		filter(!is.na(molecule)) %>%
		select(cell, molecule) %>% 
		group_by(cell) %>%
		summarise(n_umis = length(unique(molecule))))) %>%
	unnest(data) %>%
	select(-file) %>%
	mutate(repeat_family = family)

	cell_to_repeats_umis %>% write_tsv(glue("data/processed/main/features/cell_to_{family}_umis.tsv"))
	
	cell_to_repeats_umis
})

bind_rows(repeats) %>%
  write_tsv(output_file)

print(paste0("Saved in ", output_file))
