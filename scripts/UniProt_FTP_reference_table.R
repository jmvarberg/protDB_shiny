#read in and clean up the uniprot database table
uniprot_ids_reference <- data.table::fread("./UniProt_species_dbs.rtf")
uniprot_ids_reference <- uniprot_ids_reference |> 
  janitor::clean_names() |> 
  dplyr::mutate(species_name = stringr::str_remove(species_name, "\\\\"),
                superregnum = snakecase::to_title_case(superregnum),
                fasta_file = paste0(superregnum, "/", proteome_id, "/", proteome_id, "_", tax_id, ".fasta.gz")) |> 
  tidyr::separate(col = species_name, into = c("genus", "species"), sep = " ", remove=FALSE) |> 
  dplyr::mutate(short_species = paste(stringr::str_sub(genus, 1, 1), species, sep="_")) |> 
  dplyr::rename(num_entries_main_fasta = number_1,
                num_entries_additional_fasta = number_2,
                num_entries_gene2acc_map = number_3) |> 
  dplyr::mutate(entry_number = dplyr::row_number()) |> 
  dplyr::select(species_name, everything())
saveRDS(uniprot_ids_reference, "./uniprot_ids_reference_table.Rds")