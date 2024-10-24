#read in and clean up the uniprot database table
temp_uniprot <- tempfile()
curl::curl_download(url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README", destfile = temp_uniprot)

readme <- readLines(temp_uniprot)
first_entry <- min(grep(pattern = "^UP", readme, ignore.case = FALSE, value = FALSE))
last_entry <- max(grep(pattern = "^UP", readme, ignore.case = FALSE, value = FALSE))
num_rows <- last_entry - (first_entry-1)
uniprot_ids_reference <- data.table::fread(temp_uniprot, skip = (first_entry - 2), nrows = num_rows, header = T)
uniprot_ids_reference <- uniprot_ids_reference |> 
    janitor::clean_names() |> 
    dplyr::mutate(species_name = stringr::str_remove(species_name, "\\\\"),
                  superregnum = snakecase::to_title_case(superregnum),
                  fasta_file = paste0(superregnum, "/", proteome_id, "/", proteome_id, "_", tax_id, ".fasta.gz"),
                  additional_fasta_file = paste0(superregnum, "/", proteome_id, "/", proteome_id, "_", tax_id, "_additional.fasta.gz")) |> 
    tidyr::separate(col = species_name, into = c("genus", "species"), sep = " ", remove=FALSE) |> 
    dplyr::mutate(short_species = paste(stringr::str_sub(genus, 1, 1), species, sep="_")) |> 
    dplyr::rename(num_entries_main_fasta = number_1,
                  num_entries_additional_fasta = number_2,
                  num_entries_gene2acc_map = number_3) |> 
    dplyr::mutate(entry_number = dplyr::row_number()) |> 
    dplyr::select(species_name, everything())
#saveRDS(uniprot_ids_reference, "./databases/uniprot_ids_reference_table.Rds")
fst::write.fst(uniprot_ids_reference, "./databases/uniprot_ids_reference_table.fst")