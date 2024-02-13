#download and format the NCBI refseq table

curl::curl_download(url="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt", destfile = "./NCBI_refseq_summary_table.txt")

#read in and clean up formatting
ncbi_refseq_table <- data.table::fread("./NCBI_refseq_summary_table.txt", skip=1)
ncbi_refseq_table <- ncbi_refseq_table |> 
    janitor::clean_names() |> 
    dplyr::filter(refseq_category %in% c("reference genome", "representative genome")) |> 
    dplyr::select(organism_name, infraspecific_name, taxid, species_taxid, seq_rel_date, ftp_path) |> 
    dplyr::mutate(entry_number = dplyr::row_number()) |> 
    tidyr::separate(col = organism_name, into = c("genus", "species"), sep = " ", remove=FALSE) |> 
    dplyr::mutate(short_species = paste(stringr::str_sub(genus, 1, 1), species, sep="_"))
saveRDS(ncbi_refseq_table, "./ncbi_reference_table.Rds")
