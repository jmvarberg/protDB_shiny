library(shiny)
library(gridlayout)
library(DT)
library(plotly)
library(dplyr)
library(Biostrings)
library(stringi)
library(stringr)
library(shinyalert)
library(shinybusy)


#Global:

#load the uniprot ID table
uniprot_ids_reference <- readRDS("./uniprot_ids_reference_table.Rds")
ncbi_ids_reference <- readRDS("./ncbi_reference_table.Rds")
#get today's data in YYYY-MM-DD
date <- as.character(Sys.Date())

#function to read in FASTA file and convert to a data frame
process_fasta <- function(fasta) {
    
    fasta_raw <- Biostrings::readAAStringSet(fasta)
    df <- data.frame(seq_name = names(fasta_raw),
                     sequence = paste(fasta_raw))
}

#list of standard amino acids
standard_aas <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

#function to check whether a AA sequence contains 
test_for_ns_aas <- function(x) {
    
    check <- sapply(strsplit(x, "")[[1]], FUN = match, table = standard_aas)
    num_nas <- sum(is.na(check))
    if(num_nas != 0) {
        #print("Sequence contains non-standard amino acids")
        return(FALSE)
    }
    
    else( 
        return(TRUE)
    )
    
}

outpath <- getwd()

#NA default values for number reversed, number shuffled, contaminants added
num_reversed = "NA"
num_shuffled = "NA"
num_contaminants = "NA"
num_NS_removed = "NA"

ui <- navbarPage(
  title = "Protein Database Dashboard",
  selected = "Database Generation",
  header = add_busy_spinner(spin = "fading-circle"),
  collapsible = TRUE,
  theme = bslib::bs_theme(),
  tabPanel(
    title = "Database Generation",
    grid_container(
      layout = c(
        "userOptions tablePreview",
        "existing    Output      "
      ),
      row_sizes = c(
        "1.35fr",
        "0.65fr"
      ),
      col_sizes = c(
        "295px",
        "1fr"
      ),
      gap_size = "10px",
      grid_card(
        area = "userOptions",
        title = "Processing Options",
        radioButtons(
          inputId = "sourceSelection",
          label = "FASTA Source Options",
          choices = list(
            `Retrieve from external site` = "external",
            `User provided file` = "user"
          )
        ),
        fileInput("fasta",
          "Choose FASTA File",
          multiple = FALSE,
          accept = c(".fasta")
        ),
        radioButtons(
          inputId = "dbSource",
          label = "Select FASTA Source",
          choices = list(
            NCBI = "ncbi",
            UniProt = "uniprot"
          )
        ),
        actionButton(
          inputId = "submitGenerate",
          label = "Submit",
          width = "100%"
        ),
        checkboxInput(
          inputId = "nonredundant",
          label = "Make non-redundant",
          value = TRUE,
          width = "100%"
        ),
        checkboxInput(
          inputId = "contaminants",
          label = "Add common contaminants",
          value = TRUE,
          width = "100%"
        ),
        checkboxInput(
          inputId = "reverseSeqs",
          label = "Add reversed decoy sequences",
          value = FALSE,
          width = "100%"
        ),
        checkboxInput(
          inputId = "shuffleSeqs",
          label = "Add shuffled decoy sequences",
          value = TRUE,
          width = "100%"
        ),
        checkboxInput(
          inputId = "nonStandard",
          label = "Remove sequences containing non-standard amino acids?",
          value = FALSE,
          width = "100%"
        )
      ),
      grid_card(
        area = "tablePreview",
        title = "Database Reference Table",
        DT::dataTableOutput(
          outputId = "myTable",
          width = "100%"
        )
      ),
      grid_card(
        area = "existing",
        title = "Processing Log",
        textOutput(outputId = "textLog")
      ),
      grid_card(
        area = "Output",
        title = "Output Database",
        DTOutput(
          outputId = "outDB",
          width = "100%"
        ),
        downloadButton(
          "downloadData",
          "Download Database"
        )
      )
    )
  ),
  tabPanel(
    title = "Database Concatenation",
    grid_container(
      layout = c(
        "starting      concatInputTable ",
        "adding        concatInput2Table",
        "concatMessage concatDBout      "
      ),
      row_sizes = c(
        "1.66fr",
        "2.08fr",
        "1.58fr"
      ),
      col_sizes = c(
        "1fr",
        "1fr"
      ),
      gap_size = "10px",
      grid_card(
        area = "starting",
        title = "Select input FASTA file to modify:",
        fileInput("fastaInput",
          "Choose FASTA File",
          multiple = FALSE,
          accept = c(
            ".txt",
            ".fasta"
          )
        )
      ),
      grid_card(
        area = "adding",
        title = "Sequences to add: (choose file or paste input)",
        fileInput("fastaConcat",
          "Choose FASTA File",
          multiple = FALSE,
          accept = c(
            ".txt",
            ".fasta"
          )
        ),
        textInput(
          inputId = "myTextInput",
          label = "Copy/Paste Input:",
          value = "",
          width = "100%",
          placeholder = "Paste FASTA Sequences to add here..."
        )
      ),
      grid_card(
        area = "concatMessage",
        title = "Processing Log:",
        textOutput(outputId = "concatMessage")
      ),
      grid_card(
        area = "concatInputTable",
        DTOutput(
          outputId = "concatFile1",
          width = "100%"
        )
      ),
      grid_card(
        area = "concatInput2Table",
        DTOutput(
          outputId = "concatFile2",
          width = "100%"
        )
      ),
      grid_card(
        area = "concatDBout",
        DTOutput(
          outputId = "concatDB",
          width = "100%"
        )
      )
    )
  ),
  tabPanel(
    title = "Database Subsetting",
    grid_container(
      layout = c(
        "starting      subsettingInput",
        "subsetIDs     area5          ",
        "subsetMessage subsetDB       "
      ),
      row_sizes = c(
        "1.56fr",
        "1.77fr",
        "1.67fr"
      ),
      col_sizes = c(
        "0.64fr",
        "1.14fr"
      ),
      gap_size = "10px",
      grid_card(
        area = "starting",
        title = "Select input FASTA file to modify",
        fileInput("subfastaInput",
          "Choose FASTA File",
          multiple = FALSE,
          accept = c(
            ".txt",
            ".fasta"
          )
        )
      ),
      grid_card(
        area = "subsetIDs",
        title = "Upload File with IDs (or paste input):",
        fileInput("subfastaKeep",
          "Choose File:",
          multiple = FALSE,
          accept = c(
            ".txt",
            ".csv",
            ".fasta"
          )
        ),
        textInput(
          inputId = "subsetUserInput",
          label = "Copy/Paste IDs for filtering:",
          value = "",
          width = "100%",
          placeholder = "Paste list of IDs to filter from FASTA"
        ),
        radioButtons(
          inputId = "subsetOptions",
          label = "Provided IDs should be:",
          choices = list(
            Kept = "keep",
            Removed = "remove"
          ),
          width = "100%"
        )
      ),
      grid_card(
        area = "subsettingInput",
        title = "Input FASTA:",
        DTOutput(
          outputId = "subsetInput",
          width = "100%",
          height = "100%"
        )
      ),
      grid_card(
        area = "subsetMessage",
        title = "Processing Log:",
        textOutput(outputId = "textOutputSubset")
      ),
      grid_card(
        area = "area5",
        DTOutput(
          outputId = "subsetKeep",
          width = "100%",
          height = "100%"
        )
      ),
      grid_card(
        area = "subsetDB",
        DTOutput(
          outputId = "subsetDBout",
          width = "100%"
        )
      )
    )
  )
)

server <- function(input, output, session) {
    options(shiny.maxRequestSize=200*1024^2) #sets maxiumum file upload size to 200 MB
    options(shiny.error = browser)
    
    observe({
        if(input$dbSource == "ncbi") {
            output$myTable <- DT::renderDataTable(ncbi_ids_reference, server=TRUE, selection = "single", options = list(
                scrollY="true",
                scrollX="400px",
                pageLength = 10,
                lengthMenu = c(10, 25, 50, 100),
                dom = 'Blfrtip'
            )
            )
        }
        
        if(input$dbSource == "uniprot") {
            output$myTable <- DT::renderDataTable(uniprot_ids_reference, server=TRUE, selection = "single", options = list(
                scrollY="true",
                scrollX="400px",
                pageLength = 10,
                lengthMenu = c(10, 25, 50, 100),
                dom = 'Blfrtip'
            )
            )
        }
        
        if(!is.null(input$subfastaInput)) {
            output$subsetInput <- DT::renderDataTable(process_fasta(fasta = input$subfastaInput$datapath), options = list(
                scrollY="true",
                scrollX="400px",
                pageLength = 10,
                lengthMenu = c(10, 25, 50, 100),
                dom = 'Blfrtip'
            )
            )
        }
        
    })
    
    observeEvent(input$submitGenerate, {
        
        entry_id = input$myTable_rows_selected
        
        #if user selected external site and Uniprot, then download the FASTA file from UniProt FTP site
        if(input$sourceSelection == "external" & input$dbSource == "uniprot") {
            
            #get the link to the fasta file from the reference table
            fastaFile <- uniprot_ids_reference |> 
                dplyr::filter(entry_number == entry_id) |> 
                dplyr::pull(fasta_file)
            fastaFile <- paste0("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/", fastaFile)
            
            #get the genus_species name
            genus_spec <- uniprot_ids_reference |> 
                dplyr::filter(entry_number == entry_id) |> 
                dplyr::pull(short_species)
            
            zipped <- basename(fastaFile)
            unzipped <- stringr::str_remove(basename(fastaFile), ".gz")
            fastaName <- paste0(genus_spec, "_", "UniProt_", stringr::str_remove(basename(fastaFile), ".fasta.gz"))
            
            #if neither the gzipped or unzipped FASTAs exist in the downloads subdirectory, download and unzip
            if(!file.exists(paste0("./downloads/", zipped)) & !file.exists(paste0("./downloads/", unzipped))) {
                
                curl::curl_download(fastaFile, destfile = paste0("./downloads/", zipped))
                R.utils::gunzip(paste0("./downloads/", zipped))
            }
            
            #if the gzipped exists in the downloads directory, unzip
            if(file.exists(paste0("./downloads/", zipped)) & !file.exists(paste0("./downloads/", unzipped))) {
                
                R.utils::gunzip(paste0("./downloads/", zipped))
                
            }
            
            unzipped <- paste0("./downloads/", unzipped)
        }
        
        #repeat for NCBI option
        if(input$sourceSelection == "external" & input$dbSource == "ncbi") {
            
            #get the link to the fasta file at NCBI from the reference table
            fastaFile <- ncbi_ids_reference |> 
                dplyr::filter(entry_number == entry_id) |> 
                dplyr::pull(ftp_path)
            
            #add extension to get the protein FASTA file at the ftp link
            fastaFile <- paste0(fastaFile, "/", basename(fastaFile), "_protein.faa.gz")
            
            #get the genus_species name
            genus_spec <- ncbi_ids_reference |> 
                dplyr::filter(entry_number == entry_id) |> 
                dplyr::pull(short_species)
            
            
            zipped <- basename(fastaFile)
            unzipped <- stringr::str_remove(basename(fastaFile), ".gz")
            fastaName <- paste0(genus_spec, "_", "NCBI_", stringr::str_remove(basename(fastaFile), ".faa.gz"))
            
            #if neither the gzipped or unzipped FASTAs exist in the directory, download and unzip
            if(!file.exists(paste0("./downloads/", zipped)) & !file.exists(paste0("./downloads/", unzipped))) {
                
                curl::curl_download(fastaFile, destfile = paste0("./downloads/", zipped))
                R.utils::gunzip(paste0("./downloads/", zipped))
            }
            
            #if the gzipped exists in the downloads subdirectory, unzip
            if(file.exists(paste0("./downloads/", zipped)) & !file.exists(paste0("./downloads/", unzipped))) {
                
                R.utils::gunzip(paste0("./downloads/", zipped))
                
            }
            
            unzipped <- paste0("./downloads/", unzipped)
        }
        
        #handle user provided FASTA as input
        if(input$sourceSelection == "user") {
            
            print(input$fasta$datapath)
            print(input$fasta$name)
            unzipped = input$fasta$datapath
            fastaFile <- input$fasta$name
            fastaName <- stringr::str_remove(fastaFile, ".fasta")
            
        }
        
        #begin downstream processing
        #read the FASTA file in using Biostrings readAAStringSet
        
        fasta_raw <- Biostrings::readAAStringSet(filepath = unzipped)
        
        #coerce into a data frame with the names and AA sequences
        fasta_df <- data.frame(seq_name = names(fasta_raw),
                               sequence = paste(fasta_raw))
        
        num_total_prots <- nrow(fasta_df)
        
        num_NS <- num_total_prots - sum(sapply(fasta_df$sequence, test_for_ns_aas))
        
        #add flag onto seq_name for entries with non-standard amino acids
        fasta_df <- fasta_df |> 
            dplyr::rowwise() |> 
            dplyr::mutate(seq_name = if_else(test_for_ns_aas(sequence), seq_name, paste0(seq_name, "_NONSTANDARD")))
        
        #handle non-redundant
        if(input$nonredundant) {
            
            fasta_df <- fasta_df |> 
                dplyr::distinct(sequence, .keep_all = TRUE)
            
            num_nr_prots <- nrow(fasta_df)
            fastaName <- paste0(fastaName, "_NR")
            
        }
        
        #handle reversed decoy sequences
        if(input$reverseSeqs) {
            
            decoy_df <- fasta_df |> 
                dplyr::mutate(seq_name = paste0("Reverse_", seq_name),
                              sequence = stringi::stri_reverse(sequence))
            
            #add to forward data frame
            fasta_df <- fasta_df |> dplyr::bind_rows(decoy_df)
            fastaName <- paste0(fastaName, "_Reversed")
            
            #get number of reversed sequences
            num_reversed <- nrow(fasta_df |> dplyr::filter(str_detect(seq_name, "Reverse")))
        }
        
        #handle contaminants
        if(input$contaminants) {
            
            #read in the contaminants FASTA file 
            cont_fasta <- Biostrings::readAAStringSet("./NR-Contaminants_2020-12-11.fasta")
            cont_df <- data.frame(seq_name = names(cont_fasta),
                                  sequence = paste(cont_fasta))
            num_contaminants <- nrow(cont_df)
            
            #add to the fasta data frame
            fasta_df <- fasta_df |> dplyr::bind_rows(cont_df)
            fastaName <- paste0(fastaName, "_wCont")
        }
        
        
        #handle shuffling sequences
        if(input$shuffleSeqs) {
            
            shuffled_df <- fasta_df |> 
                dplyr::mutate(seq_name = paste0("SHUFFLED_", stringr::word(seq_name, 1), " ", "FALSE POSITIVE"),
                              sequence = stringi::stri_rand_shuffle(sequence))
            
            #combine the fasta file with the shuffled data frame
            fasta_df <- fasta_df |> 
                dplyr::bind_rows(shuffled_df)
            fastaName <- paste0(fastaName, "_SHUFF")
            
            #number of shuffled proteins
            num_shuffled <- nrow(fasta_df |> dplyr::filter(str_detect(seq_name, "SHUFFLED")))
        }
        
        #Handling to remove sequences with non-standard amino acids
        if(input$nonStandard) {
            
            num_starting_prots <- nrow(fasta_df)
            
            fasta_df <- fasta_df |> 
                dplyr::filter(!stringr::str_detect(seq_name, "NONSTANDARD"))
            
            fastaName <- paste0(fastaName, "_nsAA_REMOVED")
            
            num_NS_removed <- num_starting_prots - nrow(fasta_df)
            
        }
        
        #write the final FASTA file out
        #convert the object back from df to XStringSet
        
        output$outDB <- DT::renderDataTable(fasta_dff, server=TRUE, selection = "single", options = list(
            scrollY="true",
            scrollX="400px",
            pageLength = 10,
            lengthMenu = c(10, 25, 50, 100),
            dom = 'Blfrtip'
        ))
        
        fasta_dff <- Biostrings::AAStringSet(fasta_df$sequence)
        names(fasta_dff) <- fasta_df$seq_name
        
        # Downloadable csv of selected dataset ----
        output$downloadData <- downloadHandler(
            filename = function() {
                paste0(fastaName, "_", date, ".fasta")
            },
            content = function(file) {
                write.table(fasta_dff, file, row.names = FALSE)
            }
        )
        
        
        output$textOutput <- renderText({paste0("Input FASTA file downloaded from ", fastaFile, " on ", date, ". The processed FASTA file was saved at the following path: ", paste0(outpath, fastaName, "_", date, ".fasta"), "\n",
                                                "Total number of proteins in original FASTA: ", num_total_prots, ".\n",
                                                "Total number of non-redundant proteins: ", num_nr_prots, ".\n", 
                                                "Total number of contaminants added: ", num_contaminants, ".\n", 
                                                "Total number of reversed proteins: ", num_reversed, ". \n",
                                                "Total number of shuffled proteins: ", num_shuffled, ". \n", 
                                                "Total number of sequences containing non-standard amino acids in the input FASTA file: ", num_NS, ". \n", 
                                                "Total number of sequences containing non-standard amino acids in the input FASTA file that were removed: ", num_NS_removed, ". \n", 
                                                "Total number of proteins in final output FASTA: ", nrow(fasta_df), ".")})
        
        sink(paste0("./", fastaName, "_", date, "_log.txt"))
        cat(
            paste0("Input FASTA file downloaded from ", fastaFile, " on ", date, ". The processed FASTA file was saved at the following path: ", paste0("./", fastaName, "_", date, ".fasta"), "\n",
                   "Total number of proteins in original FASTA: ", num_total_prots, ".\n",
                   "Total number of non-redundant proteins: ", num_nr_prots, ".\n", 
                   "Total number of contaminants added: ", num_contaminants, ".\n", 
                   "Total number of reversed proteins: ", num_reversed, ". \n",
                   "Total number of shuffled proteins: ", num_shuffled, ". \n", 
                   "Total number of sequences containing non-standard amino acids in the input FASTA file: ", num_NS, ". \n", 
                   "Total number of sequences containing non-standard amino acids in the input FASTA file that were removed: ", num_NS_removed, ". \n", 
                   "Total number of proteins in final output FASTA: ", nrow(fasta_df), ".")
        )
        sink()
    }
    )
    
    observe({
        if(!is.null(input$fastaInput)) {
            print(input$fastaInput)
            output$concatFile1 <- DT::renderDataTable(process_fasta(fasta = input$fastaInput$datapath), options = list(scrollX=TRUE, scrollY = "250px"))
        }
        
        if(!is.null(input$fastaConcat)) {
            print(input$fastaConcat)
            output$concatFile2 <- DT::renderDataTable(process_fasta(fasta = input$fastaConcat$datapath), options = list(scrollX=TRUE, scrollY = "250px"))
        }
        if(!input$myTextInput == "") {
            print(input$myTextInput)
            sink("./temp.fasta")
            cat(input$myTextInput)
            sink()
            output$concatFile2 <- DT::renderDataTable(process_fasta(fasta = "./temp.fasta"), options = list(scrollX=TRUE))
        }
    }
    )
    
    observeEvent(input$submitConcat, {
        
        outname = input$concatName
        print(length(outname))
        
        if(length(outname) == 0) {
            # Show a modal when the button is pressed
            shinyalert("Oops!", "Please provide a output file name.", type = "error")
        }
        
        
        
        if(!is.null(input$fastaInput) & !is.null(input$fastaConcat)) {
            
            input_data <- Biostrings::readAAStringSet(input$fastaInput$datapath)
            number_orig_prots <- length(input_data)
            add_data <- Biostrings::readAAStringSet(input$fastaConcat$datapath) 
            number_added_prots <- length(add_data)
            combined_data <- c(input_data, add_data)
            num_final_prots <- length(combined_data)
            Biostrings::writeXStringSet(combined_data, filepath = paste0(outpath, input$concatName))
            
            output$concatMessage <- renderText({paste0("Concatinated FASTA file completed. Find at the following path: ", paste0(outpath, input$concatName, "\n"), ". Number of proteins in original FASTA: ", number_orig_prots, ". Number of proteins in second FASTA: ", number_added_prots, ". Number of proteins in final output FASTA: ", num_final_prots)
            })
            
        }
        
    }
    )
    
    observeEvent(input$submitSubset, {
        
        unzipped = input$subfastaInput$datapath
        fastaFile <- unzipped
        fastaName <- stringr::str_remove(basename(unzipped), ".fasta")
        
        fasta_raw <- Biostrings::readAAStringSet(unzipped)
        
        #coerce into a data frame with the names and AA sequences
        fasta_df <- data.frame(seq_name = names(fasta_raw),
                               sequence = paste(fasta_raw))
        
        num_total_prots <- nrow(fasta_df)
        
        #read in the IDs to filter out
        keep = input$subfastaKeep$datapath
        keep_ids <-data.table::fread(keep, header=F) |> dplyr::pull(V1)
        
        num_keep_ids <- length(keep_ids)
        
        #filter only keep rows that contain the keep IDs list
        filt_df <- fasta_df |> 
            dplyr::filter(grepl(paste(keep_ids, collapse="|"), seq_name))
        
        num_out_prots <- nrow(filt_df)
        
        out_df <- Biostrings::AAStringSet(filt_df$sequence)
        names(out_df) <- filt_df$seq_name
        Biostrings::writeXStringSet(out_df, filepath = paste0(outpath, input$subsetOutname))
        
        output$textOutputSubset <- renderText({paste0("Subset FASTA file completed. Find at the following path: ", paste0(outpath, input$subsetOutname, "\n"), ". Number of proteins in original FASTA: ", num_total_prots, ". Number of proteins in filtering set: ", num_keep_ids, ". Number of proteins in final output FASTA: ", num_out_prots)
            
        })
    }
    )
}


shinyApp(ui, server)
