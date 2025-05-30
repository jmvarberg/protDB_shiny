options(repos = BiocManager::repositories())
library(shiny)
library(DT)
library(dplyr)
library(Biostrings)
library(stringi)
library(stringr)
library(shinyalert)
library(shinybusy)
library(fst)


#Global:
ncbi_ids_reference <- fst::read.fst("./databases/ncbi_reference_table.fst")
uniprot_ids_reference <- fst::read.fst("./databases/uniprot_ids_reference_table.fst")

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

#function to check whether a AA sequence contains. This is really slow/inefficient right now, look at options for improving. 
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

#NA default values for number reversed, number shuffled, contaminants added
num_nr_prots = "NA"
num_reversed = "NA"
num_shuffled = "NA"
num_contaminants = "NA"
num_NS <- "NA"
num_NS_removed = "NA"
num_total_prots_canonical = "NA"
num_total_prots_iso = "NA"
num_total_combined = "NA"
num_total_prots = "NA"
fastaFile = "NA"
isoformFastaFile = "NA"

# Set up the UI for downloading and manipulating FASTA database files-----------
ui <- navbarPage(
    title = "Protein Database Dashboard",
    selected = "FASTA Generator",
    header = add_busy_spinner(spin = "fading-circle"),
    collapsible = TRUE,
    theme = bslib::bs_theme(),
    
    #Set up tab for downloading from external databases-------------------------
    tabPanel(
        title = "FASTA Generator",
        
        #set up row with two columns for input and database preview
        fluidRow(
            
            #first column - for user input
            column(3,
                   
                   #Buttons to select UniProt or NCBI reference database--------
                   radioButtons(
                       inputId = "dbSource",
                       label = "Select FASTA Source",
                       choices = list(
                           NCBI = "ncbi",
                           UniProt = "uniprot",
                           User = "user"
                       )
                   ),
                   
                   br(), 
                   
                   uiOutput("userUpload"), 
                   
                   br(), 
                   
                   #Add Checkboxes to select the options for FASTA manipulations
                   #conditional box if uniprot selected
                   conditionalPanel(
                       condition = "input.dbSource == 'uniprot'",
                       checkboxInput(
                           inputId = "canonical_only",
                           label = "Use Database with ONLY canonical sequences? If selected, isoforms will not be included.",
                           value = FALSE,
                           width = "100%"
                       )
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
                       inputId = "nonStandard_check",
                       label = "Check for sequences containing non-standard amino acids? (*this is slow)",
                       value = FALSE,
                       width = "100%"
                   ),
                   checkboxInput(
                       inputId = "nonStandard",
                       label = "Remove sequences containing non-standard amino acids?",
                       value = FALSE,
                       width = "100%"
                   ),
                   
                   br(),
                   
                   #Add action button to submit the request to download/process
                   actionButton(
                       inputId = "submitGenerate",
                       label = "Submit",
                       width = "100%"
                   ),
                   
                   br(),
                   
                   uiOutput("downloadData"),
                   
                   br(), 
                   
                   uiOutput("downloadLog")
            ),
            
            #Second column for displaying the reference tables
            column(9,
                   DT::dataTableOutput(
                       outputId = "myTable",
                       width = "95%"
                   )
            )
        )
    ),
    
    #Second tab for Fasta Manipulations (filtering, concatenation)
    tabPanel(
        title = "FASTA Manipulator",
        
        #set up row with two columns for input and database preview
        fluidRow(
            
            #first column - for user input
            column(3,
                   
                   #Buttons to select manipulation type--------
                   radioButtons(
                       inputId = "manipulationType",
                       label = "Select FASTA Manipulation Type",
                       choices = list(
                           Filter = "filter",
                           Concatenate = "concatenate"
                       )
                   ),
                   
                   #File Upload for FASTA file
                   fileInput(inputId = "userFasta_manipulate", 
                             label = "Choose FASTA to Manipulate", 
                             multiple = FALSE, 
                             accept = ".fasta"),
                   
                   #conditionally show second upload for concatenation
                   conditionalPanel(
                       condition = "input.manipulationType == 'concatenate'",
                       fileInput(inputId = "fastaToConcat", 
                                 label = "Choose Second FASTA to Concatenate with First", 
                                 multiple = FALSE, 
                                 accept = ".fasta")
                   ),
                   
                   #conditionally allow string input for filtering
                   #conditionally show second upload for concatenation
                   conditionalPanel(
                       condition = "input.manipulationType == 'filter'",
                       textInput(inputId = "filterString", 
                                 label = "Enter string/text in FASTA header to use for filtering.",
                                 placeholder = "Fragment, isoform, MOUSE_, etc."
                       )
                   ),
                   
                   br(),
                   
                   #Display text with number of entries in FASTA containint the provided string
                   conditionalPanel(
                       condition = "input.manipulationType == 'filter'",
                       textOutput("rowCountText")
                       ),
                   
                   #Add action button to submit the request to download/process
                   actionButton(
                       inputId = "submitManipulate",
                       label = "Submit for Manipulation",
                       width = "100%"
                   ),
                   br(),
                   
                   uiOutput("downloadManipulationDB"),
                   
                   br(), 
                   
                   uiOutput("downloadManipulationLog")
            ),
            
            #Second column for displaying the reference tables
            column(9,
                   DT::dataTableOutput(
                       outputId = "modFastaInput",
                       width = "95%"
                   )
            )
        )
    )
)


#Define server-side processing
server <- function(input, output, session) {
    options(shiny.maxRequestSize=200*1024^2) #sets maximum file upload size to 200 MB
    options(shiny.error = browser)
    
    #When NCBI is selected, then show the table of all of the available FASTA from NCBI as a searchable DT::dataTable
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
        
        #When Uniprot is selected, then show the table of all of the available FASTA from Uniprot as a searchable DT::dataTable
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
        
    })
    
    #When user provided FASTA is selected, generate the upload button
    observe({
        
        #handle user provided fasta file input----------------------------------
        if (input$dbSource == "user") {
            
            #make the download buttons appear if the job is submitted and completed
            output$userUpload <- renderUI({
                fileInput(inputId = "userFasta", label = "Choose FASTA", multiple = FALSE, accept = ".fasta")
            })
            
            
        }
        
    })
    
    #Reactive processing - once submit button is generated, the highlighted/selected line in the FASTA table
    #is used as user selection for downstream processing according to selected options
    
    observeEvent(input$submitGenerate, {
        
        #if user selected Uniprot, then download the FASTA file from UniProt FTP site
        if(input$dbSource == "uniprot") {
            
            #Get the row ID of the selected species to process
            entry_id = input$myTable_rows_selected
            
            #get the genus_species name
            genus_spec <- uniprot_ids_reference |> 
                dplyr::filter(entry_number == entry_id) |> 
                dplyr::pull(short_species)
            
            #logic to handle selection of canonical vs. with isoforms
            if(input$canonical_only) {
                
                #get the link to the fasta file from the reference table
                fastaFile <- uniprot_ids_reference |> 
                    dplyr::filter(entry_number == entry_id) |> 
                    dplyr::pull(fasta_file)
                fastaFile <- paste0("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/", fastaFile)
                
                fasta_zipped <- basename(fastaFile)
                unzipped <- stringr::str_remove(basename(fastaFile), ".gz")
                fastaName <- paste0(genus_spec, "_", "UniProt_", stringr::str_remove(basename(fastaFile), ".fasta.gz"))
                
                #download the fasta.gz as a tempfile
                temp <- tempfile()
                download.file(fastaFile, temp)
                
                fasta_raw <- Biostrings::readAAStringSet(filepath = temp)
                print("FASTA completed downloading/uploading and is ready to process...")
                
                #coerce into a data frame with the names and AA sequences
                fasta_df <- data.frame(seq_name = names(fasta_raw),
                                       sequence = paste(fasta_raw))
                
                num_total_prots_canonical <- nrow(fasta_df)
                print(paste("Total number of entries in the provided FASTA: ", num_total_prots_canonical, sep = ""))
                
            } else {
                
                #if user selected to also include the isoforms in "additional.fasta", first check whether there are any entries in the additional fasta. can use info in "num_entries_additional_fasta" to check
                additional_check <- uniprot_ids_reference |> 
                    dplyr::filter(entry_number == entry_id) |> 
                    dplyr::pull(num_entries_additional_fasta)
                
                print("Number of entries in additional fasta: ")
                print(additional_check)
                
                if(additional_check == 0) {
                    
                    shinyalert::shinyalert(title = "Additional FASTA with isoforms doesn't exist.", 
                                           text = paste0("Selected species does not have any isoforms provided in additional fasta file. Please check box to only use 'canonical' and resubmit."),
                                           type = "error")
                } 
                
                if(additional_check > 0) {
                    
                    #get the link to the isoform fasta file from the reference table
                    isoformFastaFile <- uniprot_ids_reference |> 
                        dplyr::filter(entry_number == entry_id) |> 
                        dplyr::pull(additional_fasta_file)
                    isoformFastaFile <- paste0("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/", isoformFastaFile)
                    
                    #get the link to the fasta file from the reference table
                    fastaFile <- uniprot_ids_reference |> 
                        dplyr::filter(entry_number == entry_id) |> 
                        dplyr::pull(fasta_file)
                    fastaFile <- paste0("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/", fastaFile)
                    
                    canonical_fasta_zipped <- basename(fastaFile)
                    canonical_unzipped <- stringr::str_remove(basename(fastaFile), ".gz")
                    canonical_fastaName <- paste0(genus_spec, "_", "UniProt_", stringr::str_remove(basename(fastaFile), ".fasta.gz"))
                    
                    #download the fasta.gz as a tempfile
                    temp_canonical <- tempfile()
                    download.file(fastaFile, temp_canonical)
                    
                    #repeat for isoform version
                    iso_fasta_zipped <- basename(isoformFastaFile)
                    iso_unzipped <- stringr::str_remove(basename(isoformFastaFile), ".gz")
                    iso_fastaName <- paste0(genus_spec, "_", "UniProt_", stringr::str_remove(basename(isoformFastaFile), "_additional.fasta.gz"))
                    
                    #download the fasta.gz as a tempfile
                    temp_iso <- tempfile()
                    download.file(isoformFastaFile, temp_iso)
                    
                    fastaName <- canonical_fastaName #give the name to be used to just match the canonical. This will be used to add NR, SHUFF etc downstream.
                
                    #read in and process fasta files
                    canonical_raw <- Biostrings::readAAStringSet(filepath = temp_canonical)
                    isoform_raw <- Biostrings::readAAStringSet(filepath = temp_iso)
                    print("FASTA files completed downloading/uploading and are ready to process...")
                    
                    #coerce into a data frame with the names and AA sequences
                    canonical_fasta_df <- data.frame(seq_name = names(canonical_raw),
                                                     sequence = paste(canonical_raw))
                    
                    num_total_prots_canonical <- nrow(canonical_fasta_df)
                    print(paste("Total number of entries in the canonical FASTA: ", num_total_prots_canonical, sep = ""))
                    
                    #coerce into a data frame with the names and AA sequences
                    iso_fasta_df <- data.frame(seq_name = names(isoform_raw),
                                               sequence = paste(isoform_raw))
                    
                    num_total_prots_iso <- nrow(iso_fasta_df)
                    print(paste("Total number of entries in the isoform FASTA: ", num_total_prots_iso, sep = ""))
                    
                    #concatenate into single fasta file for next processing steps
                    fasta_df <- canonical_fasta_df |> 
                        dplyr::bind_rows(iso_fasta_df)
                    
                    num_total_combined <- nrow(fasta_df)
                    print(paste("Total number of entries in the combined FASTA: ", num_total_combined, sep = ""))
                }
                
                
            }
            
            
        }
        
        #repeat for NCBI option
        if(input$dbSource == "ncbi") {
            
            #Get the row ID of the selected species to process
            entry_id = input$myTable_rows_selected
            
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
            
            #download the fasta.gz as a tempfile
            temp <- tempfile()
            download.file(fastaFile, temp)
            
            #read in and process fasta file
            fasta_raw <- Biostrings::readAAStringSet(filepath = temp)
            print("FASTA completed downloading/uploading and is ready to process...")
            
            #coerce into a data frame with the names and AA sequences
            fasta_df <- data.frame(seq_name = names(fasta_raw),
                                   sequence = paste(fasta_raw))
            
            num_total_prots <- nrow(fasta_df)
            num_total_prots_canonical <- num_total_prots
            print(paste("Total number of entries in the provided FASTA: ", num_total_prots, sep = ""))
            
        }
        
        #repeat for user upload
        if(input$dbSource == "user") {
            
            print(input$userFasta$datapath)
            print(input$userFasta$name)
            temp <-  input$userFasta$datapath
            fastaFile <- input$userFasta$name
            fastaName <- stringr::str_remove(fastaFile, ".fasta")
            
            #read in and process fasta
            #read in and process fasta file
            fasta_raw <- Biostrings::readAAStringSet(filepath = temp)
            print("FASTA completed downloading/uploading and is ready to process...")
            
            #coerce into a data frame with the names and AA sequences
            fasta_df <- data.frame(seq_name = names(fasta_raw),
                                   sequence = paste(fasta_raw))
            
            num_total_prots <- nrow(fasta_df)
            num_total_prots_canonical <- num_total_prots
            print(paste("Total number of entries in the provided FASTA: ", num_total_prots, sep = ""))
            
        }
        
        #begin downstream processing
        
        if(input$nonStandard_check) {
            
            print("Checking for entries with non-standard amino acids...")
            num_NS <- num_total_prots - sum(sapply(fasta_df$sequence, test_for_ns_aas))
            #add flag onto seq_name for entries with non-standard amino acids
            fasta_df <- fasta_df |> 
                dplyr::rowwise() |> 
                dplyr::mutate(seq_name = if_else(test_for_ns_aas(sequence), seq_name, paste0(seq_name, "_NONSTANDARD")))
        }
        
        #handle non-redundant
        if(input$nonredundant) {
            
            print("Finding and removing redundant entries...")
            
            fasta_df <- fasta_df |> 
                dplyr::distinct(sequence, .keep_all = TRUE)
            
            num_nr_prots <- nrow(fasta_df)
            fastaName <- paste0(fastaName, "_NR")
            
        }
        
        #handle contaminants
        if(input$contaminants) {
            
            print("Adding common contaminants")
            
            #read in the contaminants FASTA file 
            #cont_fasta <- Biostrings::readAAStringSet("./NR-Contaminants_2020-12-11.fasta")
            cont_fasta <- Biostrings::readAAStringSet("./NIHMS1873978-supplement-Contaminant_FASTA.fasta")
            cont_df <- data.frame(seq_name = names(cont_fasta),
                                  sequence = paste(cont_fasta))
            num_contaminants <- nrow(cont_df)
            
            #add to the fasta data frame
            fasta_df <- fasta_df |> dplyr::bind_rows(cont_df)
            fastaName <- paste0(fastaName, "_wCont")
        }
        
        #handle reversed decoy sequences
        if(input$reverseSeqs) {
            
            print("Generating reversed decoy sequences...")
            
            decoy_df <- fasta_df |> 
                dplyr::mutate(seq_name = paste0("Reverse_", seq_name),
                              sequence = stringi::stri_reverse(sequence))
            
            #add to forward data frame
            fasta_df <- fasta_df |> dplyr::bind_rows(decoy_df)
            fastaName <- paste0(fastaName, "_Reversed")
            
            #get number of reversed sequences
            num_reversed <- nrow(fasta_df |> dplyr::filter(str_detect(seq_name, "Reverse")))
        }
        
        #handle shuffling sequences
        if(input$shuffleSeqs) {
            
            print("Generating shuffled decoy sequences...")
            
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
            
            print("Removing entries with non-standard amino acids...")
            
            num_starting_prots <- nrow(fasta_df)
            
            fasta_df <- fasta_df |> 
                dplyr::filter(!stringr::str_detect(seq_name, "NONSTANDARD"))
            
            fastaName <- paste0(fastaName, "_nsAA_REMOVED")
            
            num_NS_removed <- num_starting_prots - nrow(fasta_df)
            
        }
        
        #write the final FASTA file out
        #convert the object back from df to XStringSet
        
        print("Generating output FASTA and processing log...")
        
        fasta_dff <- Biostrings::AAStringSet(fasta_df$sequence)
        names(fasta_dff) <- fasta_df$seq_name
        
        #make the download buttons appear if the job is submitted and completed
        output$downloadData <- renderUI({
            req(input$submitGenerate, exists("fasta_dff"))
            downloadButton(outputId = "downloadDataDB", label = "Download Database")
        })
        
        output$downloadLog<- renderUI({
            req(input$submitGenerate, exists("fasta_dff"))
            downloadButton(outputId = "downloadDataLog", label = "Download Log")
        })
        
        # Download FASTA file of processed database ----------------------------
        output$downloadDataDB <- downloadHandler(
            filename = function() {
                paste0(fastaName, "_", date, ".fasta")
            },
            content = function(file) {
                Biostrings::writeXStringSet(x = fasta_dff, filepath = file)
            }
        )
        
        #Download txt file of processing log -----------------------------------
        output$downloadDataLog<- downloadHandler(
            filename = function() {
                paste0(fastaName, "_", date, "_processing_log.txt")
            },
            content = function(file) {
                sink(file)
                cat(
                    paste0("Canonical FASTA file downloaded/uploaded from ", fastaFile, " on ", date, ".\n",
                           "Isoform FASTA file downloaded/uploaded from ", isoformFastaFile, " on ", date, ".\n",
                           "Total number of proteins in canonical FASTA: ", num_total_prots_canonical, ".\n",
                           "Total number of proteins in isoform FASTA: ", num_total_prots_iso, ".\n",
                           "Total number of proteins in combined canonical/isoform FASTA: ", num_total_combined, ".\n",
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
    })
    
    #reactive display of uploaded FASTA for manipulation panel.
    observe({
        
        req(input$userFasta_manipulate)
        
        #get path to uploaded file from temp name
        #print(input$userFasta_manipulate$datapath)
        #print(input$userFasta_manipulate$name)
        temp <-  input$userFasta_manipulate$datapath
        fastaFile <- input$userFasta_manipulate$name
        fastaName <- stringr::str_remove(fastaFile, ".fasta")
        
        #read in the FASTA and convert to data frame
        fasta <- Biostrings::readAAStringSet(filepath = temp)
        
        #coerce into a data frame with the names and AA sequences
        fasta_df <- data.frame(seq_name = names(fasta),
                               sequence = paste(fasta))
        
        output$modFastaInput <- DT::renderDataTable(fasta_df, server=TRUE, selection = "single", options = list(
            scrollY="true",
            scrollX="400px",
            pageLength = 10,
            lengthMenu = c(10, 25, 50, 100),
            dom = 'Blfrtip'
        )
        )
        
        #display how many rows contain the provided string
        # Reactive text output to display the count of rows containing the filter string
        output$rowCountText <- renderText({
            req(input$userFasta_manipulate, input$filterString) # Ensure input is available
            filtered_count <- sum(grepl(input$filterString, fasta_df$seq_name, ignore.case = TRUE))
            paste("Number of entries containing the string to filter:", filtered_count)
        })
    })
    
    
    #reactive processing - manipulation. Once submit button pushed for manipulation, perform the desired manipulation
    observeEvent(input$submitManipulate, {
        
        if(input$manipulationType == "filter") {
            
            #test functionality
            testFasta <- "./NIHMS1873978-supplement-Contaminant_FASTA.fasta"
            testFasta <- Biostrings::readAAStringSet(filepath = testFasta)
            #coerce into a data frame with the names and AA sequences
            fasta_df <- data.frame(seq_name = names(testFasta),
                                   sequence = paste(testFasta))

            rowsInput <- nrow(fasta_df)
            fasta_df_filtered <- fasta_df |> dplyr::filter(!stringr::str_detect(seq_name, "PIG"))

            rowsOutput <- nrow(fasta_df_filtered)
            rowsFiltered <- rowsInput - rowsOutput
            
            #Set input/filtered row counts to NULL
            rowsInput <- NULL
            rowsFiltered <- NULL
            rowsOutput <- NULL
            
            
            req(input$userFasta_manipulate, input$filterString)
            
            #read in the uploaded fasta file and make data frame
            temp <-  input$userFasta_manipulate$datapath
            fastaFile <- input$userFasta_manipulate$name
            fastaName <- stringr::str_remove(fastaFile, ".fasta")
            
            #read in the FASTA and convert to data frame
            fasta <- Biostrings::readAAStringSet(filepath = temp)
            
            #coerce into a data frame with the names and AA sequences
            fasta_df <- data.frame(seq_name = names(fasta),
                                   sequence = paste(fasta))
            
            rowsInput <- nrow(fasta_df)
            
            #filter out headers that contain the string provided
            fasta_df_filtered <- fasta_df %>% dplyr::filter(!stringr::str_detect(seq_name, input$filterString))
            
            rowsOutput <- nrow(fasta_df_filtered)
            rowsFiltered <- rowsInput - rowsOutput
            
            
            #generate output with download buttons
            fasta_df_filtered_out <- Biostrings::AAStringSet(fasta_df_filtered$sequence)
            names(fasta_df_filtered_out) <- fasta_df_filtered$seq_name
            
            print(fasta_df_filtered_out)
            print(rowsInput)
            print(rowsFiltered)
            print(rowsOutput)
            
            #make the download buttons appear if the job is submitted and completed
            output$downloadManipulationDB <- renderUI({
                req(input$submitManipulate, exists("fasta_df_filtered_out"))
                downloadButton(outputId = "downloadManipulatedFasta", label = "Download Manipulated FASTA")
            })
            
            output$downloadManipulationLog<- renderUI({
                req(input$submitManipulate, exists("fasta_df_filtered_out"))
                downloadButton(outputId = "downloadManipulatedLog", label = "Download FASTA Manipulation Log")
            })
            
            # Download FASTA file of processed database ----------------------------
            output$downloadManipulatedFasta <- downloadHandler(
                filename = function() {
                    paste0(fastaName, "_", date, "_filtered.fasta")
                },
                content = function(file) {
                    Biostrings::writeXStringSet(x = fasta_df_filtered_out, filepath = file)
                }
            )
            
            #Download txt file of processing log -----------------------------------
            output$downloadManipulatedLog<- downloadHandler(
                filename = function() {
                    paste0(fastaName, "_", date, "_filtered_processing_log.txt")
                },
                content = function(file) {
                    sink(file)
                    cat(
                        paste0( "Total number of proteins in input FASTA: ", rowsInput, ".\n",
                               "Total number of proteins in filtered output FASTA: ", rowsOutput, ".\n",
                               "Total number of proteins that contained the filtering string/text: ", rowsFiltered, ".\n",
                               "String/Text used for filtering: ", input$filterString)
                    )
                    sink()
                }
            )
            
        }
        
        if(input$manipulationType == "concatenate") {
            
            rowsOriginal <- NULL
            rowsSecond <- NULL
            rowsConcatenated <- NULL
            
            
            #read in first fasta and make data frame
            req(input$userFasta_manipulate, input$fastaToConcat)
            
            #read in the uploaded fasta file and make data frame
            temp <-  input$userFasta_manipulate$datapath
            fastaFile <- input$userFasta_manipulate$name
            fastaName <- stringr::str_remove(fastaFile, ".fasta")
            
            #read in the FASTA and convert to data frame
            fasta <- Biostrings::readAAStringSet(filepath = temp)
            
            #coerce into a data frame with the names and AA sequences
            fasta_df <- data.frame(seq_name = names(fasta),
                                   sequence = paste(fasta))
            
            rowsOriginal <- nrow(fasta_df)
            
            
            #read in second fasta and make data frame
            #read in the uploaded fasta file and make data frame
            temp_2 <-  input$fastaToConcat$datapath
            fastaFile_2 <- input$fastaToConcat$name
            fastaName_2 <- stringr::str_remove(fastaFile_2, ".fasta")
            
            #read in the FASTA and convert to data frame
            fasta_2 <- Biostrings::readAAStringSet(filepath = temp_2)
            
            #coerce into a data frame with the names and AA sequences
            fasta_df_2 <- data.frame(seq_name = names(fasta_2),
                                   sequence = paste(fasta_2))
            
            rowsSecond <- nrow(fasta_df_2)
            
            
            #combine with bind_rows
            fasta_concatenated <- dplyr::bind_rows(fasta_df, fasta_df_2)
            
            rowsConcatenated <- nrow(fasta_concatenated)
            
            
            #generate output database and log
            #generate output with download buttons
            fasta_concatenated_out <- Biostrings::AAStringSet(fasta_concatenated$sequence)
            names(fasta_concatenated_out) <- fasta_concatenated$seq_name
            
            #make the download buttons appear if the job is submitted and completed
            output$downloadManipulationDB <- renderUI({
                req(input$submitManipulate, exists("fasta_concatenated_out"))
                downloadButton(outputId = "downloadConcatDB", label = "Download Manipulated FASTA")
            })
            
            output$downloadManipulationLog<- renderUI({
                req(input$submitManipulate, exists("fasta_concatenated_out"))
                downloadButton(outputId = "downloadConcatLog", label = "Download FASTA Manipulation Log")
            })
            
            # Download FASTA file of processed database ----------------------------
            output$downloadConcatDB <- downloadHandler(
                filename = function() {
                    paste0(fastaName, "_", date, "_concatenated.fasta")
                },
                content = function(file) {
                    Biostrings::writeXStringSet(x = fasta_concatenated_out, filepath = file)
                }
            )
            
            #Download txt file of processing log -----------------------------------
            output$downloadConcatLog<- downloadHandler(
                filename = function() {
                    paste0(fastaName, "_", date, "_concatenated_processing_log.txt")
                },
                content = function(file) {
                    sink(file)
                    cat(
                        paste0( "Total number of proteins in first input FASTA: ", rowsOriginal, ".\n",
                                "Total number of proteins in second input FASTA: ", rowsSecond, ".\n",
                                "Total number of proteins in concatenated FASTA: ", rowsConcatenated, ".\n",
                                "First FASTA File Name: ", fastaName, ".\n", 
                                "Second FASTA File Name: ", fastaName_2, ".")
                    )
                    sink()
                }
            )
            
            
            
            
            
        }
        
        
        
    })
    
}


shinyApp(ui, server)







