library(shiny)
library(DT)
library(plotly)
library(dplyr)
library(Biostrings)
library(stringi)
library(stringr)
library(shinyalert)
library(shinybusy)


#Global:

ncbi_ids_reference <- readRDS("./databases/ncbi_reference_table.Rds")
uniprot_ids_reference <- readRDS("./databases/uniprot_ids_reference_table.Rds")

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
                       width = "100%"
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
    
    #When user provied FASTA is selected, generate the upload button
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
            
            #download the fasta.gz as a tempfile
            temp <- tempfile()
            download.file(fastaFile, temp)
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
        }
        
        #repeat for user upload
        if(input$dbSource == "user") {
        
            print(input$userFasta$datapath)
            print(input$userFasta$name)
            temp <-  input$userFasta$datapath
            fastaFile <- input$userFasta$name
            fastaName <- stringr::str_remove(fastaFile, ".fasta")
        }
        
        #begin downstream processing
        #read the FASTA file in using Biostrings readAAStringSet
        
        fasta_raw <- Biostrings::readAAStringSet(filepath = temp)
        print("FASTA completed downloading/uploading and is ready to process...")
        
        #coerce into a data frame with the names and AA sequences
        fasta_df <- data.frame(seq_name = names(fasta_raw),
                               sequence = paste(fasta_raw))
        
        num_total_prots <- nrow(fasta_df)
        print(paste("Total number of entries in the provided FASTA: ", num_total_prots, sep = ""))
        
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
            cont_fasta <- Biostrings::readAAStringSet("./NR-Contaminants_2020-12-11.fasta")
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
                    paste0("Input FASTA file downloaded/uploaded from ", fastaFile, " on ", date, "\n",
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
    })
}


shinyApp(ui, server)







