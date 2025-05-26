# protDB_shiny
This is a simple Shiny dashboard for generating and manipulating FASTA files for database searches for mass spectrometry analysis.

Available features include: 

* Download reference FASTA files from Uniprot and NCBI databases via FTP
* Searchable tables to find reference FASTA file at each location
* Upload of user-provided FASTA files to be processed
* Manipulations of FASTA files, including detection and removal of redundant entries, generation of reversed and/or shuffled sequences (for use as decoys for FDR control), addition of common contaminant proteins, and detection and optional removal of entries containing non-standard amino acids.

If no manipulation of the FASTA file is desired, this dashboard can be used to download the FASTA of choice by unchecking all options for manipulations/filtering.

Lastly, in addition to generating the output FASTA file, a text log file is also generated that captures the date and origin of the downloaded FASTA file, as well as basic summary statistics about the number of entries before and after any manipulations are performed.

The dashboard can be accessed at https://jmvarberg.shinyapps.io/protDB_generator/.

Databases last updated on May 26, 2025.


## Updates

* 10/22/24 - Updated UniProt database table and UI to allow option to download either the FASTA file containing only canonical sequences, or FASTA containing both canonical and isoform sequences. Prior to this update, only canonical sequences were downloaded when using UniProt option.
* 05/25/25 - Pushed update with additional FASTA Manipulation Tab that allows for interactive searching/filtering based on strings, as well as simple concatenation of two user-provided FASTA files.