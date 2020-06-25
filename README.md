# BLAntarctic v0.1



![image-20200625124510323](/Users/diegoruiz/Library/Application Support/typora-user-images/image-20200625124510323.png)





## Requirements

- **Python Standard Library modules:** os, sys, getopt, datetime, tkinter, textwrap, subprocess, shutil,  re, pickle.

- **BioPython modules:** Seq, SeqIO, Blast, ExPASy, Align, Phylo.

- **Other Python modules:** matplotlib.

- **Other software requirements:** muscle, NCBI BLAST

- **Databases:** prosite.dat and prosite.doc **must be manually included** in the package folder before executing the program (not uploaded to this repository due to size limitations).

  

## Usage

```bash
python BLAntractic.py -q ../queries.fa -s ../subjects_folder/ -n project_name
```

Where:

- **queries.fa** is a single file in (multi)FASTA format, containing the sequence(s) to be used as query/queries.
- **subjects_folder** is a directory containing one or several GeneBank files (.gbff) of the genomes against which the BLAST analysis is going to be conducted.
- **project_name** is an arbitrary string, that will be included in all directories generated during the execution of the program.

Further options are:

- **–h** to show a help.
- **––eval** followed by a <u>float</u>, to restrict the minimum e-value of BLAST results. Default: 10e-5.
- **––ident** followed by a <u>float</u>, to set a minimum identity threshold for the BLAST results. Default: 30%.
- **––cov** followed by a <u>float</u>, to set a minimum coverage threshold for the BLAST results. Default: 50%.
- **––exclude** in order to exclude from the search those PROSITE domains marked with `/SKIP-FLAG`, which are commonly found post-translational modifications in the majority of sequences. Removing these domains from the analysis may improve visibility of other more relevant domains.

## Input

- Multiple **query** proteins must be contained in a single MULTIFASTA file (several FASTA files are not supported).
- **GeneBank** file of each subject genome must be provided and grouped in a <u>single folder</u>, whose path will be given as an input.



## Analysis

An overview of the workflow underlying the execution of BLAntarctic would include the following steps:

1. Processing the **input files** and generating suitable files for a BLAST analysis.
2. **BLAST** of the query proteins against a combined MULTIFASTA with all CDSs extracted from the GeneBank files.
3. BLAST results are stored in memory, as **instances** of the class 'Protein'. All information obtained for each protein in further analysis steps will be included in these instances. This allows for faster execution and plotting, since data are stored in memory and the <u>number and frequency of disk accesses is reduced</u>.
4. 



## Output

BLAntarctic generates two folders each time, always in the parent directory containing the package folder, i.e., output folders are created <u>outside</u> the package folder. 

The output files tree has the following structure:

```
parent_directory
├── BLAntractic_package
│		├── prosite.dat
│		├── prosite.doc
│		├── BLAntarctic.py
│		└── ... (all other modules)
│
├── name_data_1
│		├── raw_GBs
│		│		├── Genome1.gbff
│		│		├── Genome2.gbff
│		│		└── ...
│		├── genomes_multifasta.fa
│		└── queries.fa
│
└── name_results_1
		├── unaligned_matches
		│		├── Q1_query1_matches.fa
		│		└── ...
		├── aligned_matches
		│		├── Q1_query1_aligned.fa
		│		└── ...
		├── trees_nw
		│		├── Q1_tree.nw
		│		└── ...
		├── trees_plot
		│		├── Q1_tree_plot.pdf
		│		└── ...
		├── protein_domains
		│		├── Q1_query1_domains.txt
	  │		└── ...
		├── plots
		│		├── Q1_query1.png
		│		└── ...
		├── blast_results.txt
		├── filtered_blast_results.txt
 		└── project.bapj

```



- The **name_data** directory contains the original input used for the analysis (queries and subjects). Both the original GeneBank files ("raw_GBs" folder) and a compound multifasta derived from them ("genomes_multifasta.fa") are included in this directory.

- The **name_results** directory contains all the outputs generated during the execution.

  - The **raw BLAST results** ("<u>blast_results.txt</u>"), as well as the resulting file after applying an identity percentage filter ("<u>filtered\_blast\_results.txt</u>"). Those filtered results are further rearranged in separate files (one for each query), that are stored in the "<u>unaligned_matches</u>" folder.
  - The **alignments** resulting from MUSCLE, contained in the "<u>aligned_matches</u>" folder.
  - The **phylogenetic trees** built from the alignments previously mentioned, both in Newick format ("<u>trees_nw</u>" folder) and their corresponding plots in *.pdf files ("<u>trees_plot</u>" folder).
  - Lists of the **conserved domains** found in the matching proteins for each query, in the "<u>protein_domains</u>" folder.
  - Combined **plots** for the results derived from each query, in the "<u>plots</u>" folder. These plots represent the BLAST results, multiple alignment and protein domains in a single diagram.
  - A **BLAntractic project file (*.bapj)** that can be opened anytime with the `-o` option, and gives access to the interactive results without having to run the whole analysis.

  

