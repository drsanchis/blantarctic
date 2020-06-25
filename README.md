![](https://github.com/drsanchis/blantarctic/blob/master/images/header.png)


This package is capable of performing a tandem bioinformatic analysis (including BLAST, multiple alignment, phylogeny and conserved domains characterization) comparing one or more query proteins with the protein sequences contained in several GeneBank files.

Furthermore, this particular package:

- Stores intermediate data in the form of **Python objects** (instances), which potentially reduces the number of disk accesses during the execution.
- Generates a **comprehensive summary plot** that combines the alignment view with the conserved domains.
- Gives access to an **interactive version of those plots**, in which domain names are linked to their representation, and a pop-up window can be opened with further information about a certain domain in its context, as well as the protein in which it was found.

 

## Requirements

- **Python Standard Library modules:** os, sys, getopt, datetime, tkinter, textwrap, subprocess, shutil,  re, pickle. *(Are normally included in the standard Python installation.)*

- **BioPython modules:** Seq, SeqIO, Blast, ExPASy, Align, Phylo.

- **Other Python modules:** matplotlib.

- **Other software requirements:** muscle, NCBI BLAST

- **Databases:** prosite.dat and prosite.doc **must be manually included** in the package folder before executing the program (not uploaded to this repository due to size limitations).


## Usage
No instalation is required, just download and run.

```bash
python BLAntractic.py -q ../queries.fa -s ../subjects_folder/ -n project_name
```

Where:

- **queries.fa** is a single file in (multi)FASTA format, containing the sequence(s) to be used as query/queries.
- **subjects_folder** is a directory containing one or several GeneBank files (.gbff) of the genomes against which the BLAST analysis is going to be conducted.
- **project_name** is an arbitrary string, that will be included in all directories generated during the execution of the program.

Further options are:

- **–h** to show help.

- **--eval** followed by a <u>float</u>, to restrict the minimum e-value of BLAST results. Default: 10e-5.

- **--ident** followed by a <u>float</u>, to set a minimum identity threshold for the BLAST results. Default: 30%.

- **--cov** followed by a <u>float</u>, to set a minimum coverage threshold for the BLAST results. Default: 50%.

- **--exclude** in order to exclude from the search those PROSITE domains marked with `/SKIP-FLAG`, which are commonly found post-translational modifications in the majority of sequences. Removing these domains from the analysis may improve visibility of other more relevant domains.

- **-o** followed by the path to a BLAnatrctic project file (\*.bapj) to open it and directly access the interactive plots.
  

## Input

- Multiple **query** proteins must be contained in a single MULTIFASTA file (several FASTA files are not supported).
- **GeneBank** file of each subject genome must be provided and grouped in a <u>single folder</u>, whose path will be given as an input.



## Analysis

An overview of the workflow underlying the execution of BLAntarctic would include the following steps:

1. Processing the **input files** and generating suitable files for a BLAST analysis.
2. **BLAST** of the query proteins against a combined MULTIFASTA with all CDSs extracted from the GeneBank files.
3. BLAST results are stored in memory, as **instances** of the class 'Protein'. All information obtained about each protein in further analysis steps will be included in these instances. This allows for faster execution and plotting, since data are stored in memory and the <u>number and frequency of disk accesses is reduced</u>.
4. The **original sequence** of each BLAST match is retrieved from the input files. We opted for the original, complete sequences because their alignment and domains will arguably provide more meaningful insights into the differences between the original queries and the (hopefully) cold-adapted matches in the genomes from the Antarctic.
5. A multiple alignment is conducted on these sequences with the **MUSCLE** algorithm.
6. That multiple alignment is used to build a **phylogenetic tree**, using the MUSCLE algorithm as well.
7. **Conserved protein domains** contained in the PROSITE database are searched within the matching sequences retrieved form the BLAST analysis.
8. A **static plot** is generated, representing an alignment of the query and the complete matches. Additionally, relevant protein domains are marked on the query and subjects.
9. An **interactive version** of that same plot is presented. **Hovering** over the domains shows a label with its name and accession number. By **clicking** on the domains, a **pop-up window** can be opened, with relevant information about the domain itself, as well as about the protein in which it is located and the genome from which it was retrieved.



## Output

BLAntarctic generates two folders each time, always in the parent directory containing the package folder, i.e., output folders are created <u>outside</u> the package folder. 

The output files tree has the following structure:

```
parent_directory
├── BLAntractic_package
│   ├── prosite.dat
│   ├── prosite.doc
│   ├── BLAntarctic.py
│   └── ... (all other modules)
│
├── name_data_1
│   ├── raw_GBs
│   │   ├── Genome1.gbff
│   │   ├── Genome2.gbff
│   │   └── ...
│   ├── genomes_multifasta.fa
│   └── queries.fa
│
└── name_results_1
   ├── unaligned_matches
   │   ├── Q1_query1_matches.fa
   │   └── ...
   ├── aligned_matches
   │   ├── Q1_query1_aligned.fa
   │   └── ...
   ├── trees_nw
   │   ├── Q1_tree.nw
   │   └── ...
   ├── trees_plot
   │   ├── Q1_tree_plot.pdf
   │   └── ...
   ├── protein_domains
   │   ├── Q1_query1_domains.tsv
   │   └── ...
   ├── plots
   │   ├── Q1_query1.png
   │   └── ...
   ├── blast_results.tsv
   ├── filtered_blast_results.tsv
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



## Interactive results

Besides all the files enumerated in the previous section, an interactive interface opens at the end of the execution.

#### 1. Thumbnails menu

Upon finishing the analysis, an interactive Matplotlib plot is opened with thumbnails of the graphs generated. If any of the queries generated no results in the BLAST alignment, no thumbnail will be created for that specific query.

Select a plot by **clicking on its thumbnail**.

<img src="https://github.com/drsanchis/blantarctic/blob/master/images/thumbs_menu.png" width="600">

#### 2. Interactive alignment

An interactive plot of a query's alignment is show.

- The **query** protein is depicted in **blue** at the top of the plot.
- All **subject** proteins are aligned under the query, depicted in **grey**. Their corresponding locus tag and species is noted in the right margin. 
- The actual sequence of both query and subject proteins is depicted by a wider box, while spacers (gaps) resulting from the alignment are represented by a narrower line that horizontally connects those boxes.
- Colored regions on the proteins depict conserved **domains**. 
  - **Green boxes** for common <u>phosphorylation</u> sites.
  - **Yellow boxes** for <u>myristoylation</u> sites.
  - **Orange boxes** for all <u>other</u> motifs.

By **hovering** over the domains, a text box appears with the name of the domain and is PROSITE accession number. By **clicking ** on any of those domains, a pop-up window opens with extra information.

You can **return** to the thumbnails menu anytime by clicking the "Volver" button.

<img src="https://github.com/drsanchis/blantarctic/blob/master/images/interactive_plot.png" width="600">

#### 3. Expanded info pop-up

It may take a few seconds for the information to be retrieved and for this windows to be opened. Please, be patient. It consists of two tabs, where relevant information of both the **domain** itself (including its position in the alignment) and the **protein** in which is located is presented.

<img src="https://github.com/drsanchis/blantarctic/blob/master/images/popup.png" width="800">

#### 4. Re-opening a project file

If the user wants to access the interactive view above explained **at any other moment**, **the complete analysis doesn't need to be run**. By accessing the `project.bapj` file contained in the results directory, the interactive plots can be loaded in a much faster way. In order to open the project file, use:

```bash
python BLAntarctic.py -o [path_to_project_file]
```



