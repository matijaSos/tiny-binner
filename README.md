# tiny-binner

Tiny binner is a metagenomic binning tool experimenting with a few different ideas.

## LCA - based
 
Metagenomic sample is taken as an input and present species are identified by the Lowest Common Ancestor algorithm.
For testing purposes MetaSim is used to generate syntethic metagenomic samples with chosen species and abundances.
Current functionalities are located in folder **snippets/**:

* **test_lca** - Evaluate LCA algorithm with given MetaSim FASTA file
* **test_megan** - Evaluate MEGAN with given MetaSim FASTA file
* **lca_vs_megan** - Compare LCA and MEGAN output with given MetaSim FASTA file

## Ribosomal coding sequences based

Idea here is to identify present species based on the presence of ribosomal coding sequences in the given metagenomic
sample.

**get_ribosomal_CDS_data** - Get binning output with statistics of each processing phase.
