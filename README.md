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

* **get_ribosomal_CDS_data** - Get binning output with statistics of each processing phase.

### Usage example

Download any dataset from http://www.hmpdacc.org/RSEQ/.<br/>
Align it to reference genome, for example using bwa mem algorithm. (http://bio-bwa.sourceforge.net/bwa.shtml)

If you have acquired .sam file as output, use <code>formats/sam2input.py</code> to get .in file.<br/>
If you have acquired .blast file as output, use <code>formats/blast2input.py</code> to get .in file.<br/>

Having .in file ready, run <code>python snippets/get_ribosomal_CDS_data input_file.in output_folder</code>.<br/>
Add <code>--remove_host</code> for host.
