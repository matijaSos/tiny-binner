CDS_FASTA="cds.fa"
BLAST_DB="db/cds"
READS_INPUT_FILE="reads.fa"
BLAST_OUTPUT="blast.out"
BINNER_PATH="../../metaBinner"
ALN_FILE="../scripts/test_generator/lisa.in"

makeblastdb -in $CDS_FASTA -dbtype nucl -out $BLAST_DB
blastn -query $READS_INPUT_FILE -db $BLAST_DB -out $BLAST_OUTPUT -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore'
cd $BINNER_PATH
python snippets/parse_blast_output.py $BLAST_OUTPUT $ALN_FILE 'qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore'