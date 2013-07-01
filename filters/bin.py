import data.resultdata as resdata

def bin (
         read_repository, 
         cds_repository, 
         read2cds_repository, 
         tax_tree, 
         target_organism_taxids, 
         assign_multiple_alignment_to_same_organism,
         assign_multiple_alignment_to_different_organisms,
         output=True
         ):
    '''
    :param read_repository dictionary (key: read_id, value: Read object)
    :param cds_repository dictionary (key: UnityCds, value: CdsAlignment)
    :param read2cds_repository dictionary (key: read_id, value: [CdsAlignment])
    :param tax_tree TaxTree object
    :param target_organism_taxids list of identified organism tax IDs. These
    are the organisms all reads will be mapped to.
    :param assign_multiple_alignment_to_same_organism function
    :param assign_multiple_alignment_to_different_organisms function
    :param output (boolean) if True, binning status will be output to stdio
    '''
    organisms = _create_organisms(target_organism_taxids, tax_tree)
    additional_organisms = []
    processed_reads = []

    #--------------------------------------------------#
    #--- stage1: annotate zero alignment reads --------#
    if output:
        print 'Stage 1: annotating zero alignment reads...'
    zero_alignment_reads = annotate_zero_alignment_reads(read_repository)
    processed_reads.extend(zero_alignment_reads)
    if output:
        print 'done. [ZERO_ALN_READS: %d]' % len(zero_alignment_reads)
    #--------------------------------------------------#
    #--- stage2: assign reads with only one alignment--#
    if output:
        print 'Stage 2: assigning reads with one alignment...'
    one_alignment_reads = assign_single_alignment_reads(read_repository, read2cds_repository, organisms, additional_organisms)
    processed_reads.extend(one_alignment_reads)
    if output:
        print 'done. [ONE_ALIGNMENT_READS: %d]' % num_assigned

    #--------------------------------------------------#
    #--- stage3: assign multiple alignment reads, -----#   
    #--- but to the same organism                 -----#

def _create_organisms(target_organism_taxids, tax_tree):
    '''
    Create Organism objects from list of target tax IDs.
    
    :param target_organism_taxids list of tax IDs
    :param tax_tree TaxTree object
    :rtype dictionary (key: tax_id, value Organism)
    '''
    organisms = {}
    for tax_id in target_organism_taxids:
        tax_node = tax_tree.nodes[tax_id]
        organisms[tax_id] = resdata.Organism(tax_id, tax_node.organism_name, tax_node.rank)
    return organisms    

def annotate_zero_alignment_reads (read_repository):
    '''
    Annotates zero alignment reads. Set both their binning and 
    mapping status to NO_ALIGNMENT.

    :param read_repository dictionary (key: read_id, value: Read)
    :rtype dictionary (key: read_id, value: BinnedRead)
    '''
    zero_aln_reads = {}
    for read in read_repository.values():
        if len(read.get_alignments(format=list)) == 0:
            binned_read = resdata.BinnedRead(read.id, 
                                             mapping_status=read2cds_alignment_status.NO_ALIGNMENT,
                                             binning_status=read2organism_binning_status.NO_ALIGNMENT
                                             )
            zero_aln_reads[read.id] = binned_read
    return zero_aln_reads
