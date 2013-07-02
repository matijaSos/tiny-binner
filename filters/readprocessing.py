from utils import enum
#------------------------------------------------------------------------#
#------------ READ ANNOTATIONS REGRADING [ORGANISM/ALIGNMENT] -----------#
alignment_count_status = enum (                             #0b0000000xx
                                NO_ALIGNMENT                =0b0000000001,
                                ONE_ALIGNMENT               =0b0000000010,
                                MULTIPLE_ALIGNMENTS         =0b0000000011,
                                MASK                        =0b0000000011
                              )
organism_count_status = enum  (                             #0b000000xx00
                                NO_ORGANISMS                =0b0000000000,
                                ONE_ORGANISM                =0b0000000100,
                                MULTIPLE_ORGANISMS          =0b0000001000,
                                MASK                        =0b0000001100
                              )
organism_type_status  = enum  (                             #0b0000xx0000
                                TARGET_ORGANISM             =0b0000001000,
                                NON_TARGET_ORGANISM         =0b0000010000,
                                MIXED                       =0b0000011000,
                                MASK                        =0b0000011000
                              )
coding_region_aln_count_status = enum (                     #0b00xx000000
                                NO_CODING_ALIGNMENTS        =0b0000000000,
                                ONE_CODING_ALIGNMENT        =0b0001000000,
                                MULTIPLE_CODING_ALIGNMENTS  =0b0010000000,
                                MASK                        =0b0011000000
                              )
coding_region_alignment_status = enum (                     #0bxx00000000
                                ONE_TARGET_ORGANISM         =0b0000000000,
                                MULTIPLE_TARGET_ORGANISMS   =0b0100000000,
                                NONTARGET_ORGANISMS         =0b1000000000,
                                MIXED_ORGANISMS             =0b1100000000,
                                MASK                        =0b1100000000
                              )

def annotate_reads(all_nonhost_reads,read2cds_repository, tax_tree, target_organisms):
    '''
    :param all_nonhost_reads list of reads that have been assessed not 
    to belong to a potential host
    :param tax_tree TaxTree object
    :param target_organisms list of tax IDs of organisms that have been
    determined to be potential patogens
    '''

    for read in all_nonhost_reads:
        alignments = read.get_alignments(format=list)
        alignment_num = len(alignments)
        # Case 1: read has no alignments
        if alignment_num == 0:
            mark_zero_alignment_read(read)
        elif alignment_num == 1:
            mark_single_alignment_read(read, read2cds_repository, tax_tree, target_organisms)
        else:
            mark_multiple_alignment_read(read, read2cds_repository, tax_tree, target_organisms)

def mark_zero_alignment_read(read):
    status = alignment_count_status.NO_ALIGNMENT
    read.set_status(status)

def mark_single_alignment_read(read, read2cds_repository, tax_tree, target_organisms):
    alignment = read.get_alignments(format=list)[0]
    # alignment count status
    aln_count = alignment_count_status.ONE_ALIGNMENT
    # organism count status
    org_count = organism_count_status.ONE_ORGANISM
    # organism type status
    mapped_to_target = False
    for tax_id in target_organisms:
        if tax_tree.is_child(alignment.tax_id, tax_id):
            mapped_to_target = True
            break
    if mapped_to_target:
        org_type = organism_type_status.TARGET_ORGANISM
    else:
        org_type = organism_type_status.NON_TARGET_ORGANISM
    # coding region alignment count status
    if not read2cds_repository.has_key(read.id):
        cds_aln_count = coding_region_aln_count_status.NO_CODING_ALIGNMENTS
    else:
        cds_aln_count = coding_region_aln_count_status.ONE_CODING_ALIGNMENT
    # coding region organism type
    if mapped_to_target:
        cds_aln_type = coding_region_alignment_status.ONE_TARGET_ORGANISM
    else:
        cds_aln_type = coding_region_alignment_status.NONTARGET_ORGANISMS

    status = aln_count | org_count | org_type | cds_aln_count | cds_aln_type
    read.set_status(status)

def mark_multiple_alignment_read(read, read2cds_repository, tax_tree, target_organisms):
    alignments = read.get_alignments(format=list)
    # alignment count status
    aln_count = alignment_count_status.MULTIPLE_ALIGNMENTS
    # organism count
    organisms = set()
    for alignment in alignments:
        organisms.add(alignment.tax_id)
    if len(organisms) == 1:
        org_count = organism_count_status.ONE_ORGANISM
    else:
        org_count = organism_count_status.MULTIPLE_ORGANISMS
    # organism type 
    children_count = get_child_count(organisms, target_organisms, tax_tree)
    if children_count == len(alignments):
        org_type = organism_type_status.TARGET_ORGANISM
    elif children_count == 0:
        org_type = organism_type_status.NON_TARGET_ORGANISM
    else:
        org_type = organism_type_status.MIXED
    # coding region alignment type
    cds_aln_type = 0
    if not read2cds_repository.has_key(read.id):
        cds_aln_count = coding_region_aln_count_status.NO_CODING_ALIGNMENTS
    else:
        cdss = read2cds_repository[read.id]
        if len(cdss) == 1:
            cds_aln_count = coding_region_aln_count_status.ONE_CODING_ALIGNMENT
        else:
            cds_aln_count = coding_region_aln_count_status.MULTIPLE_CODING_ALIGNMENTS
        # coding region organism type
        coding_alns_tax_ids = set()
        for cds in cdss:
            coding_alns_tax_ids.add(cds.cds.taxon)
        children_count = get_child_count(coding_alns_tax_ids, target_organisms, tax_tree)
        if children_count == len(coding_alns_tax_ids):
            if children_count == 1:
                cds_aln_type = coding_region_alignment_status.ONE_TARGET_ORGANISM
            else:
                cds_aln_type = coding_region_alignment_status.MULTIPLE_TARGET_ORGANISMS
        elif children_count == 0:
            cds_aln_type = coding_region_alignment_status.NONTARGET_ORGANISMS
        else:
            cds_aln_type = coding_region_alignment_status.MIXED_ORGANISMS

    status = aln_count | org_count | org_type | cds_aln_count | cds_aln_type
    read.set_status(status)



def get_child_count(tax_ids, target_organisms, tax_tree):
    children = 0
    for tax_id in tax_ids:
        child = False
        for target_taxid in target_organisms:
            if tax_tree.is_child(tax_id, target_taxid):
                child = True
                break
        if child:
            children += 1
    return children
