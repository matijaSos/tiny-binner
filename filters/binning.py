import data.resultdata as resdata
import filters.readprocessing as rstate
from ncbi.taxonomy.ranks import ranks as tax_ranks
from utils.location import Location

def bin_reads (
         reads, 
         cds_repository, 
         read2cds_repository, 
         tax_tree, 
         target_organism_taxids, 
         assign_multiple_alignment_to_same_organism,
         assign_multiple_alignment_to_different_organisms,
         output=True
         ):
    '''
    :param reads                    list of Read objects
    :param cds_repository           dictionary (key: UnityCds, value: CdsAlignment)
    :param read2cds_repository      dictionary (key: read_id, value: [CdsAlignment])
    :param tax_tree                 TaxTree object
    :param target_organism_taxids   list of identified organism tax IDs. These
                                    are the organisms all reads will be mapped to.

    :param assign_multiple_alignment_to_same_organism function
    :param assign_multiple_alignment_to_different_organisms function
    :param output (boolean) if True, binning status will be output to stdio
    '''
    organisms = _create_organisms(target_organism_taxids, tax_tree)

    for read in reads:
        '''
        print bin(read.status)
        for aln in read.get_alignments():
            print aln.tax_id,
        print 
        print 'zero aln read:', rstate.is_zero_alignment_read(read.status)
        print 'one aln read: ', rstate.is_single_alignment_read(read.status)
        print 'n aln read:   ', rstate.is_multiple_alignment_read(read.status)
        print 'target orgs:  ', rstate.is_mapped_to_target_organisms(read.status)
        print 'mixd orgs:    ', rstate.is_mapped_to_mixed_organisms(read.status)
        '''

    # STEP 0:
    #           | target | nontarget |
    #    coding |   0    |     X     |
    # noncoding |   0    |     X     |
        if rstate.is_zero_alignment_read(read.status)\
        or rstate.is_mapped_to_nontarget_organisms(read.status):
            # maybe write down additional organisms
            # skip this one
            # print '\tApplying step (1).'
            continue

    # STEP 1:
    #           | target | nontarget |
    #    coding |   1    |     0     |
    # noncoding |   0    |     0     |
        elif rstate.is_single_alignment_read(read.status)\
        and  rstate.is_mapped_to_single_coding_region(read.status)\
        and  rstate.is_mapped_to_coding_regions_of_single_target_organism(read.status):
            
            # print '\tApplying step (2).'
            # add gene to organism
            # add read to organism
            target_alignment = read.get_alignments(format=list)[0]
            target_cds = read2cds_repository[read.id]
            target_organism_taxid = determine_target_organism(
                                        target_alignment.tax_id, 
                                        target_organism_taxids, 
                                        tax_tree)
            add_cds_to_organism(organisms[target_organism_taxid], read, target_alignment) # Why not target_cds?

    # STEP 2:
    #           | target | nontarget |
    #    coding |   N    |     X     |
    # noncoding |   X    |     X     |
        elif rstate.is_mapped_to_coding_regions_of_multiple_target_organisms(read.status)\
        or   rstate.is_mapped_to_coding_regions_of_mixed_organisms(read.status):
            
            # print '\tApplying step (3).'
            # if organisms are related, assign to lowest - not implemented
            # if not, assign to best score
            # add read to organism
            # add gene to organism
            target_alignments = extract_coding_target_alignments(
                read.get_alignments(format=list),
                read2cds_repository,
                target_organism_taxids,
                tax_tree)
            best_alignment = find_best_alignment(target_alignments)
            target_organism_taxid = determine_target_organism(
                                        best_alignment.tax_id,
                                        target_organism_taxids,
                                        tax_tree)
            add_cds_to_organism(organisms[target_organism_taxid], read, best_alignment) # Why not CDS but whole alignment?

    # STEP 3:
    #           | target | nontarget |
    #    coding |   0    |     0     |
    # noncoding |   1+   |     X     |
        elif   rstate.is_not_mapped_to_coding_region(read.status)\
        and (rstate.is_mapped_to_target_organisms(read.status)\
            or rstate.is_mapped_to_mixed_organisms(read.status)):
            continue
            # assign to best score
            # add read to organism

            # print '\tApplying step (4).'
            cds_num = 0
            for aln in read.get_alignments(format=list):
                cds_num += len(aln.aligned_cdss)
            # print cds_num

            target_alignments = extract_noncoding_target_alignments(
                read.get_alignments(format=list),
                read2cds_repository,
                target_organism_taxids,
                tax_tree)
            best_alignment = find_best_alignment(target_alignments)
            target_organism_taxid = determine_target_organism(
                                        best_alignment.tax_id,
                                        target_organism_taxids,
                                        tax_tree)
            add_read_to_organism(organisms[target_organism_taxid], read, best_alignment)
        else:
            continue
            # print read.id, read.status
    return organisms

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

def determine_target_organism(alignment_tax_id, target_organism_taxids, tax_tree):
    if alignment_tax_id in target_organism_taxids:
        return alignment_tax_id
    parent_taxids = []
    for target_taxid in target_organism_taxids:
        if tax_tree.is_child(alignment_tax_id, target_taxid):
            parent_taxids.append(target_taxid)
    # if only one parent taxid present, return that one
    if len(parent_taxids) == 1:
        return parent_taxids[0]
    # if multiple parent taxids present, return the lowest rank organism (closest)
    else:
        taxids_by_rank = sorted(parent_taxids, 
                                key = lambda taxid: tax_ranks[tax_tree.nodes[taxid].rank])
        return taxids_by_rank[-1]

def add_cds_to_organism(organism, read, target_alignment):
    target_cdss = target_alignment.aligned_cdss
    assert(len(target_cdss) >= 1)
    binned_read = resdata.BinnedRead(read.id)
    if len(target_cdss) == 1:
        # do stuffs
        (target_cds, intersection) = target_cdss[0]
    else:
        # Get CDSs
        cdss = []
        for (cds, intersection) in target_cdss:
            cdss.append(cds)
        # Sort by length
        sorted_cdss = sorted(cdss, key = lambda cds: Location.from_location_str(cds.location).length())
        # Take the one with the biggest length
        target_cds = sorted_cdss[-1]

    if organism.contains_identified_coding_region(target_cds):
        identified_cds = organism.identified_coding_regions[target_cds]
        identified_cds.add_binned_read(binned_read)
    else:
        identified_cds = resdata.IdentifiedCds(target_cds)
        identified_cds.add_binned_read(binned_read)
        organism.add_identified_coding_region(identified_cds)

def find_best_alignment(target_alignments):
    # print 'Target alns:', len(target_alignments)
    sorted_alignments = sorted(target_alignments, key = lambda aln: aln.score)
    return sorted_alignments[-1]

def extract_coding_target_alignments(
                read_alignments,
                read2cds_repository,
                target_organism_taxids,
                tax_tree):
    coding_target_alignments = []
    for alignment in read_alignments:
        is_target = False
        if len(alignment.aligned_cdss) == 0:
            continue
        else:
            if alignment.tax_id in target_organism_taxids:
                is_target = True
            else:
                for tax_id in target_organism_taxids:
                    if tax_tree.is_child(alignment.tax_id, tax_id):
                        is_target = True
        if is_target:
            coding_target_alignments.append(alignment)
    return coding_target_alignments

def extract_noncoding_target_alignments(
                read_alignments,
                read2cds_repository,
                target_organism_taxids,
                tax_tree):
    noncoding_target_alignments = []
    for alignment in read_alignments:
        print alignment.tax_id
        is_target = False
        if alignment.tax_id in target_organism_taxids:
            is_target = True
        else:
            for tax_id in target_organism_taxids:
                if tax_tree.is_child(alignment.tax_id, tax_id):
                    is_target = True
        if is_target:
            print is_target
            noncoding_target_alignments.append(alignment)
    return noncoding_target_alignments

def add_read_to_organism(organism, read, target_alignment):
    binned_read = resdata.BinnedRead(read.id)
    organism.add_read_aligned_to_noncoding_region(binned_read)
