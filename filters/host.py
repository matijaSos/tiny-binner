from operator import attrgetter

def filter_potential_hosts_alignments(reads, tax2category, potential_hosts, delete_host_alignments, filter_unassigned, unassigned_taxid=None):
    '''
    Method for marking (or deleting) potential host alignments. 
    Some organisms, however, do not fall either in potential hosts, nor in 
    microbes. Such organisms are some artificially created DNA strands or 
    unclassified environmental samples. For such organisms you can choose
    whether to filter them out or not.

    :param reads list or iterator of Read objects
    :param tax2category dict(key=taxid, value=category_taxid). Categories are
    listed in ncbi/taxonomy/organisms
    :param potential_hosts list of potential host organism taxids
    :param delete_host_alignments (boolean) if True, host alignments are 
    set to None. If False, their potential_host status is set to True/False.
    :param filter_unassigned (boolean) if True, unassigned organisms are
    filtered out.
    :param unassigned_taxid (int) which taxid means organism doesn't fall
    into neither category.
    '''

    #---- HOW TO FILTER READ (DELETE/MARK) -----#
    filter_alignment = determine_filtering_method (delete_host_alignments)
    #---- WHAT TO DO WITH UNASSIGNED TAXIDS ----# 
    if filter_unassigned:
        potential_hosts.append(unassigned_taxid)

    #-------------- FILTERING ------------------#
    for read in reads:
        alignments = read.get_alignments(format=list)
        potential_host_indexes = []

        for i in range(0, len(alignments)):
            read_alignment = alignments[i]

            if read_alignment.tax_id is None:
                potential_host_indexes.append(i)

            taxid_category = tax2category.get(read_alignment.tax_id, unassigned_taxid)
            if taxid_category in potential_hosts:
                potential_host_indexes.append(i)
            else:
                read_alignment.potential_host = False

        if delete_host_alignments:
            for i in reversed(potential_host_indexes):
                alignments.pop(i)
        else:
            for i in potential_host_indexes:
                alignments[i].potential_host = True

        read.set_alignments(alignments)

    return reads



def filter_potential_host_reads(reads, tax2category, potential_hosts, delete_reads, filter_unassigned, unassigned_taxid, find_host_status, percentage=0.5):
    '''
    Method for marking (or deleting) potential host reads. 
    Some organisms, however, do not fall neither in potential hosts, nor in 
    microbes. Such organisms are some artificially created DNA strands or 
    unclassified environmental samples. For such organisms you can choose
    whether to filter them out or not.

    :param reads list or iterator of Read objects
    :param tax2category dict(key=taxid, value=category_taxid). Categories are
    listed in ncbi/taxonomy/organisms
    :param potential_hosts list of potential host organism taxids
    :param delete_host_alignments (boolean) if True, host alignments are 
    set to None. If False, their potential_host status is set to True/False.
    :param filter_unassigned (boolean) if True, unassigned organisms are
    filtered out, meaning they are treated the same as potential host organisms.
    :param unassigned_taxid (int) which taxid means organism doesn't fall
    into neither category.
    :find_host_status (function) method for filtering host reads. User can choose between
    is_best_score_host, perc_of_host_alignments_larger_than and are_all_alignments_host
    '''

    #---- WHAT TO DO WITH UNASSIGNED TAXIDS ----# 
    if filter_unassigned:
        potential_hosts.append(unassigned_taxid)
    #---- CHECK FILTERING METHOD IS VALID  -----#
    if find_host_status not in (is_best_score_host, perc_of_host_alignments_larger_than, are_all_alignments_host):
        raise ValueError('''Supplied find_host_status method does not match any of the options for filtering: \
            ... is_best_score_host, perc_of_host_alignments_larger_than and are_all_alignments_host''')

    if isinstance(reads, iter):
        reads = list(reads)

    # 1. Identify possible host reads
    potential_host_indexes = []
    for i in range(0, len(reads)):
        read = reads[i]
        is_host = find_host_status(read.get_alignments(), tax2category)
        if is_host:
            potential_host_indexes.append(i)
        else:
            read.potential_host = False
    
    # 2.  Mark / delete reads said to be potential hosts
    if delete_reads:
        for i in reversed(potential_host_indexes):
            reads.pop(i)
    else:
        for i in potential_host_indexes:
            reads[i].potential_host = True
    return reads

def is_best_score_host(read_alignments, tax2category, potential_hosts):
    largest_score_alignment = max(read_alignments, key=attrgetter('score'))
    category = tax2category.get(largest_score_alignment.tax_id, None)
    if category is None or category in potential_hosts:
        return True
    else: return False

def perc_of_host_alignments_larger_than(read_alignments, tax2category, potential_hosts, percentage=0.5):
    host_aln_count = 0
    aln_count = 0
    for read_alignment in read_alignments:
        category = tax2category.get(read_alignment.tax_id, None)
        if category is None or category in potential_hosts:
            host_aln_count += 1
        aln_count += 1

    if float(host_aln_count)/aln_count >= percentage:
        return True
    return False

def are_all_alignments_host(read_alignments, tax2category, potential_hosts):
    if perc_of_host_alignments_larger_than(read_alignments, tax2category, potential_hosts, percentage=1.0):
        return True
    return False

 
    


def determine_filtering_method(delete_host):
    if delete_host:
        def filter_alignment(x_list, new_x, status):
            if not status:
                new_x.potential_host = False
                x_list.append(new_x)
            # else do nothing...
    else:
        def filter_alignment(x_list, new_x, status):
            new_x.potential_host = status
            x_list.append(new_x)

    return filter_alignment
