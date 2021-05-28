import pysam, os
from collections import Counter
import numpy as np


def get_ref_coverage_array(bamfilepath, query_positions, cache_dict):
    '''
    Get read counts for a given position or positions if range is > 1
    :param samfile: PySAM connection to BAM file
    :param query_positions: query nucleotide positions in the pileup
    :param cache_site_coverages_dict: cached coverages information for all or some positions (speed)
    :return: dictionary of positions and read counts
    '''
    coverages_list = [] #dict.fromkeys(query_positions)

    print("INFO: Getting coverage information from {}".format(bamfilepath))
    index_file = os.path.abspath(bamfilepath) + '.bai'
    if not os.path.exists(index_file):
        print("Creating bam index for {}".format(bamfilepath))
        pysam.index(bamfilepath, index_file)
    bamfile = pysam.AlignmentFile(bamfilepath, "rb")

    for position in query_positions:
        if bamfilepath in cache_dict.keys() and cache_dict[bamfilepath].get(str(position)):
            coverages_list.append(cache_dict[bamfilepath].get(str(position)))
            #print("Using cache for position {}".format(position))
        else:
            #print("Getting coverage for position {}".format(position))
            coverage_tuple = bamfile.count_coverage(bamfile.get_reference_name(0),
                                                   quality_threshold=1, start=position-1,stop=position,
                                                   read_callback="nofilter")
            coverages_list.append(np.sum(coverage_tuple))

    return coverages_list

