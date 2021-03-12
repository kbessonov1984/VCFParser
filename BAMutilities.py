import pysam, os
from collections import Counter
import numpy as np


def get_ref_coverage_array(bamfilepath, query_positions):
    '''
    Get read counts for a given position or positions if range is > 1
    :param samfile: PySAM connection to BAM file
    :param query_positions: query nucleotide positions in the pileup

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
        coverage_tuple = bamfile.count_coverage(bamfile.get_reference_name(0),
                                               quality_threshold=1, start=position-1,stop=position,
                                               read_callback="nofilter")

        coverages_list.append(np.sum(coverage_tuple))
    return coverages_list

