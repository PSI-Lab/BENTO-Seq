import pysam, re
import numpy as np
from copy import copy
from collections import Counter, OrderedDict
from read_distribution import ReadDistributionCollection

has_junction = re.compile(r'(\d+)N')                # Read has junction
soft_clipping_left = re.compile(r'^(\d+)S(\d+)M')   # Read has soft-clipping on left side
soft_clipping_right = re.compile(r'(\d+)M(\d+)S$')  # Read has soft-clipping on right side
indel = re.compile(r'(\d+)[ID]')                    # Read contains indels
indel_at_ss_left = re.compile(r'(\d+)[ID](\d+)N')   # Read has indel left of SS
indel_at_ss_right = re.compile(r'(\d+)N(\d+)[ID]')  # Read has indel right of SS 
indel_right = re.compile(r'(\d+)M(\d+)([ID])')      # Indel right of match
merge_cigar = re.compile(r'(\d+)M(\d+)M')           # Merge consecutive matches
find_junctions = re.compile(r'(\d+)M(\d+)N')        # Find splice junctions
last_match = re.compile(r'(\d+)M$')                 # Find last aligned segment

def filter_junction_reads(read, max_edit_distance=2, max_mismatches=1):
    # Skip reads without junctions
    if not has_junction.search(read.cigarstring): return None

    # Extract [nN]M tag
    tags = {key: value for key, value in read.tags}
    if 'NM' in tags:
        mapper = 'TopHat'
        edit_distance = tags['NM']
    elif 'nM' in tags:
        mapper = 'STAR'
        edit_distance = tags['nM']
    else:
        raise ValueError("Incompatible BAM/SAM format: "
                         "optional TAG [Nn]M is not present.")

    # Skip if number of mismatches is greater than allowed
    if tags['NH'] > max_mismatches: return None

    cigar = copy(read.cigarstring)

    # Count soft clipping towards the edit distance
    m = soft_clipping_left.search(cigar)
    if m:
        edit_distance += int(m.groups()[0])
        tmp = sum(map(int, m.groups()))
        cigar = soft_clipping_left.sub('%dM' % tmp, cigar)

    m = soft_clipping_right.search(cigar)
    if m:
        edit_distance += int(m.groups()[1])
        tmp = sum(map(int, m.groups()))
        cigar = soft_clipping_right.sub('%dM' % tmp, cigar)

    # Count indels for STAR input
    if mapper == 'STAR':
        for m in indel.finditer(cigar):
            edit_distance += int(m.groups()[0])

    # Skip if edit distance greater than allowed
    if edit_distance > max_edit_distance: return None

    # Skip if there are indels right at the splice junction
    if indel_at_ss_left.search(cigar) or indel_at_ss_right.search(cigar):
        return None

    # Pre-process indels to properly find read positions relative to junctions
    if indel.search(cigar):
        m = indel_right.search(cigar)
        while m:
            indel_type = m.groups()[2]

            if indel_type == 'I':
                tmp = int(m.groups()[0])
            elif indel_type == 'D':
                tmp = sum(map(int, m.groups()[:2]))
            else:
                raise ValueError
            cigar = indel_right.sub('%dM' % tmp, cigar)
            m = indel_right.search(cigar)

        m = merge_cigar.search(cigar)
        while m:
            tmp = sum(map(int, m.groups()))
            cigar = merge_cigar.sub('%dM' % tmp, cigar)
            m = merge_cigar.search(cigar)

    assert not merge_cigar.search(cigar)

    strand = '-' if read.is_reverse else '+'
    pos = [read.pos]
    n_junc = 0
    for m in find_junctions.finditer(cigar):
        n_junc += 1
        len_match = int(m.groups()[0])
        len_junction = int(m.groups()[1])
        end_prev = pos[-1] + len_match   # End of the preceding aligned segment
        start_next = end_prev + len_junction   # Start of next aligned segment
        pos += [end_prev, start_next]

    m = last_match.search(cigar)
    pos += [pos[-1] + int(m.groups()[0])]

    junction_read = {
        'id': read.qname,
        'edit_distance': edit_distance,
        'n_mismatches': tags['NH'],
        'chr': read.tid,
        'strand': strand,
        'n_junc': n_junc,
        'pos': pos
    }
    return junction_read

def update_read_distribution(read_distribution, junction_read):
    pos = junction_read['pos']
    segments = [end - start for start, end in zip(pos[::2], pos[1::2])]
    for i, start, end in zip(range(junction_read['n_junc']), pos[1::2], pos[2::2]):
        junc_key = (junction_read['chr'], start, end)

        read_start_pos = -sum(segments[:(i+1)]) # Read start position relative to splice junction
        read_distribution.inc(junc_key, read_start_pos)

    return read_distribution
    
                
def process_samfile(filename):
    if filename.endswith('.bam'):
        file_flags = 'rb'
    elif filename.endswith('.sam'):
        file_flags = 'r'
    else:
        raise ValueError("Filename must have a .sam or .bam extension: %s" % filename)
    
    samfile = pysam.Samfile(filename, file_flags, check_header=False)
    read_distribution = None
    for read in samfile:
        if read_distribution is None:
            read_distribution = ReadDistributionCollection(read.rlen)
        junction_read = filter_junction_reads(read)
        if not junction_read: continue
        read_distribution = update_read_distribution(read_distribution, junction_read)

    read_distribution.sort()

    return read_distribution

if __name__ == '__main__':
    read_distribution = process_samfile('/home/hannes/git/bootstrap-tool/examples/STAR_chr21.bam')
    pass
