import pysam, re
import numpy as np
from copy import copy
from read_distribution import ReadDistribution

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

def read_distribution_for_junction(samfile, junction,
                                   max_edit_distance=2,
                                   max_num_mapped_loci=1):
    chromosome, junction_start, junction_end = junction
    read_length = samfile.next().rlen
    read_distribution = ReadDistribution(chromosome, junction_start, junction_end, read_length)
    
    for read in samfile.fetch(chromosome, junction_start, junction_start + 1):
        # Skip reads without junctions
        if not has_junction.search(read.cigarstring): continue
        
        read_junctions = [(read.blocks[i][1], read.blocks[i + 1][0])
                          for i in range(len(read.blocks) - 1)]

        if (junction_start, junction_end) not in read_junctions: continue

        pos = read.pos

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

        # Skip if the number of loci the read maps to is greater than allowed
        if tags['NH'] > max_num_mapped_loci: continue

        cigar = copy(read.cigarstring)

        # Count soft clipping towards the edit distance
        m = soft_clipping_left.search(cigar)
        if m:
            edit_distance += int(m.groups()[0])
            tmp = sum(map(int, m.groups()))
            cigar = soft_clipping_left.sub('%dM' % tmp, cigar)
            pos -= int(m.groups()[0])

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
        if edit_distance > max_edit_distance: continue

        # Skip if there are indels right at the splice junction
        if indel_at_ss_left.search(cigar) or indel_at_ss_right.search(cigar):
            continue

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
        blocks = [pos]
        n_junc = 0
        for m in find_junctions.finditer(cigar):
            n_junc += 1
            len_match = int(m.groups()[0])
            len_junction = int(m.groups()[1])
            end_prev = blocks[-1] + len_match      # End of the preceding aligned segment
            start_next = end_prev + len_junction   # Start of next aligned segment
            blocks += [end_prev, start_next]

        m = last_match.search(cigar)
        blocks += [blocks[-1] + int(m.groups()[0])]

        block_sizes = [end - start for start, end in zip(blocks[::2], blocks[1::2])]
        junction_idx = read_junctions.index((junction_start, junction_end))
        rel_pos = -sum(block_sizes[:(junction_idx + 1)])
        
        read_distribution.inc(rel_pos, read.qname)

    return read_distribution
