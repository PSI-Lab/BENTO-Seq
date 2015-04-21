import re
from copy import copy
from collections import Counter

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

class ReadDistribution(object):
    """This class represents a distribution of reads across a splice junction.

    **Parameters:**

    chromosome : string

    start : int
        Start of the splice junction, must be 0-based coordinate.

    end : int
        End of the splice junction, must be 0-based coordinate.

    read_length : int

    read_distribution : dict (optional)

        Read distribution as a dictionary in the format ``{pos1:
        counts1, pos2: counts2, ...}``, where the keys are position
        relative to the splice junction, *i.e.* ``-10`` would be the
        position 10nt upstream of the splice junction.

    """
    
    __slots__ = ['_counter', '_read_ids', 'chromosome', 'start', 'end', 'read_length']

    def __init__(self, chromosome, start, end, read_length, read_distribution=None):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.read_length = read_length
        self._counter = Counter(read_distribution)
        self._read_ids = {}

    @property
    def is_empty(self):
        """Returns whether the read distribution contains any reads.
        """
    
        return not bool(self._counter)

    def to_list(self, min_overhang=0):
        """Return the read distribution as a list given the minimum overhang.

        **Parameters:**

        min_overhang : int (default=0)

            The minimum overhang for which to compute the read
            distribution. *E.g.* if the read length is 75 and
            ``min_overhang=5``, then the read distribution is computed
            for positions ``[-70, -69, ..., -5]``.

        **Returns:**

        read_distribution : list
            Read distribution in the format ``[(pos1, counts1), (pos2,
            counts2), ...]``.
        """
    
        return [(pos, self._counter[pos]) for pos in self.get_positions(min_overhang)]

    def __getitem__(self, pos):
        """Return the number of reads at ``pos``.
        """
    
        return self._counter[pos]

    def __setitem__(self, pos, value):
        self._counter[pos] = value

    def get_positions(self, min_overhang):
        """Get all possible mapping positions for a read length and
        minimum overhang.

        """
        start = -(self.read_length - min_overhang)
        end = - min_overhang + 1
        return range(start, end)

    def inc(self, rel_pos, read):
        """Add a read to the distribution and keep a reference to it.

        **Parameters:**

        rel_pos : int
            Mapping position of the read relative to the splice junction.

        read : :py:class:`pysam.AlignedRead`
            Reference to the read.

        """
    
        if rel_pos not in self._read_ids:
            self._read_ids[rel_pos] = []
        self._read_ids[rel_pos].append(read)
        self._counter[rel_pos] += 1

    @classmethod
    def from_junction(cls, bamfiles, junction,
                      max_edit_distance=2,
                      max_num_mapped_loci=1):
        """Build the read distribution from a BAM-file.

        **Parameters:**

        bamfile : list of :py:class:`pysam.Samfile`

        junction : tuple
            Tuple in the format ``(chromosome, start, end)``.

        max_edit_distance : int (default=2)

        max_num_mapped_loci : int (default=1)

        **Returns:**

        read_distribution : :class:`bs_psi.read_distribution.ReadDistribution`    

        """


        if not isinstance(bamfiles, (list, tuple)):
            bamfiles = [bamfiles]
        chromosome, junction_start, junction_end = junction
        read_length = bamfiles[0].next().rlen
        read_distribution = cls(chromosome, junction_start, junction_end, read_length)

        for bamfile in bamfiles:
            for read in bamfile.fetch(chromosome, junction_start, junction_start + 1):
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

                read_distribution.inc(rel_pos, read)

        return read_distribution
