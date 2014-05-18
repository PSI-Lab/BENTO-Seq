from collections import Counter

class ReadDistribution(object):
    __slots__ = ['_counter', '_read_ids', 'chromosome', 'start', 'end', 'read_length']

    def __init__(self, chromosome, start, end, read_length, read_distribution=None):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.read_length = read_length
        self._counter = Counter(read_distribution)
        self._read_ids = {}

    def to_list(self, min_overhang=0):
        return [(pos, self._counter[pos]) for pos in self.get_positions(min_overhang)]

    def __getitem__(self, key):
        return self._counter[key]

    def __setitem__(self, key, value):
        self._counter[key] = value

    def get_positions(self, min_overhang):
        """Get all possible mapping positions for a read length and
        minimum overhang.

        """
        start = -(self.read_length - min_overhang)
        end = - min_overhang + 1
        return range(start, end)

    def inc(self, rel_pos, read_id):
        if rel_pos not in self._read_ids:
            self._read_ids[rel_pos] = []
        self._read_ids[rel_pos].append(read_id)
        self._counter[rel_pos] += 1
