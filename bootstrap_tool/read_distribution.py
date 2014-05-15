from collections import OrderedDict, Counter

class ReadDistribution(object):
    __slots__ = ['_counter', 'read_length']

    def __init__(self, read_length, read_distribution=None):
        self.read_length = read_length
        self._counter = Counter(read_distribution)

    def to_list(self, min_overhang=0):
        start = -(self.read_length - min_overhang)
        end = -min_overhang
        return [(pos, self._counter[pos]) for pos in range(start, end)]

    def __getitem__(self, key):
        return self._counter[key]

    def __setitem__(self, key, value):
        self._counter[key] = value


class ReadDistributionCollection(object):
    __slots__ = ['_dict', 'read_length']
    
    def __init__(self, read_length):
        self._dict = OrderedDict()
        self.read_length = read_length

    def __getitem__(self, key):
        return self._dict[key]

    def __setitem__(self, key, value):
        self._dict[key] = ReadDistribution(self.read_length)
        for pos, reads in zip(value):
            self._dict[key][pos] = reads

    def inc(self, key, pos, value=1):
        if not key in self._dict:
            self._dict[key] = ReadDistribution(self.read_length)
        self._dict[key][pos] += value

    def sort(self):
        self._dict = OrderedDict([(k, self._dict[k]) for k in sorted(self._dict.keys())])

    def to_file(self):
        yield "#read_length:%d" % read_length
        yield "#junction_ID\tpositions\tread_distribution\n"
        for key, value in self._dict.iteritems():
            key_str = '%d:%d:%d' % key
            positions, reads = zip(*value.to_list(0))
            positions_str = ','.join(map(str, positions))
            reads_str = ','.join(map(str, reads))
            yield '\t'.join((key_str, positions_str, reads_str)) + '\n'

    @classmethod
    def from_file(cls, iterable):
        read_length = int(iterable.next().rstrip().split(':'))[1])
        read_distribution = cls(read_length)
        for line in iterable:
            if line.startswith('#'): continue
            key_str, position_string, reads_str = line.rstrip().split('\t')
            key = tuple(key_str.split(':'))
            positions = map(int, positions_str.split(','))
            reads = map(int, reads_str.split(','))
            read_distribution[key] = zip(positions, reads)
        return read_distribution

    def get_positions(self, min_overhang):
        """Get all possible mapping positions for a read length and
        minimum overhang.

        """
        
        start = -(self.read_length - min_overhang)
        end = -min_overhang
        return range(start, end)