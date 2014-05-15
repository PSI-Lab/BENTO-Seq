class AltSpliceEvent(object):
    def __init__(self, event_type, event_id, chromosome, strand, exons):
        event_type = event_type.upper()
        
        if event_type not in ('CAS', 'A5SS', 'A3SS', 'MXE', 'AFE', 'ALE', 'SPR'):
            raise ValueError("Unknown alternative splicing event type: %s" % str(event_type))

        if event_type == 'MXE' and len(exons) != 4:
            raise ValueError("Incorrect number of exons: %d (must be 4)." % len(exons))
        elif len(exons) != 3:
            raise ValueError("Incorrect number of exons: %d (must be 3)." % len(exons))

        for i, exon in enumerate(exons):
            if len(exon) != 2:
                raise ValueError("Exon #%d has wrong length: %d (must be 2)." % len(exon))

        if strand not in ('-', '+'):
            raise ValueError("Unknown strand type: %s (must be '+' or '-')." % str(strand))

        self.event_type = event_type
        self.event_id = event_id
        self.chromosome = chromosome
        self.strand = strand
        self.exons = exons

        if strand == '-':
            self.exons = [(-e[0], -e[1]) for e in self.exons]

        # Swap exon end and start if in wrong order
        self.exons = [(e[1], e[0]) if e[0] > e[1] else e for e in self.exons]

        for i in range(len(self.exons) - 1):
            if self.exons[i][1] > self.exons[i + 1][0]:
                raise ValueError("Exons must be listed from 5' to 3' "
                                 "on the transcribed strand.")

        if strand == '+':
            self.junctions = [(self.chromosome, self.exons[0][1], self.exons[1][0]),
                              (self.chromosome, self.exons[1][1], self.exons[2][0]),
                              (self.chromosome, self.exons[0][1], self.exons[2][0])]
        elif strand == '-':
            self.junctions = [(self.chromosome, -self.exons[1][0], -self.exons[0][1]),
                              (self.chromosome, -self.exons[2][0], -self.exons[1][1]),
                              (self.chromosome, -self.exons[2][0], -self.exons[0][1])]
        else:
            raise ValueError

        self.exons_lengths = [e[1] - e[0] for e in self.exons]

        if self.event_type == 'CAS':
            if not (self.exons[1][0] > self.exons[0][1] and
                    self.exons[2][0] > self.exons[1][1]):
                raise ValueError("Event is not a valid CAS event.")
        elif self.event_type == 'A5SS':
            if not (self.exons[0][0] == self.exons[1][0] and
                    self.exons[2][0] > self.exons[1][1]):
                raise ValueError("Event is not a valid A5SS event.")
        elif self.event_type == 'A3SS':
            if not (self.exons[1][1] == self.exons[2][1] and
                    self.exons[1][0] > self.exons[0][1]):
                raise ValueError("Event is not a valid A3SS event.")
        elif self.event_type == 'MXE':
            if strand == '+':
                self.junctions.append(
                    (self.chromosome, self.exons[1][1], self.exons[3][0]))
                self.junctions.append(
                    (self.chromosome, self.exons[2][1], self.exons[3][0]))
            elif strand == '-':
                self.junctions.append(
                    (self.chromosome, -self.exons[3][0], -self.exons[1][1]))
                self.junctions.append(
                    (self.chromosome, -self.exons[3][0], -self.exons[2][1]))

            if not(self.exons[1][0] > self.exons[0][1] and
                   self.exons[2][0] > self.exons[1][1] and
                   self.exons[3][0] > self.exons[2][1]):
                raise ValueError("Event is not a valid MEX event.")
            
        elif self.event_type == 'AFE':
            if not (self.exons[2][0] > self.exons[0][1] and
                    self.exons[2][0] > self.exons[1][1]):
                raise ValueError("Event is not a valid AFE event.")
        elif self.event_type == 'ALE':
            if not (self.exons[1][0] > self.exons[2][1] and
                    self.exons[2][0] > self.exons[0][1]):
                raise ValueError("Event is not a valid ALE event.")
        elif self.event_type == 'SPR':
            if not (self.exons[1][0] > self.exons[0][1] and
                    self.exons[2][0] > self.exons[0][1]):
                raise ValueError("Event is not a valid SPR event.")
        else:
            raise ValueError

    def build_read_distribution(self, read_distributions, min_overhang=5):
        positions = read_distributions.get_positions(min_overhang)
        n_positions = len(positions)
        read_length = read_distributions.read_length

        self.junction_read_distributions = []
        for junction in self.junctions:
            positions, reads = zip(
                *read_distributions[junction].to_list(min_overhang)))
            if strand == '-': reads = reads[::-1]
            self.junction_read_distributions.append(reads)

        if self.event_type == 'CAS':
            reads[0] = self.trim_reads(
                reads[0], read_length, min_overhang, self.exons_lengths[1], True)
            
            self.reads_inc = (reads[0], reads[1])
            self.reads_exc = (reads[2],)
        elif self.event_type == 'A5SS':
            self.reads_inc = (reads[1],)
            self.reads_exc = (reads[2],)
        elif self.event_type == 'A3SS':
            self.reads_inc = (reads[0],)
            self.reads_exc = (reads[2],)
        elif self.event_type == 'MXE':
            reads[0] = self.trim_reads(
                reads[0], read_length, min_overhang, self.exons_lengths[1], True)
            reads[2] = self.trim_reads(
                reads[2], read_length, min_overhang, self.exons_lengths[2], True)
            self.reads_inc = (reads[0], reads[3])
            self.reads_exc = (reads[1], reads[4])
        elif self.event_type == 'AFE':
            reads[2] = self.trim_reads(
                reads[2], read_length, min_overhang, self.exons_lengths[0], False)
            reads[1] = self.trim_reads(
                reads[1], read_length, min_overhang, self.exons_lengths[1], False)
            
            self.reads_inc = (reads[1],)
            self.reads_exc = (reads[2],)
        elif self.event_type == 'ALE':
            reads[0] = self.trim_reads(
                reads[0], read_length, min_overhang, self.exons_lengths[1], True)
            reads[2] = self.trim_reads(
                reads[2], read_length, min_overhang, self.exons_lengths[2], True)
            
            self.reads_inc = (reads[0],)
            self.reads_exc = (reads[2],)
        elif self.event_type == 'SPR':
            self.reads_inc = (reads[0],)
            self.reads_exc = (reads[2],)
        else:
            raise ValueError

    def trim_reads(self, reads, read_length, min_overhang, exon_length, from_end=False):
        """Avoid double-counting reads spanning two splice junctions. """

        trim = read_length - min_overhang - exon_length
        if trim > 0:
            if from_end:
                reads = reads[:-trim]
            else:
                reads = reads[-trim:]
        return reads

        