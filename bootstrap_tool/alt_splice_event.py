import numpy as np
import pysam
from .bamfile_parsing import read_distribution_for_junction
from .bootstrap import gen_pdf

class AltSpliceEvent(object):
    def __init__(self, event_type, event_id, chromosome, strand, exons, one_based_pos=False):
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

        if one_based_pos:
            self.exons = [(e[0] - 1, e[1]) for e in self.exons]

        if strand == '-':
            self.exons = [(-e[0], -e[1]) for e in self.exons]

        # Swap exon end and start if in wrong order
        self.exons = [(e[1], e[0]) if e[0] > e[1] else e for e in self.exons]

        for i in range(len(self.exons) - 1):
            if not ((self.exons[1][0] > self.exons[0][0] or
                     self.exons[1][1] > self.exons[0][1]) and
                    (self.exons[2][0] > self.exons[1][0] or
                     self.exons[2][1] > self.exons[1][1])):
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

    def build_read_distribution(self, samfile, min_overhang=5,
                                max_edit_distance=2,
                                max_num_mapped_loci=1):
        self.junction_read_distributions = []
        for junction in self.junctions:
            read_distribution = \
                read_distribution_for_junction(
                    samfile, junction,
                    max_edit_distance,
                    max_num_mapped_loci)
            read_length = read_distribution.read_length

            positions, reads = zip(
                *read_distribution.to_list(min_overhang))
            if self.strand == '-': reads = reads[::-1]
            self.junction_read_distributions.append(reads)

        if self.event_type == 'CAS':
            self.junction_read_distributions[0] = self.trim_reads(
                self.junction_read_distributions[0], read_length,
                min_overhang, self.exons_lengths[1])
            
            self.reads_inc = self.junction_read_distributions[0] + \
                              self.junction_read_distributions[1]
            self.reads_exc = self.junction_read_distributions[2]
            
        elif self.event_type == 'A5SS':
            self.reads_inc = self.junction_read_distributions[1]
            self.reads_exc = self.junction_read_distributions[2]
            
        elif self.event_type == 'A3SS':
            self.reads_inc = self.junction_read_distributions[0]
            self.reads_exc = self.junction_read_distributions[2]
            
        elif self.event_type == 'MXE':
            self.junction_read_distributions[0] = self.trim_reads(
                self.junction_read_distributions[0], read_length,
                min_overhang, self.exons_lengths[1])
            self.junction_read_distributions[2] = self.trim_reads(
                self.junction_read_distributions[2], read_length,
                min_overhang, self.exons_lengths[2])
            self.reads_inc = self.junction_read_distributions[0] + \
                              self.junction_read_distributions[3]
            self.reads_exc = self.junction_read_distributions[1] + \
                              self.junction_read_distributions[4]
            
        elif self.event_type == 'AFE':
            self.junction_read_distributions[2] = self.trim_reads(
                self.junction_read_distributions[2], read_length,
                min_overhang, self.exons_lengths[0])
            self.junction_read_distributions[1] = self.trim_reads(
                self.junction_read_distributions[1], read_length,
                min_overhang, self.exons_lengths[1])
            
            self.reads_inc = self.junction_read_distributions[1]
            self.reads_exc = self.junction_read_distributions[2]
            
        elif self.event_type == 'ALE':
            self.junction_read_distributions[0] = self.trim_reads(
                self.junction_read_distributions[0], read_length,
                min_overhang, self.exons_lengths[1])
            self.junction_read_distributions[2] = self.trim_reads(
                self.junction_read_distributions[2], read_length,
                min_overhang, self.exons_lengths[2])
            
            self.reads_inc = self.junction_read_distributions[0]
            self.reads_exc = self.junction_read_distributions[2]
            
        elif self.event_type == 'SPR':
            self.reads_inc = self.junction_read_distributions[0]
            self.reads_exc = self.junction_read_distributions[2]
            
        else:
            raise ValueError

    def trim_reads(self, reads, read_length,
                   min_overhang, exon_length):
        """Avoid double-counting reads spanning two splice junctions. """

        event_type = self.event_type

        if event_type in ('CAS', 'MXE'):
            trim = read_length - 2 * min_overhang - exon_length
        elif event_type in ('AFE', 'ALE'):
            trim = read_length - min_overhang - exon_length
        else:
            raise ValueError
        
        if trim > 0:
            if event_type in ('CAS', 'MXE', 'ALE'):        
                reads = reads[:-trim]
            elif event_type == 'AFE':
                reads = reads[-trim:]
            else:
                raise ValueError
        return reads
        
    def bootstrap_event(self, n_bootstrap_samples=1000, n_grid_points=100,
                        a=1, b=1, r=0):
        reads_inc = np.array(self.reads_inc)
        reads_exc = np.array(self.reads_exc)

        n_inc = reads_inc.sum()
        n_exc = reads_exc.sum()

        p_inc = reads_inc.size
        p_exc = reads_exc.size

        min_p = min(p_inc, p_exc)
        scaled_inc = n_inc / p_inc * min_p
        scaled_exc = n_exc / p_exc * min_p

        psi_standard = (scaled_inc + 1.) / (scaled_inc + scaled_exc + 2)

        pdf, grid = gen_pdf(reads_inc, reads_exc,
                            n_bootstrap_samples, n_grid_points, a, b, r)

        psi_bootstrap = np.sum(pdf * grid)
        psi_std = np.sqrt(np.sum(pdf * np.square(grid - psi_bootstrap)))

        return n_inc, n_exc, p_inc, p_exc, psi_standard, psi_bootstrap, psi_std

def process_event_file(args):
    with open(args.event_definitions) as f:
        output_file = open(args.output_file, 'w')
        bamfile = pysam.Samfile(args.bamfile, check_header=False)

        # Write header
        output_file.write(
        '\t'.join(('#ID', 'n_inc', 'n_exc', 'p_inc', 'p_exc', 'PSI_standard',
                   'PSI_bootstrap', 'PSI_bootstrap_std')) + '\n')
    
        for line in f:
            line = line.rstrip()
            if line.startswith('#') or not line: continue
            elements = line.split('\t')
            event_type, event_id, chromosome, strand = elements[:4]
            if event_type.upper() is not 'MXE':
                exons = [tuple(map(int, e.split(':'))) for e in elements[4:7]]
            else:
                exons = [tuple(map(int, e.split(':'))) for e in elements[4:8]]
            event = AltSpliceEvent(event_type, event_id, chromosome,
                                   strand, exons,
                                   one_based_pos=args.one_based_coordinates)
            event.build_read_distribution(bamfile, args.min_overhang,
                                          args.max_edit_distance,
                                          args.max_num_mapped_loci)
            psi_event = event.bootstrap_event(args.n_bootstrap_samples,
                                              args.n_grid_points,
                                              args.a, args.b, args.r)
            output_file.write(
                '\t'.join([event.event_id] + map(str, psi_event)) + '\n'
            )
            

        output_file.close()
            
