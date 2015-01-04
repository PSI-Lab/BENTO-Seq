import numpy as np
import pysam, logging
from .read_distribution import ReadDistribution
from .bootstrap import gen_pdf
from . import BENTOSeqError

class AltSpliceEvent(object):
    """ This class represents an alternative splicing event, defined
    by a name, the type of alternative splicing event, chromosome,
    strand, and the involved exons.

    **Parameters:**

    event_type : {'CAS', 'A5SS', 'A3SS', 'MXE', 'AFE', 'ALE', 'SPR'}

    event_id : string

    chromosome : string
        The chromosome where the event is located. Must match the
        chromosome names in the BAM-file.

    strand : {'+', '-'}

    exons : array_like
        A list/tuple of exons in the format
        ``[(exon_1_start, exon_1_end), (exon_2_start, exon_2_end), ... ]``.
        The number of exons must be three, except when the event type
        is 'MXE' when the number of exons must be four. By default,
        intervals are assumed to be in Python format, *i.e.* 0-based
        and right-open. Use the ``one_based_pos`` option if your
        coordinates are 1-based and closed.

    one_based_pos : bool (default=False)
        Whether to interpret the exon coordinates as 1-based and
        closed intervals. This is often the default format for genomic
        coordinates, *e.g.* on the UCSC Genome Browser. By default,
        coordiantes are assumed to be 0-based and right-open. If this
        option is True, then exon coordinates are converted to 0-based
        indexing on initialization. All internal calculations are
        based on 0-based coordinates.
    
    """

    
    def __init__(self, event_type, event_id, chromosome, strand, exons, one_based_pos=False):
        event_type = event_type.upper()
        
        if event_type not in ('CAS', 'A5SS', 'A3SS', 'MXE', 'AFE', 'ALE', 'SPR'):
            raise BENTOSeqError("Unknown alternative splicing event type: %s" % str(event_type))

        if event_type == 'MXE' and len(exons) != 4:
            raise BENTOSeqError("Incorrect number of exons: %d (must be 4)." % len(exons))
        elif event_type != 'MXE' and len(exons) != 3:
            raise BENTOSeqError("Incorrect number of exons: %d (must be 3)." % len(exons))

        for i, exon in enumerate(exons):
            if len(exon) != 2:
                raise BENTOSeqError("Exon #%d has wrong length: %d (must be 2)." % len(exon))

        if strand not in ('-', '+'):
            raise BENTOSeqError("Unknown strand type: %s (must be '+' or '-')." % str(strand))

        self.event_type = event_type
        self.event_id = event_id
        self.chromosome = chromosome
        self.strand = strand
        self.exons = exons

        if one_based_pos:
            self.exons = [(e[0] - 1, e[1]) for e in self.exons]

        if strand == '-':
            self.exons = [(-e[1] + 1, -e[0] + 1) for e in self.exons]

        if not ((self.exons[1][0] > self.exons[0][0] or
                 self.exons[1][1] > self.exons[0][1]) and
                (self.exons[2][0] > self.exons[1][0] or
                 self.exons[2][1] > self.exons[1][1])):
            raise BENTOSeqError("Exons must be listed from 5' to 3' "
                                "on the transcribed strand.")

        if strand == '+':
            self.junctions = [(self.chromosome, self.exons[0][1], self.exons[1][0]),
                              (self.chromosome, self.exons[1][1], self.exons[2][0]),
                              (self.chromosome, self.exons[0][1], self.exons[2][0])]
        elif strand == '-':
            self.junctions = [(self.chromosome, -self.exons[1][0] + 1, -self.exons[0][1] + 1),
                              (self.chromosome, -self.exons[2][0] + 1, -self.exons[1][1] + 1),
                              (self.chromosome, -self.exons[2][0] + 1, -self.exons[0][1] + 1)]
        else:
            raise BENTOSeqError

        self.exons_lengths = [e[1] - e[0] for e in self.exons]

        if self.event_type == 'CAS':
            if not (self.exons[1][0] > self.exons[0][1] and
                    self.exons[2][0] > self.exons[1][1]):
                raise BENTOSeqError("Event is not a valid CAS event.")
        elif self.event_type == 'A5SS':
            if not (self.exons[0][0] == self.exons[1][0] and
                    self.exons[2][0] > self.exons[1][1]):
                raise BENTOSeqError("Event is not a valid A5SS event.")
        elif self.event_type == 'A3SS':
            if not (self.exons[1][1] == self.exons[2][1] and
                    self.exons[1][0] > self.exons[0][1]):
                raise BENTOSeqError("Event is not a valid A3SS event.")
        elif self.event_type == 'MXE':
            if strand == '+':
                self.junctions.append(
                    (self.chromosome, self.exons[1][1], self.exons[3][0]))
                self.junctions.append(
                    (self.chromosome, self.exons[2][1], self.exons[3][0]))
            elif strand == '-':
                self.junctions.append(
                    (self.chromosome, -self.exons[3][0] + 1, -self.exons[1][1] + 1))
                self.junctions.append(
                    (self.chromosome, -self.exons[3][0] + 1, -self.exons[2][1] + 1))

            if not(self.exons[1][0] > self.exons[0][1] and
                   self.exons[2][0] > self.exons[0][1] and
                   self.exons[3][0] > self.exons[2][1]):
                raise BENTOSeqError("Event is not a valid MEX event.")
            
        elif self.event_type == 'AFE':
            if not (self.exons[2][0] > self.exons[0][1] and
                    self.exons[2][0] > self.exons[1][1]):
                raise BENTOSeqError("Event is not a valid AFE event.")
        elif self.event_type == 'ALE':
            if not (self.exons[1][0] > self.exons[0][1] and
                    self.exons[2][0] > self.exons[0][1]):
                raise BENTOSeqError("Event is not a valid ALE event.")
        elif self.event_type == 'SPR':
            if not (self.exons[1][0] > self.exons[0][1] and
                    self.exons[2][0] > self.exons[0][1]):
                raise BENTOSeqError("Event is not a valid SPR event.")
        else:
            raise BENTOSeqError

    def build_read_distribution(self, bamfile, min_overhang=5,
                                max_edit_distance=2,
                                max_num_mapped_loci=1):

        """Build the read distribution for this event from a BAM-file.

        **Parameters:**

        bamfile : :py:class:`pysam.Samfile`
            Reference to a binary, sorted, and indexed SAM-file.

        min_overhang : int (default=5)
            Minimum overhang on either side of the splice junction
            required for counting a read towards the read distribution.

        max_edit_distance : int (default=2)
            Maximum edit distance (number of mismatches against the
            reference genome) to allow before skipping a read.

        max_num_mapped_loci : int (default=1)
            Indiciates to how many locations a read may be aligned to
            be a counted. By default, only uniquely mappable reads are
            alowed.

        """
    
        self.junction_read_distributions = []
        for junction in self.junctions:
            read_distribution = \
                ReadDistribution.from_junction(
                    bamfile, junction,
                    max_edit_distance,
                    max_num_mapped_loci)

            if read_distribution.is_empty:
                logging.debug("Event %s: No reads in %s "
                              "map to junction %s:%d:%d." %
                              (self.event_id, bamfile.filename,
                               junction[0], junction[1], junction[2]))

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
            raise BENTOSeqError

    def trim_reads(self, reads, read_length,
                   min_overhang, exon_length):
        """Avoid double-counting reads spanning two splice junctions. """

        event_type = self.event_type

        if event_type in ('CAS', 'MXE'):
            trim = read_length - 2 * min_overhang - exon_length
        elif event_type in ('AFE', 'ALE'):
            trim = read_length - min_overhang - exon_length
        else:
            raise BENTOSeqError
        
        if trim > 0:
            if event_type in ('CAS', 'MXE', 'ALE'):        
                reads = reads[:-trim]
            elif event_type == 'AFE':
                reads = reads[-trim:]
            else:
                raise BENTOSeqError
        return reads
        
    def bootstrap_event(self, n_bootstrap_samples=1000, n_grid_points=100,
                        a=1, b=1, r=0):

        """Estimate PSI (percent spliced-in) value for this event.

        **Parameters:**

        n_bootstrap_samples : int (default=1000)
            How many bootstrap samples to use for estimation of PSI.

        n_grid_points : 100 (default=100)
            How many points to use for the numerical approximation of
            the bootstrap probability density function.

        a : int (default=1)
            Bayesian pseudo-count for the inclusion reads.

        b : int (default=1)
            Bayesian pseudo-count for the exclusion reads.

        r : int (default=0)
            Bayesian pseudo-count for the normalization of the
            bootstrap probability density function.

        **Returns:**

        n_inc : int
            The number of reads mapped to the inclusion junctions.

        n_exc : int
            The number of reads mapped to the exclusion junctions.

        p_inc : int
            Number of inclusion mapping positions.

        p_exc : int
            Number of exclusion mapping positions.

        psi_standard : float
            Naive PSI estimate.

        psi_bootstrap : float
            Bootstrap estimate of PSI.

        psi_bootstrap_std : float
            Estimated standard deviation of ``psi_bootstrap``.
        """
    
        reads_inc = np.array(self.reads_inc)
        reads_exc = np.array(self.reads_exc)

        n_inc = reads_inc.sum()
        n_exc = reads_exc.sum()

        p_inc = reads_inc.size
        p_exc = reads_exc.size

        min_p = min(p_inc, p_exc)
        scaled_inc = float(n_inc) / p_inc * min_p
        scaled_exc = float(n_exc) / p_exc * min_p

        psi_standard = (scaled_inc + 1) / (scaled_inc + scaled_exc + 2)

        pdf, grid = gen_pdf(reads_inc, reads_exc,
                            n_bootstrap_samples, n_grid_points, a, b, r)

        psi_bootstrap = np.sum(pdf * grid)
        psi_std = np.sqrt(np.sum(pdf * np.square(grid - psi_bootstrap)))

        return n_inc, n_exc, p_inc, p_exc, psi_standard, psi_bootstrap, psi_std
