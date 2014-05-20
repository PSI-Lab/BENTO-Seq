#!/usr/bin/env python

import sys, argparse, pysam, logging, datetime
from bs_psi.alt_splice_event import AltSpliceEvent

# def _warning(
#     message,
#     category = UserWarning,
#     filename = '',
#     lineno = -1):
#     """Make warnings a little cleaner."""
#     print(message)

# warnings.showwarning = _warning

def run_bootstrap():
    parser = argparse.ArgumentParser()
    parser.add_argument('event_definitions',
                        help="Alternative splicing event definitions file. "
                        "See http://github.com/xxxx for an example and "
                        "formatting instructions.")
    parser.add_argument('bamfile',
                        help="A bam-file with aligned reads from either the "
                        "TopHat or the STAR read alignmnent tool. Must be a "
                        "binary and indexed bam-file (.bam file extension "
                        "and a .bai index file must be present in the "
                        "same directory.)")

    parser.add_argument('output_file', help="Name of the file where "
                        "the output should be stored. If the file "
                        "exists, it will be overwritten.")

    parser.add_argument('-1', '--one-based-coordinates',
                        action='store_true', help="Use this option when "
                        "the coordinates in your event file use "
                        "one-based indexing. Otherwise, standard Python "
                        "zero-based indexing is used, i.e. the first "
                        "element is 0 and intervals are "
                        "right-open. [1:4] --> (1, 2, 3).")

    parser.add_argument('-q', '--quiet',
                        action='store_true',
                        help="Suppress all warnings and messages")

    parser.add_argument('-nm', '--max-edit-distance',
                        help="(default=2) The maximum edit distance or "
                        "number of mismatches allowed for a read. If "
                        "the edit distance is greater than this, the "
                        "read will not be counted towards the read  "
                        "distribution. The edit distance is the number "
                        "of insertions, deletions, or substitutions "
                        "required to match the read to the consensus "
                        "sequence (reference genome). The edit "
                        "distance corresponds to the 'NM' tag in "
                        "bam-files created with TopHat and the 'nM' in "
                        "files created with SAM.", type=int, default=2)

    parser.add_argument('-nh', '--max-num-mapped-loci',
                        help="(default=1) The maximum number of loci "
                        "(genome locations) that a read "
                        "is allowed to be mapped to. If the number of "
                        "loci the read is mapped to is greater than "
                        "this, the the read will not be counted towards "
                        "the read distribution. By default, only "
                        "uniquely mapped reads are considered "
                        "(--max-num-mapped-loci=1). The number of "
                        "mapped loci is given by the 'NH' tag in the "
                        "bam-file.", type=int, default=1)

    parser.add_argument('-oh', '--min-overhang',
                        help="(default=5) The minimum overhang of a read "
                        "on either side of the splice junction to be "
                        "considered for the read distribution. By default "
                        "a read must have a 5nt overhang on either side of "
                        "the splice junction to be counted.",
                        type=int, default=5)

    parser.add_argument('-S', '--n-bootstrap-samples',
                        help="(default=1000) The number of bootstrap  "
                        "samples to draw for each event. More bootstrap "
                        "samples will better estimate the bootstrap "
                        "probability density function, but will need "
                        "more memory and will take longer to compute.",
                        type=int, default=1000)

    parser.add_argument('-G', '--n-grid-points',
                        help="(default=100) The number of points "
                        "used for the numerical approximation of the "
                        "bootstrap probability density function. More "
                        "points will better estimate the bootstrap "
                        "probability density function, but will need "
                        "more memory and will take longer to compute.",
                        type=int, default=100)

    parser.add_argument('-a', help="(default=1) Bayesian pseudo-count "
                        "for inclusion reads. See http://github.com/xxx "
                        "for details.", type=int, default=1)
    
    parser.add_argument('-b', help="(default=1) Bayesian pseudo-count "
                        "for exclusion reads. See http://github.com/xxx "
                        "for details.", type=int, default=1)

    parser.add_argument('-r', help="(default=0) Bayesian pseudo-count "
                        "for numerical integration of bootstrap "
                        "probability density function. See "
                        "http://github.com/xxx for details.", type=int,
                        default=1)

    args = parser.parse_args()
    if not args.quiet:
        logging.basicConfig(level='INFO')
    else:
        logging.basicConfig(level='ERROR')
    process_event_file(args)

def process_event_file(args):
    start_t = datetime.datetime.now()
    with open(args.event_definitions) as f:
        output_file = open(args.output_file, 'w')
        bamfile = pysam.Samfile(args.bamfile, check_header=False)

        # Write header
        output_file.write(
        '\t'.join(('#ID', 'n_inc', 'n_exc', 'p_inc', 'p_exc', 'PSI_standard',
                   'PSI_bootstrap', 'PSI_bootstrap_std')) + '\n')
    
        for n_events, line in enumerate(f):
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
        logging.info("Output written to file '%s'." % args.output_file)
        runtime = datetime.datetime.now() - start_t
        logging.info("Processed %d events in %.2f seconds." % (n_events, runtime.total_seconds()))
            
if __name__ == '__main__':
    sys.exit(run_bootstrap())
