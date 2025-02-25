import argparse
import sys
from .scripts.run_scans import Scanner
from datetime import datetime

DNA_ALPHABET = ['A', 'T', 'G', 'C', '-', '*']


def valid_alignment(alignment):
    """
    Check that the input alignment is valid
    :param alignment: a list of lists containing the sequence headers and the aligned sequences
    :return True if the alignment is valid, false otherwise
    """
    aln_len = len(alignment[0])
    for seq in alignment:
        if len(seq) != aln_len:
            print("Improper alignment! Not all alignments are the same length.")
            return False
    return True


def valid_chars(alignment):
    """
    Check that the alignment only contains valid characters
    :param alignment: a list of lists containing the sequence headers and the aligned sequences
    :return: True if the alignment contains only valid characters, False otherwise
    """
    for s in alignment:
        if not all(pos in DNA_ALPHABET for pos in s):
            print("Alignment contains invalid characters.")
            return False
    return True


def read_fasta(handle):
    """
    Converts a FASTA formatted file to a tuple containing a list of headers and sequences
    :param handle: file stream for the FASTA file
    :return: tuple of headers and sequences
    """
    result = []
    headers, seqs = [], []
    sequence, h = '', ''

    # Verifies files have the correct formatting
    for i, line in enumerate(handle):
        if line.startswith('>'):
            break
        else:
            print("No header")
            raise NameError

    # Reset pointer to beginning of file
    if hasattr(handle, 'seek'):
        handle.seek(0)

    for line in handle:
        if line.startswith('>'):
            if len(sequence) > 0:
                headers.append(h)
                seqs.append(sequence)
                sequence = ''
            h = line.strip('>\t\n\r')
        else:
            sequence += line.strip('\n\r').upper()

    # Handle the last entry
    seqs.append(sequence)
    headers.append(h)

    return headers, seqs


def parse_args():
    parser = argparse.ArgumentParser(
        description='An open source implementation of RDP5.\n'
    )

    parser.add_argument('infile',
                        help='File containing sequence alignment (FASTA or CLUSTAL) format')

    parser.add_argument('outfile',
                        help='Path to the output file')

    parser.add_argument('-cfg',
                        help='Path to file that contains parameters')

    parser.add_argument('-geneconv', '--geneconv',
                        help='Perform GeneConv analysis',
                        action='store_true')

    parser.add_argument('-bootscan', '--bootscan',
                        help='Perform Bootscan analysis',
                        action='store_true')

    parser.add_argument('-maxchi', '--maxchi',
                        help='Perform MaxChi analysis',
                        action='store_true')

    parser.add_argument('-siscan', '--siscan',
                        help='Perform Siscan analysis',
                        action='store_true')

    parser.add_argument('-chimaera', '--chimaera',
                        help='Perform Chimaera analysis',
                        action='store_true')

    parser.add_argument('-threeseq', '--threeseq',
                        help='Perform 3Seq analysis',
                        action='store_true')

    parser.add_argument('-rdp', '--rdp',
                        help='Perform RDP analysis',
                        action='store_true')

    parser.add_argument('-quiet', '--quiet',
                        help='Hide progress messages',
                        action='store_true')

    parser.add_argument('-all', '--all',
                        help='Perform all 7 analyses',
                        action='store_true')

    return parser.parse_args()


def openrdp(infile, outfile, cfg=None, run_geneconv=True, run_three_seq=True, run_rdp=True,
                      run_siscan=True, run_maxchi=True, run_chimaera=True, run_bootscan=True, quiet=False):
    # Check that the OS is valid
    platform = sys.platform
    try:
        platform.startswith("win") or platform.startswith("linux") or sys.platform == "darwin"
    except OSError:
        print("OSError: {} is not supported".format(sys.platform))

    with open(infile) as in_handle:
        if infile.endswith('.fa') or infile.endswith('.fasta'):
            names, aln = read_fasta(in_handle)

    if not valid_alignment(aln) and not valid_chars(aln):
        sys.exit(1)

    scanner = Scanner(names, infile, outfile, cfg, run_geneconv, run_three_seq, run_rdp,
                      run_siscan, run_maxchi, run_chimaera, run_bootscan, quiet)
    scanner.run_scans(aln)


def main():
    args = parse_args()

    # Retrieve arguments
    infile = args.infile
    outfile = args.outfile
    cfg = args.cfg
    if args.all:
        run_geneconv = True
        run_three_seq = True
        run_rdp = True
        run_siscan = True
        run_maxchi = True
        run_chimaera = True
        run_bootscan = True
    else:
        run_geneconv = args.geneconv
        run_three_seq = args.threeseq
        run_rdp = args.rdp
        run_siscan = args.siscan
        run_maxchi = args.maxchi
        run_chimaera = args.chimaera
        run_bootscan = args.bootscan

    openrdp(infile, outfile, cfg, run_geneconv, run_three_seq, run_rdp,
                      run_siscan, run_maxchi, run_chimaera, run_bootscan, args.quiet)


if __name__ == '__main__':
    main()
