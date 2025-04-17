#!/usr/bin/env python3
AMINOS   = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
BLOSUM62 = [[ 4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2],
            [ 0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2],
            [-2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3],
            [-1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2],
            [-2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3],
            [ 0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3],
            [-2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2],
            [-1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1],
            [-1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2],
            [-1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1],
            [-1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1],
            [-2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2],
            [-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3],
            [-1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1],
            [-1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2],
            [ 1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2],
            [ 0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2],
            [ 0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1],
            [-3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2],
            [-2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2, 7]]

# None = STOP
# 0 = (i-1, j,   k)
# 1 = (i,   j-1, k)
# 2 = (i,   j,   k-1)
# 3 = (i-1, j-1, k)
# 4 = (i-1, j,   k-1)
# 5 = (i,   j-1, k-1)
# 6 = (i-1, j-1, k-1)
def threeway_align(s1, s2, s3, indel_rate, VERBOSE=False):
    '''
    This function computes the threeway sequence alignment between strings s1, s2, and s3 using a relative indel rate of indel_rate
    :param s1: string s1
    :param s2: string s2
    :param s3: string s3
    :param indel_rate: the relative indel rate
    '''
    B = {AMINOS[i]:{AMINOS[j]:BLOSUM62[i][j] for j in range(len(AMINOS))} for i in range(len(AMINOS))}
    from math import log2
    gap = log2(0.92*indel_rate)
    two_gap = 2.*gap # I use 2*gap penalty everywhere, so avoid recomputing it each time

    # initialize (S[i][j][k] = (score,arrow))
    if VERBOSE:
        from sys import stderr; print("Initializing cube", file=stderr)
    S = [[[None for k in range(len(s3)+1)] for j in range(len(s2)+1)] for i in range(len(s1)+1)]
    S[0][0][0] = (0, None)        # corner

    # fill in cube axes (need 2*gap penalty because the character that appears has gap with each of the other 2 strings)
    if VERBOSE:
        print("Filling cube axes", file=stderr)
    for i in range(1, len(s1)+1): # s1 axis
        S[i][0][0] = (two_gap*i, 0)
    for j in range(1, len(s2)+1): # s2 axis
        S[0][j][0] = (two_gap*j, 1)
    for k in range(1, len(s3)+1): # s3 axis
        S[0][0][k] = (two_gap*k, 2)

    # fill in cube faces
    if VERBOSE:
        print("Filling cube faces", file=stderr)
    for i in range(1, len(s1)+1): # (s1,s2) face
        for j in range(1, len(s2)+1):
            S[i][j][0] = max([(S[i-1][j][0][0]+two_gap, 0), (S[i][j-1][0][0]+two_gap, 1), (S[i-1][j-1][0][0]+B[s1[i-1]][s2[j-1]]+two_gap, 3)])
    for i in range(1, len(s1)+1): # (s1,s3) face
        for k in range(1, len(s3)+1):
            S[i][0][k] = max([(S[i-1][0][k][0]+two_gap, 0), (S[i][0][k-1][0]+two_gap, 2), (S[i-1][0][k-1][0]+B[s1[i-1]][s3[k-1]]+two_gap, 4)])
    for j in range(1, len(s2)+1): # (s2,s3) face
        for k in range(1, len(s3)+1):
            S[0][j][k] = max([(S[0][j-1][k][0]+two_gap, 1), (S[0][j][k-1][0]+two_gap, 2), (S[0][j-1][k-1][0]+B[s2[j-1]][s3[k-1]]+two_gap, 5)])

    # fill in rest of cube
    if VERBOSE:
        print("Filling rest of cube", file=stderr)
    for i in range(1, len(s1)+1):
        for j in range(1, len(s2)+1):
            for k in range(1, len(s3)+1):
                S[i][j][k] = max([
                                     (S[i-1][j][k][0]+two_gap, 0),
                                     (S[i][j-1][k][0]+two_gap, 1),
                                     (S[i][j][k-1][0]+two_gap, 2),
                                     (S[i-1][j-1][k][0]+B[s1[i-1]][s2[j-1]]+two_gap, 3),
                                     (S[i-1][j][k-1][0]+B[s1[i-1]][s3[k-1]]+two_gap, 4),
                                     (S[i][j-1][k-1][0]+B[s2[j-1]][s3[k-1]]+two_gap, 5),
                                     (S[i-1][j-1][k-1][0]+B[s1[i-1]][s2[j-1]]+B[s1[i-1]][s3[k-1]]+B[s2[j-1]][s3[k-1]], 6)
                                 ])
    # backtrack to get alignments
    if VERBOSE:
        print("Backtracking to build alignment", file=stderr)
    aln_s1 = ""; aln_s2 = ""; aln_s3 = ""
    i = len(s1); j = len(s2); k = len(s3); arrow = S[i][j][k][1]
    while arrow is not None:
        if arrow == 0:
            aln_s1 += s1[i-1]; aln_s2 += '-'; aln_s3 += '-'; i -= 1
        elif arrow == 1:
            aln_s1 += '-'; aln_s2 += s2[j-1]; aln_s3 += '-'; j -= 1
        elif arrow == 2:
            aln_s1 += '-'; aln_s2 += '-'; aln_s3 += s3[k-1]; k -= 1
        elif arrow == 3:
            aln_s1 += s1[i-1]; aln_s2 += s2[j-1]; aln_s3 += '-'; i -= 1; j -= 1
        elif arrow == 4:
            aln_s1 += s1[i-1]; aln_s2 += '-'; aln_s3 += s3[k-1]; i -= 1; k -= 1
        elif arrow == 5:
            aln_s1 += '-'; aln_s2 += s2[j-1]; aln_s3 += s3[k-1]; j -= 1; k -= 1
        elif arrow == 6:
            aln_s1 += s1[i-1]; aln_s2 += s2[j-1]; aln_s3 += s3[k-1]; i -= 1; j -= 1; k -= 1
        else:
            raise ValueError("Invalid arrow: %s" % str(arrow))
        arrow = S[i][j][k][1]
    return aln_s1[::-1],aln_s2[::-1],aln_s3[::-1] # I built them backwards, so reverse them when I return

def read_FASTA(stream):
    '''
    This function reads a FASTA file from a given stream and returns a dictionary mapping identifiers to sequences
    '''
    seqs = {}; name = None; seq = ''
    for line in stream:
        l = line.strip()
        if len(l) == 0:
            continue
        if l[0] == '>':
            if name is not None:
                assert len(seq) != 0, "Malformed FASTA"
                seqs[name] = seq
            name = l[1:]
            assert name not in seqs, "Duplicate sequence ID: %s" % name
            seq = ''
        else:
            seq += l
    assert name is not None and len(seq) != 0, "Malformed FASTA"
    seqs[name] = seq
    return seqs

if __name__ == "__main__":
    '''
    This is how we handle loading the input dataset, running your function, and printing the output
    '''
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input_sequences', required=False, type=str, default='stdin', help="Input Sequences (FASTA format)")
    parser.add_argument('-r', '--indel_rate', required=True, type=float, help="Indel Rate")
    parser.add_argument('-o', '--output_alignment', required=False, type=str, default='stdout', help="Output Alignment (FASTA format)")
    parser.add_argument('-v', '--verbose', action='store_true', help="Verbose")
    args = parser.parse_args()
    assert args.indel_rate >= 0, "Relative indel rate must be non-negative"
    if args.input_sequences == 'stdin':
        from sys import stdin as infile
    else:
        infile = open(args.input_sequences)
    if args.output_alignment == 'stdout':
        from sys import stdout as outfile
    else:
        outfile = open(args.output_alignment,'w')
    seqs = read_FASTA(infile); infile.close()
    assert len(seqs) == 3, "Input FASTA file must have exactly 3 sequences"
    k1,k2,k3 = sorted(seqs.keys())
    a1,a2,a3 = threeway_align(seqs[k1], seqs[k2], seqs[k3], args.indel_rate, VERBOSE=args.verbose)
    assert len(a1) == len(a2) == len(a3), "Aligned sequences must have equal length"
    assert len(a1) != 0, "Returned empty strings"
    for seq,orig in [(a1,seqs[k1]), (a2,seqs[k2]), (a3,seqs[k3])]:
        assert seq.replace('-','') == orig, "Removing gaps from the aligned sequences must yield the original sequences"
    for ID,seq in [(k1,a1), (k2,a2), (k3,a3)]:
        outfile.write(">%s\n%s\n" % (ID,seq))
    outfile.close()
