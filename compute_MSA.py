#!/usr/bin/env python3

if __name__ == "__main__":
    '''
    This is how we handle loading the input dataset, running your function, and printing the output
    '''
    import argparse
    from treeswift import *
    from msa.todo import MSA
    from msa.utils import *

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input_sequences', required=False, type=str, default='stdin', help="Input Sequences (FASTA format)")
    parser.add_argument('-r', '--indel_rate', required=True, type=float, help="Indel Rate")
    parser.add_argument('-t', '--guidance_tree', required=True, type=str, help="Guidance tree (Newick format)")
    parser.add_argument('-o', '--output_alignment', required=False, type=str, default='stdout', help="Output Alignment (FASTA format)")
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
    tree = read_tree_newick(args.guidance_tree)
    indel_rate = args.indel_rate

    aln = MSA(seqs,tree,indel_rate)
    for ID in aln:
        outfile.write(">%s\n%s\n" % (ID,aln[ID]))
    outfile.close()
