#!/usr/bin/env python3
from treelib.todo import *
from sys import argv
from treeswift import read_tree_newick


if __name__ == "__main__":
    '''
    This is how we handle loading the input tree, running your function, and printing the output
    '''
    
    if len(argv) != 3:
        print("USAGE: %s <newick_tree> <output_file>" % argv[0]); exit(1)
    
    tree = read_tree_newick(argv[1])
    outfile = argv[2]
    d_path,d_length = compute_diameter(tree)
    
    with open(outfile,'w') as fout:
        for x in d_path:
            fout.write(x + " ")
        fout.write("\n" + str(d_length))    
