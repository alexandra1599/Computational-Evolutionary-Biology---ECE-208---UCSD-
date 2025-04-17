from treeswift import *
from sys import argv
from random import random, randint
from math import log

tree = read_tree_newick(argv[1])
nquery = int(argv[2])
outfile = argv[3]

def randomize_query(L,p):
    q = []
    for x in L:
        r = random()
        if r < p:
            q.append(x)
    return q
    
with open(outfile,'w') as fout:
    L = [node.label for node in tree.traverse_preorder()]
    n = len(L)
    p = log(n)/n
    for i in range(nquery):
        first_idx = randint(0,n-1)
        q =  [L[first_idx]]
        q += randomize_query(L[:first_idx]+L[first_idx+1:],p)           
        fout.write(' '.join(q) + "\n")
