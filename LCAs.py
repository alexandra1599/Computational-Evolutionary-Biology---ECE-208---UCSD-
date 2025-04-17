from treeswift import *
from treelib.todo import *
from sys import argv


if __name__=="__main__":
    '''
    This is how we handle loading the input tree, running your function, and printing the output
    '''
    if len(argv) != 4:
        print("USAGE: %s <newick_tree> <query_file> <output_file>" % argv[0]); exit(1)
    
    tree = read_tree_newick(argv[1])
    queryFile = argv[2]
    outFile = argv[3]

    queries = []
    with open(queryFile,'r') as fin:
        for line in fin:
            q = [x.strip() for x in line.strip().split()]
            queries.append(q)

    lcas = find_LCAs(tree,queries)

    with open(outFile,'w') as fout:
        fout.write('\n'.join(lcas))
