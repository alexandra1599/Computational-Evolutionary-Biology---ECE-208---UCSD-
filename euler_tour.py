from treeswift import *
from treelib.todo import *
from sys import argv


if __name__=="__main__":
    '''
    This is how we handle loading the input tree, running your function, and printing the output
    '''
    if len(argv) != 3:
        print("USAGE: %s <newick_tree> <output_file>" % argv[0]); exit(1)
    
    tree = read_tree_newick(argv[1])
    outfile = argv[2]
    E,F,H = euler_tour(tree)

    with open(outfile,'w') as fout:
        fout.write(' '.join(x for x in E)+"\n")
        fout.write(' '.join(x+":"+str(F[x]) for x in F)+"\n")
        fout.write(' '.join(x+":"+str(H[x]) for x in H))
