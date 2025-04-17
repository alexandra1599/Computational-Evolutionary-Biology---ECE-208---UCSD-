import unittest
from treelib.todo import *
from os.path import dirname, realpath, join, normpath
from os import listdir, remove
import signal
from treeswift import *
import platform

if platform.system() == "Windows":
    from func_timeout import func_timeout, FunctionTimedOut 

TIMEOUT = 30 # seconds
EPS = 1e-4


test_path = normpath(join(dirname(realpath(__file__)),"test_data","checking"))
test_ids = [str(x) for x in range(1,6)]
test_cases = ['100','1000','10000','100000']

def is_valid_euler(tree,E):
    visited_edges = {}
    for node in tree.traverse_preorder():
        if not node.is_root():
            visited_edges[(node.label,node.parent.label)] = False
            visited_edges[(node.parent.label,node.label)] = False
    for i,u in enumerate(E[:-1]):
        v = E[i+1]
        if visited_edges[(u,v)]:
        # this edge has been marked before, so the euler tour is invalid
            return False
        else:
            visited_edges[(u,v)] = True
    # check if any edge is missing in E        
    for e in visited_edges:
        if not visited_edges[e]:
            return False
    return True            

class Tests_euler(unittest.TestCase):
    def __sanity_check__(self,ID):
        for case in test_cases:
            tree_file = normpath(join(test_path, "test_" + ID + "_" + case + ".nwk"))
            ref_file = normpath(join(test_path, "test_" + ID + "_" + case + "_euler.txt"))
            try:
                T = read_tree_newick(tree_file)
            except:        
                self.assertTrue(False,msg="Couldn't parse the tree file " + tree_file + "!")
            
            with open(ref_file,'r') as f:
                E = f.readline().strip().split()
                # Check the validity of the Euler tour
                self.assertTrue(is_valid_euler(T,E), msg="Invalid Euler tour reported by reference file " + ref_file + "!")
                  

    def test_01_sanity(self):
        ID = '1'
        self.__sanity_check__(ID)
    
    def test_02_sanity(self):
        ID = '2'
        self.__sanity_check__(ID)
    
    def test_03_sanity(self):
        ID = '3'
        self.__sanity_check__(ID)
    
    def test_04_sanity(self):
        ID = '4'
        self.__sanity_check__(ID)
    
    def test_05_sanity(self):
        ID = '5'
        self.__sanity_check__(ID)

    def __run_test__(self,case):
        for ID in test_ids:
            tree_file = normpath(join(test_path, "test_" + ID + "_" + case + ".nwk"))
            ref_file = normpath(join(test_path, "test_" + ID + "_" + case + "_euler.txt"))

            T = read_tree_newick(tree_file)

            if platform.system() == "Windows":
                 E,F,H = func_timeout(TIMEOUT,euler_tour,args=(T,))
            else:
                signal.alarm(TIMEOUT)
                E,F,H = euler_tour(T)
                signal.alarm(0)

            with open(ref_file,'r') as fin:
                fin.readline() # skip E, because there are many possible answers so we cannot directly compare the Euler tours
                fin.readline() # skip F, because there are many possible answers so we cannot directly compare the Euler tours
                #F_true = {x:int(y) for (x,y) in [s.split(":") for s in fin.readline().strip().split()] }
                H_true = {x:int(y) for (x,y) in [s.split(":") for s in fin.readline().strip().split()] }
            
                self.assertTrue(is_valid_euler(T,E), msg="Failed test_" + ID + "_" + case + ": Invalid Euler tour!")
            F_true = {}
            for i,x in enumerate(E):
                if not x in F_true:
                    F_true[x] = i
            is_valid_F = True
            is_valid_H = True
            for x in F_true:
                self.assertTrue(F[x]==F_true[x], msg="Failed test_" + ID + "_" + case + ": Wrong F!")
                self.assertTrue(H[x]==H_true[x], msg="Failed test_" + ID + "_" + case + ": Wrong H!")                
    
    def test_06_correctness(self):
        case = '100'
        self.__run_test__(case)
    def test_07_correctness(self):
        case = '1000'
        self.__run_test__(case)
    def test_08_correctness(self):
        case = '10000'
        self.__run_test__(case)
    def test_09_correctness(self):
        case = '100000'
        self.__run_test__(case)
