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


test_path = dirname(realpath(__file__))+"/test_data/checking"
test_ids = [str(x) for x in range(1,6)]
test_cases = ['100','1000','10000','100000']

class Tests_lcas(unittest.TestCase):
    def __sanity_check__(self,ID):
        for case in test_cases:
            tree_file = test_path + "/test_" + ID + "_" + case + ".nwk" 
            query_file = test_path + "/test_" + ID + "_" + case + "_queries.txt"
            ref_file = test_path + "/test_" + ID + "_" + case + "_LCAs.txt"
            try:
                T = read_tree_newick(tree_file)
            except:        
                self.assertTrue(False,msg="Couldn't parse the tree file " + tree_file + "!")

            L = set(node.label for node in T.traverse_preorder())

            try:
                with open(query_file,'r') as f:
                    nq = 0
                    for line in f:
                        q = line.strip().split()
                        nq += 1
                        for x in q:
                            self.assertTrue(x in L, msg="Invalid query file " + query_file + "! Node label " + x + " does not exist in the tree!")
            except:
                self.assertTrue(False,msg="Couldn't parse the query file " + query_file + "!")            
           
            try:
                with open(ref_file,'r') as f:
                    lcas = f.read().strip().split()
                    self.assertTrue(len(lcas)==nq, msg="Invalid reference file " + ref_file + ": answer-query mistmatch!")
            except:
                self.assertTrue(False,msg="Couldn't parse the reference file " + ref_file + "!")            

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
            tree_file = test_path + "/test_" + ID + "_" + case + ".nwk" 
            query_file = test_path + "/test_" + ID + "_" + case + "_queries.txt"
            ref_file = test_path + "/test_" + ID + "_" + case + "_LCAs.txt"

            T = read_tree_newick(tree_file)
            Q = []
            with open(query_file,'r') as f:
                for line in f:
                    Q.append(line.strip().split())
                    
            try:
                if platform.system() == "Windows":
                    LCAs = func_timeout(TIMEOUT,find_LCAs,args=(T,Q,))
                else:
                    signal.alarm(TIMEOUT)
                    LCAs = find_LCAs(T,Q)
                    signal.alarm(0)
            except:
                self.assertTrue(False,msg="Failed test_" + ID + "_" + case + ": couldn't produce output in time!")    
            
            with open(ref_file,'r') as fin:
                LCAs_true = fin.read().strip().split()
            
            for i,(x,y) in enumerate(zip(LCAs_true,LCAs)):                
                self.assertTrue(x==y, msg="Failed test_" + ID + "_" + case + ": Wrong LCA answer for query " + str(i+1) + ". Expect: " + str(x) + ". Your answer: " + str(y))
    
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
