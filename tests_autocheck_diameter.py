import unittest
from treelib.todo import compute_diameter
from os.path import dirname, realpath, join, normpath
from os import listdir, remove
import signal
from treeswift import *
import platform

if platform.system() == "Windows":
    from func_timeout import func_timeout, FunctionTimedOut 

TIMEOUT = 60 # seconds
EPS = 1e-4


test_path = normpath(join(dirname(realpath(__file__)),"test_data","checking"))
test_ids = [str(1),str(2), str(3), str(4) ,str(5)]
test_cases = ['100','1000','10000','100000']

def compute_path_length(tree,path):
    if not path:
        return None

    curr_label = path[0]
    curr_node = None
    for node in tree.traverse_postorder():
        if node.label == curr_label:
            curr_node = node
            break
    if curr_node is None:
        return None       
    
    d = 0     
    visited = {curr_label:curr_node}       
    for curr_label in path[1:]:
        if curr_label in visited:
            return None
        if not curr_node.is_root() and curr_node.get_parent().get_label() == curr_label:
            d += curr_node.get_edge_length() 
            curr_node = curr_node.get_parent()
        else:
            found = False
            for c in curr_node.child_nodes():
                if c.get_label() == curr_label:
                    d += c.get_edge_length()
                    curr_node = c
                    found = True
                    break
            if not found:
                return None                    
        visited[curr_label] = curr_node           
    return d         

class Tests_diameter(unittest.TestCase):
    def __sanity_check__(self,ID):
        for case in test_cases:
            tree_file = normpath(join(test_path,"test_" + ID + "_" + case + ".nwk"))
            ref_file = normpath(join(test_path,"test_" + ID + "_" + case + "_diameter.txt"))
            try:
                T = read_tree_newick(tree_file)
            except:        
                self.assertTrue(False,msg="Couldn't parse the tree file " + tree_file + "!")
            try:
                with open(ref_file,'r') as f:
                    path = f.readline().strip().split()
                    d_ref = float(f.readline().strip())
                    # Check the validity of the diameter path (e.g. it must be an actual path of the tree and its length must be d.
                    d_path = compute_path_length(T,path)
                    self.assertTrue(d_path is not None and abs(d_path-d_ref) < EPS, msg="Invalid diameter path reported by reference file " + ref_file + "!")
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
            tree_file = normpath(join(test_path, "test_" + ID + "_" + case + ".nwk"))
            ref_file = normpath(join(test_path,"test_" + ID + "_" + case + "_diameter.txt"))

            T = read_tree_newick(tree_file)

            try:
                if platform.system() == "Windows":
                    path, d = func_timeout(TIMEOUT,compute_diameter,args=(T,))
                else:
                    signal.alarm(TIMEOUT)
                    path, d = compute_diameter(T)
                    signal.alarm(0)
            except:
                self.assertTrue(False,msg="Failed test_" + ID + "_" + case + ": couldn't produce output in time!")    
            
            with open(ref_file,'r') as fin:
                fin.readline()
                d_true = float(fin.readline().strip())
            
            d_path = compute_path_length(T,path)
            
            self.assertTrue(abs(d-d_true)<EPS, msg="Failed test_" + ID + "_" + case + ": Wrong diameter length!")
            self.assertTrue(d_path is not None, msg="Failed test " + ID + "_" + case + ": Invalid diameter path!")      
            self.assertTrue(abs(d_path-d_true)<EPS, msg="Failed test " + ID + "_" + case + ": Wrong diameter path!")                    
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
