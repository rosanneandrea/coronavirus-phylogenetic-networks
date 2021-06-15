import networkx as nx
import ast
import re
from copy import deepcopy
import sys
import csv

########  Convert Newick to a networkx Digraph with labels (and branch lengths)
########  Author: Remie Janssen
########  Edited by Rosanne Wallin


#Write length newick: convert ":" to "," and then evaluate as list of lists using ast.literal_eval
# Then, in each list, the node is followed by the length of the incoming arc.
# This only works as long as each branch has length and all internal nodes are labeled.
def Newick_To_Tree(newick, current_labels = dict()):
    newick=newick[:-1]
    distances = False
    if ":" in newick:
        distances = True
        if not "'" in newick and not '"' in newick:
            newick = re.sub(r"([,\(])([a-zA-Z\d]+)", r"\1'\2", newick)
            newick = re.sub(r"([a-zA-Z\d]):", r"\1':", newick)
        newick = newick.replace(":",",")
    else:
        if not "'" in newick and not '"' in newick:
            newick = re.sub(r"([,\(])([a-zA-Z\d]+)", r"\1'\2", newick)
            newick = re.sub(r"([a-zA-Z\d])([,\(\)])", r"\1'\2", newick)
    newick = newick.replace("(","[")
    newick = newick.replace(")","]")
    nestedtree = ast.literal_eval(newick)
    edges, leaves, current_labels, current_node = NestedList_To_Tree(nestedtree,1001,current_labels, distances = distances)
    tree = nx.DiGraph()        
    if distances:
        edges.append((1000,1001,0))
        tree.add_weighted_edges_from(edges,weight='length')        
    else:
        edges.append((1000,1001))
        tree.add_edges_from(edges)
    return tree, leaves, current_labels, distances


#Auxiliary function to convert list of lists to graph 
def NestedList_To_Tree(nestedList, next_node, current_labels, distances = False):
    edges = []
    leaves = set()
    top_node = next_node
    current_node = next_node+1
    if distances:
        for i in range(0,len(nestedList),2):
            t = nestedList[i]
            length = nestedList[i+1]
            if type(t)==list: #Not a leaf 
                edges.append((top_node,current_node,length))
                extra_edges, extra_leaves, current_labels, current_node = NestedList_To_Tree(t,current_node,current_labels,distances=distances)
            else: #A leaf
                if str(t) not in current_labels:
                    current_labels[str(t)] = len(current_labels) 
                edges.append((top_node,current_labels[str(t)],length))
                extra_edges = []
                extra_leaves = set([current_labels[str(t)]])
            edges = edges + extra_edges
            leaves = leaves.union(extra_leaves)
    else:
        for t in nestedList:
            if type(t)==list: 
                edges.append((top_node,current_node))
                extra_edges, extra_leaves, current_labels, current_node = NestedList_To_Tree(t,current_node,current_labels)
            else: 
                if str(t) not in current_labels:
                    current_labels[str(t)] = len(current_labels) 
                edges.append((top_node,current_labels[str(t)]))
                extra_edges = []
                extra_leaves = set([current_labels[str(t)]])
            edges = edges + extra_edges
            leaves = leaves.union(extra_leaves)
    return edges, leaves, current_labels, current_node
    
    
    
    
    
    
    
#######    
####### Methods for phylogenetic trees

class PhT:
    def __init__(self,newick = None):
        #the actual graph
        self.nw = nx.DiGraph()
        #the set of leaf labels of the network
        self.leaves = set()
        self.labels = dict()                    #label->node
        self.labels_inv = dict()                #node->label
        if newick:
            self.Tree_From_Newick(newick = newick)
            for x in self.labels:
                self.labels_inv[self.labels[x]]=x
    
    def Tree_From_Newick(self, newick = None, current_labels = dict()):
        self.nw, self.leaves, self.labels, distances = Newick_To_Tree(newick, current_labels)
        return self.labels, distances
        
    def Is_Cherry(self,x,y):
        if (not x in self.leaves) or (not y in self.leaves):
            return False
        px=-1
        py=-1
        for p in self.nw.predecessors(x):
            px=p
        for p in self.nw.predecessors(y):
            py=p
        return px==py
   


    def Clean_Node(self,v):
        if self.nw.out_degree(v)==1 and self.nw.in_degree(v)==1:
            pv=-1
            for p in self.nw.predecessors(v):
                pv=p
            cv=-1
            for c in self.nw.successors(v):
                cv=c
            self.nw.add_edges_from([(pv,cv,self.nw[pv][v])])
            if 'length' in self.nw[pv][v] and 'length' in self.nw[v][cv]:
                self.nw[pv][cv]['length']=self.nw[pv][v]['length']+self.nw[v][cv]['length']
            self.nw.remove_node(v)
            return True
        return False
 

    def reduce_pair(self,x,y):
        if not x in self.leaves or not y in self.leaves:
            return False
        px=-1
        py=-1
        for p in self.nw.predecessors(x):
            px=p
        for p in self.nw.predecessors(y):
            py=p 
        if self.Is_Cherry(x,y):
            self.nw.remove_node(x)
            self.leaves.remove(x)
            self.Clean_Node(py)
            return True
        return False

    
    
    #Reduce the tree with the temporal cps and return the pairs that are reduced
    #Assumes the tree is binary
    def Reduce_With_Temporal_CPS(self,temp_cps):
        reducing_tree = deepcopy(self)
        reduced = []
        temp_cps = [self.labels[x] for x in temp_cps]
        for x in temp_cps:
            px=-1
            for p in reducing_tree.nw.predecessors(x):
                px=p
            ys = []
            for cpx in reducing_tree.nw.successors(px):
                if cpx in reducing_tree.leaves and cpx!=x:
                    ys+=[cpx]
                    #now we assume the tree is binary, so that x has at most one sibling
                    break
            if ys:
                reducing_tree.reduce_pair(x,ys[0])
            reduced+=[[self.labels_inv[y] for y in ys]]
        return reduced
    
    
    
    



########  A class for phylogenetic networks
########

class PhN:
    def __init__(self, seq = None):
        #the actual graph
        self.nw = nx.DiGraph()
        #the set of leaf labels of the network
        self.leaves = set()
        #a dictionary giving the node for a given leaf label
        self.labels = dict()
        #the number of nodes in the graph
        self.no_nodes = 0
        #a dictionary giving the leaf label for each pendant node
        self.leaf_nodes = dict()
        if seq:
            #Creates a phylogenetic network from a cherry picking sequence:
            for pair in reversed(seq):
                self.add_pair(*pair)

            
               
    #A method for adding a pair, using the construction from a sequence
    def add_pair(self,x,y):
        if len(self.leaves)==0:
            self.nw.add_edges_from([(0,1),(1,2),(1,3)])
            self.leaves = set([x,y])
            self.labels[x]=2
            self.labels[y]=3
            self.leaf_nodes[2]=x
            self.leaf_nodes[3]=y
            self.no_nodes=4
            return True
        if y not in self.leaves:
            return False
        node_y=self.labels[y]
        if x not in self.leaves:
            self.nw.add_edges_from([(node_y,self.no_nodes),(node_y,self.no_nodes+1)])
            self.leaves.add(x)
            self.leaf_nodes.pop(self.labels[y],False)
            self.labels[y]=self.no_nodes
            self.labels[x]=self.no_nodes+1
            self.leaf_nodes[self.no_nodes]=y
            self.leaf_nodes[self.no_nodes+1]=x
            self.no_nodes+=2          
        else:
            node_x=self.labels[x]
            for parent in self.nw.predecessors(node_x):
                px = parent
            if self.nw.in_degree(px)>1:
                self.nw.add_edges_from([(node_y,px),(node_y,self.no_nodes)])
                self.leaf_nodes.pop(self.labels[y],False)
                self.labels[y]=self.no_nodes
                self.leaf_nodes[self.no_nodes]=y
                self.no_nodes+=1
            else: 
                self.nw.add_edges_from([(node_y,node_x),(node_y,self.no_nodes),(node_x,self.no_nodes+1)])
                self.leaf_nodes.pop(self.labels[x],False)
                self.leaf_nodes.pop(self.labels[y],False)
                self.labels[y]=self.no_nodes
                self.labels[x]=self.no_nodes+1
                self.leaf_nodes[self.no_nodes]=y
                self.leaf_nodes[self.no_nodes+1]=x
                self.no_nodes+=2
        return True        


    def Print_Edges(self):
        all_edges = ""
        for e in self.nw.edges:
            node0 = e[0]
            node1 = e[1]
            #put back the original labels of leaves 
            #    (the nodes are integers with the original names stored in the leaf_nodes dictionary)
            if node0 in self.leaf_nodes:
                node0 = self.leaf_nodes[node0]
            if node1 in self.leaf_nodes:
                node1 = self.leaf_nodes[node1]
            all_edges += str(node0)+" "+str(node1)+"\r\n"
        return all_edges
    



################################################################################
################################################################################
################################################################################
########                                                           #############
########                         CutTree CLASS                     #############
########                                                           #############
################################################################################
################################################################################
################################################################################


#A class that represents a network as a tree where hybrid edges have been cut at the hybrid nodes.
#Mainly used as an intermediate to find the Newick string of a network.

class CutTree:
    def __init__(self, network = None, current_node = None, leaf_labels= dict()):
         self.hybrid_nodes = dict()
         self.no_of_hybrids = 0
         self.root = None
         self.nw = deepcopy(network)
         self.current_node = current_node
         self.leaf_labels = leaf_labels
         if not self.current_node:
             self.current_node = 2*len(self.nw)
         if network:
             self.Find_Root()
             network_nodes = list(self.nw.nodes)
             for node in network_nodes:
                 if self.nw.in_degree(node)>1:
                     self.no_of_hybrids+=1
                     enumerated_parents = list(enumerate(self.nw.predecessors(node))) 
                     for i,parent in enumerated_parents:
                         if i==0:
                             self.hybrid_nodes[node]=self.no_of_hybrids
                         else:
                             self.nw.add_edges_from([(parent,self.current_node,self.nw[parent][node])])
                             self.nw.remove_edge(parent,node)
                             self.hybrid_nodes[self.current_node] = self.no_of_hybrids
                             self.current_node+=1
#             self.CheckLabelSet()

    def Find_Root(self):
        for node in self.nw.nodes:
            if self.nw.in_degree(node)==0:
                self.root = node  
                return node           

    def Newick(self,probabilities = False):
        return self.Newick_Recursive(self.root,probabilities = probabilities)+";"

    def Newick_Recursive(self,root,probabilities = False):
        if self.nw.out_degree(root)==0:
            if root in self.hybrid_nodes:
                return "#H"+str(self.hybrid_nodes[root])
            elif root in self.leaf_labels:
                return self.leaf_labels[root]
            return str(root)
        Newick = ""
        for v in self.nw.successors(root):
            Newick+= self.Newick_Recursive(v,probabilities)
            Newick+= ","
        Newick = "("+Newick[:-1]+")"
        if root in self.hybrid_nodes:
            Newick += "#H"+str(self.hybrid_nodes[root])
        return Newick




    
    
##################################################
##################################################
#####                                        #####
#####              Main stuff                #####
#####                                        #####
##################################################
##################################################

############    
###### Input


#Read the arguments

edges = False
help  = False
i = 2
out_file = None
while i < len(sys.argv):
    arg= sys.argv[i]
    if arg == "-e":
        edges = True   
    if arg == "-h":
        help = True   
    if arg == "-o":
        i+=1
        out_file = sys.argv[i]
    i += 1

#Output the help text to the terminal if no argument is given, or if the help option is chosen.
if len(sys.argv)<2 or help:
    print("call the program with `python newick_file_path sequence_file_path [-o output_filename] [-e -h] '\n the optional argument -e makes the output edges instead of newick.\n the optional argument -o makes the program write the output to the given file instead of to the terminal")
    sys.exit()

sequence_file = sys.argv[1]

#Read each line of the input file with name set by "seq_file"
f = open("./"+sequence_file, "rt")
reader = csv.reader(f, delimiter='~', quotechar='|')
for row in reader:
    seq = row[0].replace("(","").replace(")","").split(" ->  ")
f.close()
tc_cps = []
for couple in seq:
    couple = couple.split(",")
    couple = (couple[0], couple[1])
    tc_cps.append(couple)
    

"""
##################
####### Conversion

# convert newick to trees
trees = [PhT(newick=n) for n in newicks]

#Find all second elements for the pairs that are reduced by the temporal cps
second_elements_of_pairs = [set() for x in temp_cps]
for t in trees:
    second_elements_t = t.Reduce_With_Temporal_CPS(temp_cps)
    for i,ys in enumerate(second_elements_t):
        second_elements_of_pairs[i]|=set(ys)


#Convert the sequence and the second elements into a regular (tree-child) cps
tc_cps = []
for x,ys in zip(temp_cps,second_elements_of_pairs):
    for y in ys:
        tc_cps+=[(x,y)]
"""


########################################################################################################################################
#To optimally convert a temporal CPS to a TC TCS for non-binary trees, we may need to solve some kind of hitting set problem!
#   indeed, to remove x from a multifurcating cherry (x,y_1,...y_n), we need to reduce only one of (x,y_1),...,(x,y_n)
#   this holds for each tree in the input, so we need to hit some set {y_1,...,y_n} for each tree.
#   hence, to optimally convert the element x of the temp-cps to a set of tc-cps pairs, we need to solve this hitting set problem
#   because we work with binary trees, we may ignore this, as each set we need to hit has size 1.
########################################################################################################################################


#Convert the cps into a network and print the edges
network = PhN(seq=tc_cps)




#############
###### Output


if edges:
    #Print edges
    output = network.Print_Edges()
else:
    #Print newick
    cuttree = CutTree(network = network.nw, current_node = network.no_nodes, leaf_labels = network.leaf_nodes)
    output = cuttree.Newick()



if out_file:
    f= open(out_file,"w+")
    f.write(output)
    f.close()
else:
    print(output)


   


