import argparse
from collections import defaultdict, deque
from operator import concat
from functools import reduce
import numpy as np
import sys

class vertex:
    def __init__(self, index, value):
        self.index = index
        self.value = value
        self.parent = []
        self.children = []
        self.child_genomes = {}
    def add_parent(self, parent):
        try:
            self.parent.extend(parent)
        except TypeError as e:
            self.parent.append(parent)

    def add_children(self, children):
        try:
            self.children.extend(children)
        except TypeError as e:
            self.children.append(children)

    def add_child_genomes(self, child_genomes):
        for child, genomes in child_genomes.items():
            assert child in self.children
            self.child_genomes[child] = genomes
# END


unlist = lambda n_list: reduce(concat, n_list)


def topological_sort(args):
    FIN = args.input.name
    FOUT = args.output.name
    INPUTS = deque()
    with open(FIN, 'r') as f:
        for line in f:
            if line == '': continue
            line = line.strip().split(',')
            INPUTS.append(line)

    INPUTS = [['^'] + s + ['$'] for s in INPUTS]     # Add start (^) and end ($) nodes
    # Assert that the inputs do not have duplications
    for i, input in enumerate(INPUTS):
        try:
            assert len(set(input)) == len(input)
        except AssertionError as e:
            raise AssertionError(f'Sequence {i} contains duplicated nodes. Exiting.')
            sys.exit()
    # Setup 1: initial graph. This graph is directed and contain cycles

    # All edges : generated as pairs of adjacent nodes in each sequence
    edges = defaultdict(deque)
    for i, input in enumerate(INPUTS):
        for j, c in enumerate(input[:-1]):
            edges[(c, input[j+1])].append(i)
    # Map between nodeids and edges exiting from the respective nodes
    nodeout = defaultdict(deque)
    for k, v in edges.items():
        nodeout[k[0]].append([list(v), k[1]])
    # Map between nodeids and edges entering the respective nodes
    nodeins = defaultdict(deque)
    for k, v in edges.items():
        nodeins[k[1]].append([list(v), k[0]])

    # Step 2: Remove cycles. Create a Directed acyclic graph by duplicating nodes that form the cycle.
    nodelist = deque()
    # Initiate DAG
    DAG = deque()
    index = 0
    new_vertex = lambda index, value: (vertex(index, value), index+1)
    v, index = new_vertex(index, '^')
    DAG.append(v)
    v_cur = v
    children = deque()
    child_genomes = dict()
    for n in nodeout['^']:
        # v, index = new_vertex(len(DAG), n[1])
        v = vertex(len(DAG), n[1])
        DAG.append(v)
        v.add_parent(v_cur.index)
        children.append(v.index)
        child_genomes[v.index] = n[0]
        nodelist.append(v.index)
    v_cur.add_children(children)
    v_cur.add_child_genomes(child_genomes)


    # Get node with the fewest "in/entering" genomes as the next node to add in the graph
    nodeins = defaultdict(set)
    for n in nodelist:
        nodeins[n].update([j for p in DAG[n].parent for j in DAG[p].child_genomes[n]])
    toadd = lambda x: sorted(x, key=lambda n: len(nodeins[n]))

    count = 0
    while len(nodelist) > 0:
        count += 1
        if count == 100000:
            break
        if len(nodelist) == 1 and DAG[nodelist[0]].value == '$':
            break
        sorted_nodes = toadd(nodelist)
        for n in sorted_nodes:
            n_ins = nodeins[n]
            n_outs = nodeout[DAG[n].value]
            n_outs_selected = [i for i in n_outs if set(i[0]).intersection(n_ins) != set()]
            if set([j for i in n_outs_selected for j in i[0]]) == n_ins:
                v_cur = DAG[n]
                children = deque()
                child_genomes = dict()
                for n_out in n_outs_selected:
                    matched_node = [i for i in nodelist if n_out[1] == DAG[i].value]
                    if len(matched_node) > 1:
                        print('ERROR: Multiple matches found')
                    elif len(matched_node) == 1:
                        v = DAG[matched_node[0]]
                    else:
                        v = vertex(len(DAG), n_out[1])
                        DAG.append(v)
                        nodelist.append(v.index)
                    children.append(v.index)
                    child_genomes[v.index] = n_out[0]
                    v.add_parent(v_cur.index)
                    nodeins[v.index].update(n_out[0])
                v_cur.add_children(children)
                v_cur.add_child_genomes(child_genomes)
                nodelist.remove(n)
                _ = nodeins.pop(n)
                break

    # Step 3: Get topological order
    count = 0
    alignment = deque()
    children = set([DAG[0].index])
    terminal = 0
    while len(children) > 0:
        count += 1
        if count == 100:
            break
        if len(children) == 1:
            if DAG[list(children)[0]].value == '$':
                break
        for c in children:
            if len(DAG[c].parent) == 0:
                v_cur = DAG[c]
                al = np.array(['-'] * len(INPUTS))
                gens = list(v_cur.child_genomes.values())
                gens = unlist(gens)
                al[gens] = v_cur.value
                if all(al == '-'):
                    terminal += 1
                    continue
                if v_cur.value not in {'^', '$'}:
                    alignment.append(al)
                children.remove(c)
                children.update(v_cur.children)
                for i in v_cur.children:
                    DAG[i].parent.remove(c)
                break
    with open(FOUT, 'w') as f:
        for l in np.array(alignment).T:
            f.write(','.join(l) + '\n')
    return
# END
