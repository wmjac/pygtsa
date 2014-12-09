# Copyright (C) 2014 William M. Jacobs

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

import argparse
import networkx as nx
import pygtsa.cgraph as cgraph

class Assembly(nx.Graph):
    def __init__(self, list_of_edges):
        super().__init__(list_of_edges)

    @staticmethod
    def sorted_edge(edge):
        return tuple(sorted(edge))

    def adjacent_edges(self, fragment):
        '''Return the set of edges in self that are outside but adjacent to fragment'''
        neighbors = set([])
        for node in fragment.nodes():
            neighbors |= set(Assembly.sorted_edge((node, neighbor)) for neighbor in self[node].keys())
        return neighbors.difference(fragment.edges())

    def leaves(self):
        '''Return the set of leaves in self'''
        return set(Assembly.sorted_edge((node, self.neighbors(node)[0])) for node in self.nodes() if len(self[node]) == 1)

    def bridges(self):
        '''Return the set of bridges in self'''
        return cgraph.Graph(self).bridges()

    def disconnect_edge(self, edge):
        '''Remove edge and attached node(s) if a leaf'''
        degrees = tuple(len(self[edge[i]]) for i in range(2))
        self.remove_edge(*edge)
        for i in range(2):
            if degrees[i] == 1:
                self.remove_node(edge[i])

    def connect_edge(self, edge):
        '''Add edge'''
        self.add_edge(*edge)

    def print_info(self):
        print("number of nodes:", self.number_of_nodes())
        print("number of edges:", self.number_of_edges())

    def write_edges(self, stream):
        for edge in self.edges():
            stream.write("%d %d\n" % edge)

    @staticmethod
    def read_edges(stream):
        G = nx.Graph()
        for line in stream:
            if line[0] == '#':
                continue
            t = line.split()
            if len(t) == 0:
                continue
            try:
                G.add_edge(t[0], t[1])
            except (KeyError, TypeError):
                raise TypeError("ERROR: cannot interpret line in input structure file:%d" % line)
        return Assembly(nx.convert_node_labels_to_integers(G).edges())

    def write_dot(self, stream, circsize=0.6, edgelength=1., fontface='Arial', fillcolor='#cccccc'):
        stream.write("graph G {\n")
        for node in self.nodes():
            stream.write("  %d [shape=circle,height=%g,width=%g,fontname=\"%s\",style=filled,fillcolor=\"%s\"];\n" 
                    % (node + 1, circsize, circsize, fontface, fillcolor))
        for edge in self.edges():
            stream.write("  %d -- %d [len=%g];\n" % (edge[0] + 1, edge[1] + 1, edgelength))
        stream.write("}\n")

    @staticmethod
    def read_dot(self, path):
        return Assembly(nx.convert_node_labels_to_integers(nx.Graph(nx.read_dot(path))))

    @staticmethod
    def read(path):
        try:
            with open(path, 'r') as f:
                return Assembly.read_edges(f)
        except TypeError:
            pass
        try:
            return Assembly.read_dot(path)
        except TypeError:
            pass
        raise TypeError("could not interpret assembly file %s" % path)

    def enumerate_subgraphs(self, Emax=None, verbose=False):
        '''Find all connected subgraphs (with E <= Emax) of the assembly by direct enumeration'''
        subgraphs = {1 : set((edge,) for edge in self.edges())}
        if Emax == None:
            Emax = self.number_of_edges()
        for i in range(2, Emax + 1):
            subgraphs[i] = set([])
            for h in subgraphs[i-1]:
                nodes = set([edge[0] for edge in h]) | set([edge[1] for edge in h])
                for node in nodes:
                    for edge in self.edges_iter(node):
                        edge = tuple(sorted(edge))
                        if edge not in h:
                            new_edge_set = tuple(sorted(h + (edge,)))
                            subgraphs[i].add(new_edge_set)
            if verbose:
                print("subgraphs with %d edges: (n = %d)" % (i, len(subgraphs[i])))
        return subgraphs

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('structure', type=str, help="path to input structure file")
    parser.add_argument('--output-dot', metavar='PATH', type=str, default=False, help="output structure in the DOT format to PATH [False]")
    parser.add_argument('--output-edges', metavar='PATH', type=str, default=False, help="output structure edges to PATH [False]")
    clargs = parser.parse_args()

    target = Assembly.read(clargs.structure)

    target.print_info()

    if clargs.output_dot:
        print("Writing DOT file to %s" % clargs.output_dot)
        with open(clargs.output_dot, 'w') as f:
            target.write_dot(f)
    if clargs.output_edges:
        print("Writing edges to %s" % clargs.output_edges)
        with open(clargs.output_edges, 'w') as f:
            target.write_edges(f)
