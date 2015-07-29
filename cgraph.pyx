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

import math
import numpy as np
from pygtsa.histogram import EVHistogram

cimport cgraph

# Module initialization
cgraph.random_init()
# ---------------------

cdef class Graph:
    def __init__(self, assembly):
        self._graph = cgraph.graph_alloc(assembly.number_of_nodes(), \
                                         max(assembly.degree(v) for v in assembly.nodes()))
        if self._graph == NULL: raise MemoryError()
        for edge in assembly.edges():
            cgraph.graph_add_edge(self._graph, edge[0], edge[1])

    def __dealloc__(self):
        cgraph.graph_free(self._graph)

    def edges_iter(self):
        cdef size_t i, j
        for i in range(self._graph.max_nvertices):
            for j in range(self._graph.nvedges[i]):
                if self._graph.edges[i][j] > i:
                    yield (i, self._graph.edges[i][j])

    def bridges(self):
        graph_find_bridges(self._graph)
        return set((self._graph.bridges[i].v1, self._graph.bridges[i].v2) for i in range(self._graph.nbridges))

cdef histogramEV_c_to_py(assembly, histogramEV_t * h):
    hEV = EVHistogram(assembly)
    cdef int i, j
    for i in range(hEV.h.shape[0]):
        for j in range(hEV.h.shape[1]):
            E = i
            V = i - j + 1
            hEV.h[hEV._index(E, V)] = histogramEV_get(h, E, V)
    return hEV

cdef histogramEV_py_to_c(hEV, histogramEV_t * h):
    cdef int i, j
    for i in range(hEV.h.shape[0]):
        for j in range(hEV.h.shape[1]):
            E = i
            V = i - j + 1
            histogramEV_set(h, E, V, hEV.h[hEV._index(E, V)])

def wl_simulate_EV(target, finit=1., fmin=1.e-4, flatness=0.9, ncheck=100000, output='lnhEV.dat', verbose=True):
    cdef cgraph.Graph target_graph = cgraph.Graph(target)
    cdef cgraph.Graph fragment_graph = cgraph.Graph(target)

    cdef cgraph.histogramEV_t * lnhEV = cgraph.histogramEV_alloc(target_graph._graph)
    if lnhEV == NULL: raise MemoryError()
    histogramEV_fill(lnhEV, 0.);

    cdef double f = finit
    while f >= 0.5 * fmin:
        if verbose:
            print("W-L calculation: f = %g" % f);
        cgraph.wl_sample_subgraphs_EV(target_graph._graph, fragment_graph._graph, lnhEV, f, flatness, ncheck)
        py_lnhEV = histogramEV_c_to_py(target, lnhEV)
        with open(output, 'w') as out:
            out.write("# f = %g\n\n" % f)
            py_lnhEV.write(out)
        f *= 0.5

    cgraph.histogramEV_free(lnhEV)
    return py_lnhEV

def mc_sample_EV(target, fragment, py_lnhEV, nsteps):
    cdef cgraph.Graph target_graph = cgraph.Graph(target)
    cdef cgraph.Graph fragment_graph = cgraph.Graph(fragment)

    cdef cgraph.histogramEV_t * lnhEV = cgraph.histogramEV_alloc(target_graph._graph)
    if lnhEV == NULL: raise MemoryError()
    histogramEV_py_to_c(py_lnhEV, lnhEV)
    cgraph.mc_sample_subgraphs_EV(target_graph._graph, fragment_graph._graph, lnhEV, nsteps)
    cgraph.histogramEV_free(lnhEV)

def mc_simulate_dihedrals_adjacents_EV(target, py_lnhEV, qdih, nsteps, nsamples):
    cdef cgraph.Graph target_graph = cgraph.Graph(target)
    cdef cgraph.Graph fragment_graph = cgraph.Graph(target)

    cdef cgraph.histogramEV_t * lnhEV = cgraph.histogramEV_alloc(target_graph._graph)
    if lnhEV == NULL: raise MemoryError()
    histogramEV_py_to_c(py_lnhEV, lnhEV)

    dihedrals = EVHistogram(target)
    adjacents = EVHistogram(target)
    visits = EVHistogram(target, dtype=int)

    cdef int E, V, B

    for sample in range(nsamples):
        cgraph.mc_sample_subgraphs_EV(target_graph._graph, fragment_graph._graph, lnhEV, nsteps)
        E = fragment_graph._graph.nedges
        V = fragment_graph._graph.nvertices
        B = fragment_graph._graph.nbridges
        dih = qdih**(1 + B - V)
        dihedrals.inc(E, V, dih)
        adjacents.inc(E, V, fragment_graph._graph.nadjacents)
        visits.inc(E, V, 1)

    for i in range(visits.h.shape[0]):
        for j in range(visits.h.shape[1]):
            if visits.h[i,j] > 0:
                dihedrals.h[i,j] = -math.log(dihedrals.h[i,j] / visits.h[i,j])
                adjacents.h[i,j] /= visits.h[i,j]
            else:
                dihedrals.h[i,j] = -1.
                adjacents.h[i,j] = -1.

    cgraph.histogramEV_free(lnhEV)
    return dihedrals, adjacents, visits

def mc_simulate_energies_EV(target, py_lnhEV, energy_dicts, nsteps, nsamples):
    cdef cgraph.Graph target_graph = cgraph.Graph(target)
    cdef cgraph.Graph fragment_graph = cgraph.Graph(target)

    cdef cgraph.histogramEV_t * lnhEV = cgraph.histogramEV_alloc(target_graph._graph)
    if lnhEV == NULL: raise MemoryError()
    histogramEV_py_to_c(py_lnhEV, lnhEV)

    visits = EVHistogram(target, dtype=int)
    sums = [EVHistogram(target) for k in range(len(energy_dicts))]
    energy_means = [sum(energy_dicts[k].values()) / len(energy_dicts[k]) for k in range(len(energy_dicts))]
    energy_diffs = [{edge : (energy_dicts[k][edge] - energy_means[k]) for edge in energy_dicts[k].keys()} \
                    for k in range(len(energy_dicts))]

    cdef size_t E, V

    for sample in range(nsamples):
        cgraph.mc_sample_subgraphs_EV(target_graph._graph, fragment_graph._graph, lnhEV, nsteps)
        E = fragment_graph._graph.nedges
        V = fragment_graph._graph.nvertices
        for k in range(len(energy_dicts)):
            val = math.exp(sum(energy_diffs[k][edge] for edge in fragment_graph.edges_iter()))
            sums[k].inc(E, V, val)
        visits.inc(E, V, 1)

    means = EVHistogram(target)
    for i in range(visits.h.shape[0]):
        for j in range(visits.h.shape[1]):
            if visits.h[i,j] > 0:
                E = i
                means.h[i,j] = sum(math.log(sums[k].h[i,j] / visits.h[i,j]) / E + energy_means[k] \
                                   for k in range(len(energy_dicts))) / len(energy_dicts)
            else:
                means.h[i,j] = -1.

    cgraph.histogramEV_free(lnhEV)
    return means, visits
