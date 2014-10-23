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

import numpy as np

class EVHistogram(object):
    def __init__(self, graph, dtype=float):
        if graph != None:
            Emax = len(graph.edges())
            Vmax = len(graph.nodes())
            self.h = np.zeros((Emax + 1, Emax - Vmax + 2), dtype=dtype)
        else:
            self.h = None

    def _index(self, E, V):
        index = E, (E - V + 1)
        if index[1] < 0:
            raise Exception("The graph is disconnected.")
        return index

    def get(self, E, V):
        return self.h[self._index(E, V)]

    def inc(self, E, V, w=1.):
        self.h[self._index(E, V)] += w

    def write(self, stream, addrow=True):
        if all(self.h[i,-1] == -1. for i in range(self.h.shape[0])):
            addrow = False
        for i in range(self.h.shape[0]):
            for j in range(self.h.shape[1]):
                stream.write("%d %d %g\n" % (i, j, self.h[i,j]))
            if addrow: # Add extra row for gnuplot: `pm3d map corners2color c3'
                stream.write("%d %d %g\n" % (i, self.h.shape[1], -1.))
            stream.write("\n")

    @staticmethod
    def read(stream):
        data = []
        for line in stream:
            if len(line) == 1 or line[0] == '#':
                continue
            data.append([float(s) for s in line.split()])
        hist = EVHistogram(None)
        hist.h = np.zeros((int(max(d[0] for d in data)) + 1, int(max(d[1] for d in data)) + 1))
        for d in data:
            hist.h[int(d[0]),int(d[1])] = d[2]
        return hist

    def __add__(self, s):
        return self.h + s.h
    def __sub__(self, s):
        return self.h - s.h

def read_histogram(stream):
    return EVHistogram.read(stream).h

def read_meta(stream):
    tags = {line.split()[1] : line.split()[2:] for line in stream if line[0] == '#'}
    return tags
