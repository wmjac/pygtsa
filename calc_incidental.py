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

import argparse, math
import numpy as np
from pygtsa.histogram import EVHistogram, read_histogram, read_meta
from pygtsa.calc_yield import load_dihedrals, find_target_indices, lnhzEV

def multivalentQ(n1, n2, w, qdih, qcoord, frac_attractive=1., maxk=None):
    n1 = int(n1 * frac_attractive)
    n2 = int(n2 * frac_attractive)
    Q = 0.
    if maxk == None:
        maxk = min(n1, n2)
    else:
        maxk = min(n1, n2, maxk)
    for k in range(1, maxk + 1):
        prod = 1.
        for r in range(n1 - k + 1, n1 + 1):
            prod *= r
        for r in range(1, k + 1):
            prod /= r
        for r in range(n2 - k + 1, n2 + 1):
            prod *= r
        Q += (prod * k * w) / (qdih**(2 * k - 1) * qcoord)
    return Q

def incidental_VV(lnhz, adjacents, rho, w, qdih, qcoord, frac_attractive=1.):
    N = lnhz.shape[0] - lnhz.shape[1] + 2
    Zin = np.zeros((N + 1, N + 1))
    for i in range(lnhz.shape[0]):
        for j in range(lnhz.shape[1]):
            E1 = i
            V1 = E1 - j + 1
            for ii in range(lnhz.shape[0]):
                for jj in range(lnhz.shape[1]):
                    E2 = ii
                    V2 = E2 - jj + 1
                    if np.isfinite(lnhz[i,j]) and np.isfinite(lnhz[ii,jj]) and adjacents[i,j] > 0 and adjacents[ii,jj] > 0:
                        lnhzz = lnhz[i,j] + lnhz[ii,jj]
                        Q = multivalentQ(adjacents[i,j], adjacents[ii,jj], w, qdih, qcoord, frac_attractive=frac_attractive)
                        if i == j and ii == jj:
                            Zin[V1,V2] += Q * rho * math.exp(lnhzz) / 2.
                        else:
                            Zin[V1,V2] += Q * rho * math.exp(lnhzz)
    return Zin

def on_off_pathway_ratio(lnhz, incidental):
    Zid = 0.
    for i in range(lnhz.shape[0]):
        for j in range(lnhz.shape[1]):
            if np.isfinite(lnhz[i,j]):
                V = i - j + 1
                Zid += math.exp(lnhz[i,j])
    Zin = 0.
    for i in range(incidental.shape[0]):
        for j in range(incidental.shape[1]):
            Zin += incidental[i,j]
    return Zid / Zin

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('lnhEV', type=str, help="path to lnhEV file")
    parser.add_argument('dihedrals', type=str, help="path to dihedrals file")
    parser.add_argument('rho', type=float, help="dimensionless number density")
    parser.add_argument('adjacents', type=str, help="path to adjacents file")
    parser.add_argument('w', type=float, help="mean dimensionless incidental energy")
    parser.add_argument('--output', type=str, default='incidental.dat', help="path to output file [incidental.dat]")
    subparsers = parser.add_subparsers(dest='variance')
    parser_zero = subparsers.add_parser('zero-var')
    parser_const = subparsers.add_parser('const-var')
    for p in (parser_zero, parser_const):
        p.add_argument('E', type=float, help="dimensionless energy")
    parser_const.add_argument('distribution', type=str, help="path to energy distribution file")
    clargs = parser.parse_args()

    if clargs.variance == None:
        parser.print_usage()
        print("error: variance mode not specified")
        raise SystemExit

    with open(clargs.lnhEV, 'r') as f:
        lnhEV = read_histogram(f)
    with open(clargs.dihedrals, 'r') as f:
        dihedrals, dihedrals_meta = load_dihedrals(f)
    target_indices = find_target_indices(dihedrals)
    with open(clargs.adjacents, 'r') as f:
        adjacents = read_histogram(f)

    if clargs.variance == 'const-var':
        with open(clargs.distribution, 'r') as f:
            epsilon_expmean = read_histogram(f)
        epsilon_expmean -= epsilon_expmean[max(target_indices)]
    else:
        epsilon_expmean = np.zeros(dihedrals.shape)

    mu = math.log(clargs.rho)
    lnhz = lnhzEV(lnhEV, dihedrals, mu, clargs.E, epsilon_expmean, qcoord=dihedrals_meta['qcoord'])
    incidental = incidental_VV(lnhz, adjacents, clargs.rho, clargs.w, dihedrals_meta['qdih'], dihedrals_meta['qcoord'])

    on_off_ratio = on_off_pathway_ratio(lnhz, incidental)
    print("On/off pathway ratio (Zid/Zin): %g" % on_off_ratio)

    print("Writing output to %s" % clargs.output)
    with open(clargs.output, 'w') as f:
        f.write("# epsilon = %g\n" % clargs.E)
        f.write("# on/off pathway ratio = %g\n" % on_off_ratio)
        f.write("# V1 V2 Zin(V1,V2)\n")
        for i in range(1, incidental.shape[0]):
            for j in range(1, incidental.shape[1]):
                f.write("%g %g %g\n" % (i, j, incidental[i,j]))
            f.write("\n")
