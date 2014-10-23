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

def load_dihedrals(stream):
    dihedrals = read_histogram(stream)
    stream.seek(0, 0)
    dihedrals_meta = read_meta(stream)
    dihedrals_meta['qdih'] = float(dihedrals_meta['qdih'][-1])
    dihedrals_meta['qcoord'] = float(dihedrals_meta['qcoord'][-1])
    return dihedrals, dihedrals_meta

def find_target_indices(dihedrals, tol=1.e-5):
    dih_max = np.max(dihedrals)
    target_indices = []
    for i in range(dihedrals.shape[0]):
        for j in range(dihedrals.shape[1]):
            if math.fabs(dihedrals[i,j] - dih_max) < tol:
                target_indices.append((i,j))
    return target_indices

def lnhzEV(lndos, dihedrals, mu, epsilon_mean, epsilon_expmean, qcoord=4):
    lnqcoord = math.log(qcoord)
    VG = lndos.shape[0] - lndos.shape[1] + 2
    EG = lndos.shape[0] - 1
    lnhz = np.zeros(lndos.shape)
    for i in range(lndos.shape[0]):
        for j in range(lndos.shape[1]):
            E = i
            V = E - j + 1
            if E == 0 and V == 1:
                # Monomer (trivial subgraph)
                lnhz[i,j] = math.log(VG) + mu
            elif lndos[i,j] >= 0:
                # Connected subgraph
                epsilon = epsilon_mean + epsilon_expmean[i,j]
                lnhz[i,j] = (lndos[i,j]          # Number of subgraphs
                            + V * mu             # Translational entropy
                            - (V - 1) * lnqcoord # Rotational entropy loss (due to association)
                            - dihedrals[i,j]     # Dihedral entropy loss (due to bridge elimination)
                            + E * epsilon)       # Association energy
            else:
                # Impossible structure
                lnhz[i,j] = -np.inf
    return lnhz

def yield_from_lnhz(lnhz, target_indices):
    lnhzG = math.log(sum(math.exp(lnhz[target_index]) for target_index in target_indices))
    Zid = 0.
    for i in range(lnhz.shape[0]):
        for j in range(lnhz.shape[1]):
            if np.isfinite(lnhz[i,j]):
                V = i - j + 1
                Zid += math.exp(lnhz[i,j] - lnhzG)
    return 1. / Zid

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('lnhEV', type=str, help="path to lnhEV file")
    parser.add_argument('dihedrals', type=str, help="path to dihedrals file")
    parser.add_argument('rho', type=float, help="dimensionless number density")
    parser.add_argument('--output', type=str, default='yield.dat', help="path to output file [yield.dat]")
    subparsers = parser.add_subparsers(dest='variance')
    parser_zero = subparsers.add_parser('zero-var')
    parser_const = subparsers.add_parser('const-var')
    for p in (parser_zero, parser_const):
        p.add_argument('Emin', type=float, help="minimum dimensionless energy")
        p.add_argument('Emax', type=float, help="maximum dimensionless energy")
        p.add_argument('dE', type=float, help="dimensionless energy step")
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

    if clargs.variance == 'const-var':
        with open(clargs.distribution, 'r') as f:
            epsilon_expmean = read_histogram(f)
        epsilon_expmean -= epsilon_expmean[max(target_indices)]
    else:
        epsilon_expmean = np.zeros(dihedrals.shape)

    mu = math.log(clargs.rho)
    epsilon = clargs.Emin
    eta = {}
    while epsilon < clargs.Emax + clargs.dE:
        lnhz = lnhzEV(lnhEV, dihedrals, mu, epsilon, epsilon_expmean, qcoord=dihedrals_meta['qcoord'])
        eta[epsilon] = yield_from_lnhz(lnhz, target_indices)
        epsilon += clargs.dE

    print("Writing output to %s" % clargs.output)
    with open(clargs.output, 'w') as f:
        f.write("# epsilon rho yield\n")
        for epsilon in sorted(eta.keys()):
            f.write("%g %g %g\n" % (epsilon, clargs.rho, eta[epsilon]))
