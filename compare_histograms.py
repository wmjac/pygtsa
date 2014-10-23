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
import numpy as np
from pygtsa.histogram import EVHistogram

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('hist1', type=str, help="path to histogram file 1")
    parser.add_argument('hist2', type=str, help="path to histogram file 2")
    clargs = parser.parse_args()

    with open(clargs.hist1, 'r') as f:
        hist1 = EVHistogram.read(f)
    with open(clargs.hist2, 'r') as f:
        hist2 = EVHistogram.read(f)

    diff = hist2 - hist1

    print("mean(diff):", np.mean(diff))
    print("stddev(diff):", np.std(diff))
    print("max(fabs(diff)):", np.max(np.fabs(diff)))
