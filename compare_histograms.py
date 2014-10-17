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
