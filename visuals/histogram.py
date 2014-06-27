#!/usr/bin/python

''' Simple script that creates histogram.
'''

import numpy as np
import matplotlib.pyplot as plt

import sys, os

def create_histograms(cds_data_file):

    # Get file contents
    lines = [line.rstrip('\n') for line in open(cds_data_file)]

    # Extract interesting columns - mean_cov/CDS_length
    mean_cov_data   = []
    length_data     = []
    for line in lines:
        row_data = line.split()
        mean_cov    = float(row_data[3])
        length      = int(row_data[4])

        mean_cov_data.append(mean_cov)
        length_data.append(length)

    # Create histograms

    plt.hist(length_data, bins=200)
    plt.title("all CDSs length")
    plt.show()

if __name__ == '__main__':
    cds_data_file = sys.argv[1]

    create_histograms(cds_data_file)

