#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
import math

import cairo
import pylab

import sys, os

def apply_to_folder(folder_path, img_folder_path):

    # Create folder if does not exist already
    if not os.path.exists(img_folder_path):
        os.makedirs(img_folder_path)

    # Create and store graphs
    for filename in os.listdir(folder_path):
        img_name = os.path.splitext(filename)[0] + ".png"

        data_path = os.path.join(folder_path, filename)
        img_path  = os.path.join(img_folder_path, img_name)

        print(filename)
        create_plot(data_path, img_path)

def create_plot(data_path, img_path):

    # Get file contents
    lines = [line.rstrip('\n') for line in open(data_path)]


    # Extract data from header
    mean    = float(lines.pop(0))
    length  = int(lines.pop(0))
    gene    = lines.pop(0)
    product = lines.pop(0)

    # Info - textbox
    # plt.text(0.5, 0.5, 'probni tekst', transform=plt.transAxes)

    # Load coverage data
    x = []
    y = []
    for line in lines:
        row_data = line.split()
        x.append(int(row_data[0]))
        y.append(int(row_data[1]))

    # ------------ Plotting ------------- #

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Mean
    ax.axhline(y=mean, color='r')
    # Coverage
    ax.plot(x, y)
    # Info textbox
    info = "length:\t{0}\ngene:\t{1}\nprod:\t{2}\n".format(length, gene, product)
    info = info.expandtabs()

    plt.text(0.1, 0.1, info, 
             color='white',
             bbox=dict(facecolor='black', alpha=0.6, edgecolor='white'),
             transform=ax.transAxes)

    # plt.show()
    plt.savefig(img_path, format='png')

    # Clear figure before next plotting
    plt.clf()
    plt.close(fig)

if __name__ == '__main__':
    folder_path     = sys.argv[1]
    img_folder_path = sys.argv[2]

    apply_to_folder(folder_path, img_folder_path)

