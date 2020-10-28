#!/usr/bin/env python

"""@package docstring
File: gen_fluor_movies.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""

from pathlib import Path
from matplotlib.animation import FFMpegWriter
from time import time
from mpl_toolkits import mplot3d  # for vizualization
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.animation import FuncAnimation
from scipy import special
import math
import yaml
from copy import deepcopy
from numba import jit, njit, vectorize

rng = np.random.default_rng()  # initialize generator instance

SQRT2 = np.sqrt(2)


@vectorize(nopython=True)
def nberf(x):
    return math.erf(x)


class filament():

    def __init__(self, line):
        data = line.split()
        self.gid = data[1]
        dat = np.asarray(data[2:], dtype=np.double)
        self.radius = dat[0]
        self.minus_end = dat[1:4]
        self.plus_end = dat[4:]

        self.vec = self.plus_end - self.minus_end
        self.lengthxy = np.linalg.norm(self.vec[:-1])
        self.length = np.linalg.norm(self.vec)
        self.orientation = self.vec / self.length

        self.theta = -np.arctan2(self.orientation[1], self.orientation[0])


def read_dat_sylinder(fpath):
    # Read a SylinderAscii_X.dat file

    # open the file and read the lines
    with fpath.open('r') as file1:
        filecontent = file1.readlines()

        # Delete the first two lines because they dont have any data
        filecontent[0:2] = []

        # Create list of filaments
        filaments = sorted([filament(line)
                            for line in filecontent], key=lambda x: x.gid)
    return filaments


@jit(nopython=True, parallel=True)
def GaussianLine2D(x, y, A, sigma, x0, y0, L, theta):
    yprime = y0 - y
    xprime = x0 - x
    denom = SQRT2 * sigma

    term0 = (-yprime * np.cos(theta) + xprime * np.sin(-theta)) / denom
    expterm = A * np.exp(-term0 * term0)

    erf_arg = (xprime * np.cos(theta) + yprime * np.sin(-theta)) / denom
    erfterm = nberf((L / denom) + erf_arg) - nberf(erf_arg)

    z = expterm * erfterm
    return z


def make_image_bkg(box_dim, image_params):
    pixelSize = image_params['pixelsize']
    noiseStd = image_params['noisestd']
    bkglevel = image_params['bkglevel']

    numPixelsX = np.int_(np.ceil((box_dim[0]) / pixelSize))
    numPixelsY = np.int_(np.ceil((box_dim[1]) / pixelSize))

    xpixels = np.arange(0, numPixelsX)
    ypixels = np.arange(0, numPixelsY)
    [X, Y] = np.meshgrid(xpixels, ypixels)

    imagedata = (bkglevel * np.ones((numPixelsX, numPixelsY))
                 + np.random.standard_normal((numPixelsX, numPixelsY))
                 * noiseStd)
    return X * pixelSize, Y * pixelSize, imagedata


def draw_2d_gauss_filament(X, Y, image_params, fil):
    A = image_params['A']
    sigma = image_params['sigmaxy']

    imagedata = GaussianLine2D(X, Y, A, sigma,
                               fil.minus_end[0], fil.minus_end[1],
                               fil.lengthxy, fil.theta)
    return imagedata


def get_file_number(path):
    name = path.stem
    num = name.split("_")[-1]
    return int(num)


def count_fils(path):
    with path.open('r') as pf:
        for i, l in enumerate(pf, -1):  # Don't count first two lines
            pass
    return i


def create_fluor_frame(fil_dat_path, fil_idx_arr, run_params,
                       image_params, draw_func=draw_2d_gauss_filament):
    # Get filament data
    filaments = read_dat_sylinder(fil_dat_path)
    # Get dimensions of simulation
    sim_box = np.asarray(run_params['simBoxHigh'])
    # Create background to draw filaments
    X, Y, image_data = make_image_bkg(sim_box, image_params)

    for i in fil_idx_arr:
        image_data += draw_func(X, Y, image_params, filaments[i])

    return image_data


def animate(i, ax, X, Y, frames, vmax):
    ax.clear()
    pcm = ax.pcolormesh(X, Y, frames[i], cmap='gray', vmax=vmax)
    ax.set_title("Frame {}".format(i))
    return [pcm]


##########################################
if __name__ == "__main__":
    image_params = {
        'sigmaxy': np.float_(.007),
        'A': np.float_(10.0),
        'bkglevel': np.float_(0.0),
        'noisestd': np.float_(0.0),
        'pixelsize': np.float_(.01067),
        'graph_frac': np.float_(.1),
    }

    with open('TestConfig.yaml', 'r') as yf:
        run_params = yaml.safe_load(yf)

    sim_box = np.asarray(run_params['simBoxHigh'])

    result_dir = Path("./test_results/")
    fil_dat_paths = sorted(result_dir.glob("SylinderAscii*.dat"),
                           key=get_file_number)

    rng = np.random.default_rng()
    nfils = count_fils(fil_dat_paths[0])

    fil_idx_arr = rng.choice(nfils, int(nfils * image_params['graph_frac']))
    print("Total filaments: {}, filaments graphed: {} ({}%)".format(nfils, int(
        nfils * image_params['graph_frac']), image_params['graph_frac'] * 100))

    print(len(fil_idx_arr))

    X, Y, _ = make_image_bkg(sim_box, image_params)

    frames = []
    for i, fdp in enumerate(fil_dat_paths[::3]):
        t0 = time()
        frames += [create_fluor_frame(fdp,
                                      fil_idx_arr,
                                      run_params,
                                      image_params)]
        print("Frame {} created in: {:.2g} sec".format(i, time() - t0))
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect('equal')
    ani = FuncAnimation(fig, animate, len(frames),
                        fargs=(ax, X, Y, frames, 5. * image_params['A']),
                        blit=True)
    writer = FFMpegWriter(fps=10, bitrate=1800)
    ani.save("test_video.mp4", writer=writer)
