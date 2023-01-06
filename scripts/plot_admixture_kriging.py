#!/usr/bin/env python3
import argparse
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import numpy.ma as ma
from pykrige.ok import OrdinaryKriging
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep
from itertools import permutations

# The parameters are tailored to our specific case and may not give meaningful results in other cases


def is_land(x, y, land):
    """
    Function to determine if a given point is on land or not
    :param x: float, longitude
    :param y: float, latitude
    :param land:
    :return: boolean
    """
    return land.contains(sgeom.Point(x, y))


def build_mask(x, y, land):
    """
    Create a mask of points that are on land and should be plotted
    :param x: np.array, longitude values
    :param y: np.array, latitude values
    :param land:
    :return: np.array, boolean
    """
    # mesh
    xv, yv = np.meshgrid(x, y, indexing='xy')
    # initialize mask
    mask = np.zeros(xv.shape, dtype=bool)
    for i in range(x.shape[0]):
        for j in range(y.shape[0]):
            # check if point is land
            mask[j, i] = is_land(xv[j, i], yv[j, i], land)
    return mask


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--qfile', help='List of .Q files to consider')
    parser.add_argument('-f', '--fam', help='.fam file to get populations')
    parser.add_argument('-p', '--path', help='Path to save Figure to. default=results/', default='results/')
    parser.add_argument('--prefix',
                        help='Prefix for figure file name. Figure will be saved at {path}{prefix}_admixture_krigin_k{k}.pdf')
    parser.add_argument('-c', '--coords', help='File with individual geographic coordinates')
    parser.add_argument('-k', type=int, help='K to plot')
    parser.add_argument('--vmin', type=float, help='Only admixture proportions greater vmin. default=0.5', default=0.5)
    args = parser.parse_args()
    data_file = args.qfile
    coords = args.coords
    fam_file = args.fam
    K = args.k
    vmin = args.vmin
    if K > 4:
        raise AttributeError("Only K<= is supported")
    land_shp_fname = shpreader.natural_earth(resolution='50m',
                                             category='physical', name='land')

    land_geom = unary_union(list(shpreader.Reader(land_shp_fname).geometries()))
    # create land object
    land = prep(land_geom)

    # read data
    df = pd.read_csv(data_file, sep=' ', header=None)
    fam = pd.read_csv(fam_file, sep='\t', names=['fid', 'iid'], usecols=[0, 1])
    df['population'] = fam.fid.values
    coords = pd.read_csv(coords, sep='\t', header=None, names=['fid', 'lat', 'long'], usecols=[0, 2, 3])
    coords.drop_duplicates('fid', inplace=True)
    df = df.join(coords.set_index('fid'), on='population')

    # longitude and latidude coords of interest
    long_coords = np.arange(-25, 60, 0.15, dtype=float)
    lat_coords = np.arange(-40, 50, 0.15, dtype=float)
    # color maps
    cmaps = [
        sns.color_palette("blend:white,blue", as_cmap=True),
        sns.color_palette("blend:white,purple", as_cmap=True),
        sns.color_palette("blend:white,green", as_cmap=True),
        sns.color_palette("blend:white,orange", as_cmap=True),
    ]

    projection = ccrs.PlateCarree()

    fig = plt.figure(figsize=(9, 16))
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    ax.set_extent([-25, 60, -40, 50], crs=ccrs.PlateCarree())

    # run kriging for each K
    perms = set(permutations([i for i in range(K)]))
    for perm in perms:
        for i in perm:
            cmap = cmaps[i]
            proportions = df.loc[:, i].values
            ok = OrdinaryKriging(df.long.values,
                                 df.lat.values,
                                 proportions, variogram_model='linear',
                                 coordinates_type='geographic', nlags=2)
            z, ss = ok.execute('grid', long_coords, lat_coords)
            # get mask
            mask = build_mask(long_coords, lat_coords, land)
            # normalize data
            data = (z.data - z.data.min()) / (z.data.max() - z.data.min())
            # update mask
            mask[data < vmin] = False
            masked_data = ma.masked_array(data, ~mask)
            # plot
            ax.imshow(np.flip(masked_data, 0), extent=[-25, 60, -40, 50], cmap=cmap,
                      interpolation='none', vmax=1, vmin=vmin, alpha=1 / len(perms))
    # formatting
    ax.add_feature(cfeature.COASTLINE, color='black')
    ax.add_feature(cfeature.BORDERS, linestyle=':', color='black')
    # plot populations
    ax.scatter(coords.long, coords.lat, color='black', s=8)
    ax.axis('off')
    fig.savefig(f"{args.path}{args.prefix}_admixture_kriging_k{K}.pdf", bbox_inches='tight', dpi=600)


if __name__ == '__main__':
    main(sys.argv[1:])
