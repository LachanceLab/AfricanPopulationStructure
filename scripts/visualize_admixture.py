#!/usr/bin/env python
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import sys
import argparse
import seaborn as sns
import numpy.ma as ma
from pykrige.ok import OrdinaryKriging
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep
from itertools import permutations


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


def sort_individuals(df, K):
    """
    Sort individuals by modal ancestry components and it's intensity
    :param df: pd.DataFrame, data frame containing ancestry porportions
    :param K: int, K used to run admixture
    :return: pd.DataFrame, sorted df
    """
    # get mean ancestry per component for each population
    mean_ancestry_pops = df.groupby('population').mean()
    # identify modal component for each population, i.e., the component with the highest ancestry proportion
    modal_clusters_pops = mean_ancestry_pops.idxmax(axis=1).sort_values()
    individual_orders = []
    # iterate over all k
    for i in range(K):
        # get all populations for which the current component is model
        pops = modal_clusters_pops[modal_clusters_pops == i].index.values.tolist()
        # sort the population by they mean ancestry proportions of the modal component
        pops_sorted = mean_ancestry_pops.loc[pops, i].sort_values(ascending=False).index.values
        inds = []
        # sort the individuals within a population by ancestry proportion
        for pop in pops_sorted:
            inds.extend(df.loc[df.population == pop, i].sort_values(ascending=False).index.values.tolist())
        individual_orders.extend(inds)

    df = df.loc[individual_orders, :]
    return df


def get_best_k_and_plot_best_k(qfiles, Kvals, fam, log_dir, colors, path_to_coords, vmin, prefix):
    """
    Find best K across different ADMIXTURE runs using different K and plot admixture plot of best K separately
    :param qfiles: list, list of paths to aligned .Q files
    :param Kvals: list-like, K values that were used
    :param fam: df, pd.DataFrame of plink .fam file
    :param log_dir: str, path to directory with log files of ADMIXTURE runs
    :param colors: list, colors used for plotting
    :param path_to_coords: str, path to file with individual geographic coordinates
    :param vmin: float, minimum ancestry proportion to be considered during Kriging interpolation
    :param prefix: str, prefix of output figure "_admixture_plot_best_K.pdf" will be appended
    :return: array-like, array-like; order of individuals, cross-validation error
    """
    min_k = min(Kvals)
    # get CV error of ADMIXTURE runs
    cv_error = np.zeros(len(Kvals))
    logfiles = glob.glob(log_dir + '*.log')
    # read CV error
    for log in logfiles:
        k = int(log.split('.')[-2])
        with open(log, 'r') as fp:
            for line in fp:
                if line.startswith("CV error"):
                    cv_error[k - min_k] = float(line.split(':')[1])
                    break
            fp.close()
    # get best k
    best_k = Kvals[np.argmin(cv_error)]
    df = pd.read_csv(qfiles[np.argmin(cv_error)], sep=' ', header=None, index_col=None)
    df = df.iloc[:, 6:]
    df.columns = [i for i in range(best_k)]
    df['population'] = fam.fid.values
    # sort individuals
    df_sorted = sort_individuals(df, best_k)
    # get order of idnviduals
    individual_order = df_sorted.index.values
    plot_best_k(df_sorted, best_k, min(cv_error), colors, prefix)
    # plot_krigin_best_k(df, best_k, path_to_coords, vmin, colors, prefix)
    return individual_order, cv_error


def plot_krigin_best_k(df, K, path_to_coords, vmin, colors, prefix):
    """
    Interpolate admixture proportion in space using Kriging method
    :param df: pd.DataFrame, sorted DataFrame with ancestry proportions
    :param K: int, K used for admixture run
    :param path_to_coords:
    :param vmin: float, minimum admixture proportion required for Kriging method
    :param colors: list, colors used for plotting
    :param prefix: str, prefix for figure
    """
    # The parameters are tailored to our specific case and may not give meaningful results in other cases
    if K > 4:
        raise AttributeError("Only K<= 4 is supported")
    land_shp_fname = shpreader.natural_earth(resolution='50m',
                                             category='physical', name='land')

    land_geom = unary_union(list(shpreader.Reader(land_shp_fname).geometries()))
    # create land object
    land = prep(land_geom)

    # read data
    coords = pd.read_csv(path_to_coords, sep='\t', header=None, names=['fid', 'lat', 'long'], usecols=[0, 2, 3])
    coords.drop_duplicates('fid', inplace=True)
    df = df.join(coords.set_index('fid'), on='population')

    # longitude and latidude coords of interest
    long_coords = np.arange(-25, 60, 0.15, dtype=float)
    lat_coords = np.arange(-40, 50, 0.15, dtype=float)
    # color maps
    cmaps = [
        sns.color_palette(f"blend:white,{colors[0]}", as_cmap=True),
        sns.color_palette(f"blend:white,{colors[1]}", as_cmap=True),
        sns.color_palette(f"blend:white,{colors[2]}", as_cmap=True),
        sns.color_palette(f"blend:white,{colors[3]}", as_cmap=True),
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
    fig.savefig(f"{prefix}_admixture_kriging_best_K.pdf", bbox_inches='tight', dpi=600)


def plot_best_k(df, K, cv, colors, prefix):
    """
    Create horizontal bar plot of ADMIXTURE results
    :param df: pd.DataFrame, sorted DataFrame with ancestry proportions
    :param K: int, K used for admixture run
    :param colors: list, list of colors to use
    :param cv: float, CV error of ADMIXTURE run
    :param prefix: str, prefix for figure
    """
    fig, ax = plt.subplots(figsize=(2, 8))
    populations = df.population.values
    y_coords = [0]
    y_ticks = []
    y_labels = []
    prev_pop = populations[0]
    pop_coord_start = 0
    cumulative_padding = 0
    padding = 0.15
    for i, pop in enumerate(populations[1:], 1):
        if prev_pop != pop:
            # get ytick positions and labels
            y_ticks.append((pop_coord_start + i + cumulative_padding - 1) / 2)
            y_labels.append(prev_pop)
            prev_pop = pop
            # add white space between populations
            cumulative_padding += padding
            pop_coord_start = i + cumulative_padding
        # determine y coords --> add white space between populations
        y_coords.append(i + cumulative_padding)
    y_ticks.append((pop_coord_start + i + cumulative_padding) / 2)
    y_labels.append(prev_pop)
    prev_k = []
    # do plotting
    for i in range(K):
        if i > 0:
            ax.barh(y_coords, df.loc[:, i], left=df.loc[:, prev_k].sum(axis=1), height=1, color=colors[i])
        else:
            ax.barh(y_coords, df.loc[:, i], height=1, color=colors[i])
        prev_k.append(i)
    # formatting
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([x.replace('Sabue', 'Chabu') for x in y_labels], fontsize=5)
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.invert_yaxis()
    ax.set_xlabel(f"K={K}\nCV={round(cv, 4)}", fontsize=8)
    ax.set_ylim([len(populations) + cumulative_padding - 0.5, -0.5])
    ax.set_xlim([0, 1])
    fig.savefig(f"{prefix}_admixture_plot_best_K.pdf", bbox_inches='tight', dpi=600)


def plot_all_k(qfiles, Kvals, cv_error, fam, individual_order, colors, prefix):
    """
    Plot admixture proportions across all K
    :param qfiles: list, list of paths to aligned .Q files
    :param Kvals: list-like, K values that were used
    :param cv_error: array-like, cross-validation errors of individual runs
    :param fam: pd.DataFrame, df of plink .fam file
    :param individual_order: list, order of individudals in which to plot them to guarantee consistency
    :param colors: list, colors used for plotting
    :param prefix: str, prefix of output figure "_admixture_plots.pdf" will be appended
    """
    fig, ax = plt.subplots(1, len(qfiles), figsize=(8, 11))
    plt.subplots_adjust(wspace=0.15)
    min_k = min(Kvals)
    for qfile, K, cv in zip(qfiles, Kvals, cv_error):
        df = pd.read_csv(qfile, sep=' ', header=None, index_col=None)
        df = df.iloc[:, 6:]
        df.columns = [i for i in range(K)]
        df['population'] = fam.fid.values
        df_sorted = df.loc[individual_order]
        if cv == min(cv_error):
            best_k = True
        else:
            best_k = False
        plot_admixture_proportions(df_sorted, K, ax[K - min_k], cv, colors, min_k, best_k)
    fig.savefig(prefix + "_admixture_plots.pdf", bbox_inches='tight', dpi=600)
    plt.close()


def plot_admixture_proportions(df, K, ax, cv, colors, min_k, best_k=False):
    """
    Helper to plot admixture proportions of for specific value of K
    :param df: pd.DataFrame, sorted DataFrame with ancestry proportions
    :param K: int, K used for admixture run
    :param ax: ax object to plot on
    :param cv: float, CV error of ADMIXTURE run
    :param colors: list, list of colors to use
    :param min_k: int, minimum k that was used in any run
    :param best_k: boolean, indicating whether it is the best K, if yes plot x label in bold
    :return: ax
    """
    populations = df.population.values
    y_coords = [0]
    y_ticks = []
    y_labels = []
    prev_pop = populations[0]
    pop_coord_start = 0
    cumulative_padding = 0
    padding = 0.15
    for i, pop in enumerate(populations[1:], 1):
        if prev_pop != pop:
            # get ytick positions and labels
            y_ticks.append((pop_coord_start + i + cumulative_padding - 1) / 2)
            y_labels.append(prev_pop)
            prev_pop = pop
            # add white space between populations
            cumulative_padding += padding
            pop_coord_start = i + cumulative_padding
        # determine y coords --> add white space between populations
        y_coords.append(i + cumulative_padding)
    y_ticks.append((pop_coord_start + i + cumulative_padding) / 2)
    y_labels.append(prev_pop)
    prev_k = []
    # do plotting
    for i in range(K):
        if i > 0:
            ax.barh(y_coords, df.loc[:, i], left=df.loc[:, prev_k].sum(axis=1), height=1, color=colors[i])
        else:
            ax.barh(y_coords, df.loc[:, i], height=1, color=colors[i])
        prev_k.append(i)
    # formatting
    if K > min_k:
        ax.set_yticks([])
    else:
        ax.set_yticks(y_ticks)
        ax.set_yticklabels([x.replace('Sabue', 'Chabu') for x in y_labels], fontsize=5)
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.invert_yaxis()
    if best_k:
        ax.set_xlabel(f"K={K}\nCV={round(cv, 4)}", fontsize=5, fontweight='bold')
    else:
        ax.set_xlabel(f"K={K}\nCV={round(cv, 4)}", fontsize=5)
    ax.set_ylim([len(populations) + cumulative_padding - 0.5, -0.5])
    ax.set_xlim([0, 1])
    return ax


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_dir',
                        help='Input directory of "converted" files, i.e., .Q files aligned with CLUMPAK')
    parser.add_argument('-l', '--logfile_dir',
                        help='Directory containing log files from admixture runs to extract CV error.')
    parser.add_argument('-f', '--fam', help='Plink fam file to extract individual and family IDs')
    parser.add_argument('-c', '--coords', help='File with individual geographic coordinates')
    parser.add_argument('--vmin', type=float,
                        help='Only admixture proportions greater vmin are considered for Kriging method. default=0.5',
                        default=0.5)
    parser.add_argument('-o', '--output_prefix', help='Output prefix')
    args = parser.parse_args()
    input_dir = args.input_dir
    log_dir = args.logfile_dir
    if not input_dir.endswith('/'):
        input_dir += '/'
    if not log_dir.endswith('/'):
        log_dir += '/'
    # get input q files
    aligned_q_files = glob.glob(input_dir + "*.converted")
    # read fam file
    fam = pd.read_csv(args.fam, sep='\t', names=['fid', 'iid'], usecols=[0, 1])
    # get k values
    Kvals = [int(f.split('.')[-3]) for f in aligned_q_files]
    # sort based on K
    aligned_q_files = np.array(aligned_q_files)[np.argsort(Kvals)]
    Kvals = np.sort(Kvals)
    max_k = max(Kvals)

    # initialize colors
    colors = ["purple", "orange", "green", "blue"]
    # pad colors
    while max_k > len(colors):
        colors.append((np.random.random(), np.random.random(), np.random.random()))

    individual_order, cv_error = get_best_k_and_plot_best_k(aligned_q_files, Kvals, fam, log_dir, colors, args.coords,
                                                            args.vmin, args.output_prefix)
    plot_all_k(aligned_q_files, Kvals, cv_error, fam, individual_order, colors, args.output_prefix)


if __name__ == '__main__':
    main(sys.argv[1:])
