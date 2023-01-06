#!/usr/bin/env python3
"""
Code partially taken from
https://github.com/NovembreLab/feems/blob/main/docsrc/notebooks/cross-validation.ipynb
"""

import numpy as np
import pkg_resources
from sklearn.impute import SimpleImputer
from pandas_plink import read_plink
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from feems.utils import prepare_graph_inputs
from feems import SpatialGraph, Viz
from feems.cross_validation import run_cv
import sys
import argparse
import pandas as pd

# change matplotlib fonts
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = "Arial"


def create_spatial_graph(data_path, prefix, coords_file):
    """
    Create SpatialGraph object
    :param data_path: str, data directory
    :param prefix: str, prefix for plink files in data directory
    :param coords_file: str, path to file with sampling coordinates
    :return: SpatialGraph
    """
    data_path_feems = pkg_resources.resource_filename("feems", "data/")
    # read the genotype data and mean impute missing data
    (bim, fam, G) = read_plink("{}/{}".format(data_path, prefix))
    imp = SimpleImputer(missing_values=np.nan, strategy="mean")
    genotypes = imp.fit_transform((np.array(G)).T)

    # setup graph
    coords = pd.read_csv(coords_file, sep='\t', header=None,
                         names=['fid_c', 'iid_c', 'latitude', 'longitude'])
    coords = fam.set_index('iid').join(coords.set_index('iid_c')).loc[:, ["longitude", "latitude"]].values
    grid_path = "{}/grid_100.shp".format(data_path_feems)  # path to discrete global grid

    # graph input files
    _, edges, grid, _ = prepare_graph_inputs(coord=coords,
                                             ggrid=grid_path,
                                             translated=False,
                                             buffer=0,)

    # # construct spatial graph object
    sp_graph = SpatialGraph(genotypes, coords, grid, edges, scale_snps=True)
    return sp_graph


def cross_validate(sp_graph):
    """
    Perform cross validation to find optimal lambda
    :param sp_graph: SpatialGraph
    :return: float, lambda
    """
    # reverse the order of lambdas and alphas for warmstart
    lamb_grid = np.geomspace(1e-6, 1e2, 20)[::-1]

    # run cross-validation
    cv_err = run_cv(sp_graph, lamb_grid, n_folds=sp_graph.n_observed_nodes, factr=1e10)

    # average over folds
    mean_cv_err = np.mean(cv_err, axis=0)

    # argmin of cv error
    lamb_cv = float(lamb_grid[np.argmin(mean_cv_err)])
    return lamb_cv


def fit_final_graph(sp_graph, lamb_cv, prefix, output_path):
    """
    Fit final SpatialGraph using optimal lambda inferred from cross-validation
    :param sp_graph: SpatialGraph
    :param lamb_cv: float, lambda
    :param output_path: str, directory where to save output figure named feems_plot.pdf
    """
    sp_graph.fit(lamb_cv)
    projection = ccrs.PlateCarree()
    # Plot the FEEMS result
    fig = plt.figure(figsize=(9, 16), dpi=600)
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    ax.set_extent([-25, 60, -40, 50], crs=ccrs.PlateCarree())
    v = Viz(ax, sp_graph, projection=projection, edge_width=.5,
            edge_alpha=1, edge_zorder=100, sample_pt_size=20,
            obs_node_size=7.5, sample_pt_color="black",
            cbar_font_size=5, cbar_ticklabelsize=4)
    v.draw_map()
    v.draw_edges(use_weights=True)
    v.draw_obs_nodes(use_ids=False)
    v.draw_edge_colorbar()
    fig.savefig(output_path + prefix + '_feems_plot.pdf', bbox_inches='tight', dpi=600)
    

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--data_path', help='Data directory. default=./data/', default='./data/')
    parser.add_argument('-p', '--prefix', help='Prefix to plink files')
    parser.add_argument('-c', '--coords', help='Path to file coordinates of samples. '
                                               'default=data/population_geo_coords.tab',
                        default='data/population_geo_coords.tab')
    parser.add_argument('-o', '--output_path', help='Output directory. Will save Figure called feems_plot.pdf '
                                                    'in this direcotry. default=./results/', default='./results/')
    args = parser.parse_args()
    data_path = args.data_path
    prefix = args.prefix
    coords = args.coords
    output_path = args.output_path
    if not data_path.endswith('/'):
        data_path += '/'
    if not output_path.endswith('/'):
        output_path += '/'
    sp_graph = create_spatial_graph(data_path, prefix, coords)
    lamb_cv = cross_validate(sp_graph)
    fit_final_graph(sp_graph, lamb_cv, prefix, output_path)


if __name__ == '__main__':
    main(sys.argv[1:])
