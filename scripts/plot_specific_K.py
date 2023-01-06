#!/usr/bin/env python3
import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


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


def plot_admixture_proportions(df, K, prefix, figure_path='results/'):
    """
    Create horizontal bar plot of ADMIXTURE results
    :param df: pd.DataFrame, sorted DataFrame with ancestry proportions
    :param K: int, K used for admixture run
    :param prefix: str, prefix for figure
    :param figure_path: str, path where to save figure
    """
    colors = ["blue", "purple", "green", 'orange']
    # pad colors
    while K > len(colors):
        colors.append((np.random.random(), np.random.random(), np.random.random()))
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
    ax.set_yticklabels(y_labels, fontsize=4)
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.invert_yaxis()
    ax.set_xlabel(f"K={K}", fontsize=8)
    ax.set_ylim([len(populations) + cumulative_padding - 0.5, -0.5])
    ax.set_xlim([0, 1])
    fig.savefig(f"{figure_path}{prefix}_admixture_plot_k{K}.pdf", bbox_inches='tight', dpi=600)


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--qfile', help='List of .Q files to consider')
    parser.add_argument('-f', '--fam', help='.fam file to get populations')
    parser.add_argument('-p', '--path', help='Path to save Figure to. default=results/', default='results/')
    parser.add_argument('--prefix',
                        help='Prefix for figure file name. Figure will be saved at {path}{prefix}_admixture_plot_k{k}.pdf')
    parser.add_argument('-k', type=int, help='K to plot')

    args = parser.parse_args()
    figure_path = args.path
    if not figure_path.endswith('/'):
        figure_path += '/'
    df = pd.read_csv(args.qfile, sep=' ', header=None)
    fam = pd.read_csv(args.fam, sep='\t', names=['fid', 'iid'], usecols=[0, 1])
    df['population'] = fam.fid.values
    df_sorted = sort_individuals(df, args.k)
    plot_admixture_proportions(df_sorted, args.k, args.prefix, figure_path)


if __name__ == '__main__':
    main(sys.argv[1:])
