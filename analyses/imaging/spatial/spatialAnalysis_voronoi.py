#!/usr/bin/env python
# coding: utf-8

"""
Spatial cell-cell interaction analysis using Voronoi diagrams.
Use random permutation to calculate z-scores of interaction
frequencies.

Requirements:
    * Python >= 3.7.0
    * scipy >= 1.4.1
    * seaborn >= 0.10.1
    * pandas >= 1.0.3
    * matplotlib >= 3.2.1
    * numpy >= 1.19.0


Copyright (c) 2022 Dietmar Rieder <dietmar.rieder@i-med.ac.at>
MIT License <http://opensource.org/licenses/MIT>

"""

RELEASE = False
__version_info__ = (
    "0",
    "1",
)
__version__ = ".".join(__version_info__)
__version__ += "-dev" if not RELEASE else ""


# In[1]:
# Import libraries
#

import os
import argparse
import numpy as np
import pickle
import pandas as pd
from multiprocessing import Pool
from scipy.spatial import Voronoi
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

import matplotlib
matplotlib.use('agg')

import celltypedefs as ctdef

# In[2]:
# sample assignment
sample_map = {
    "17_8146": "CRC02",
    "17_9065": "CRC03",
    "17_10053": "CRC04",
    "18_02076": "CRC13",
    "18_12804": "CRC26",
    "18_12816": "CRC26LM",
}

# In[10]:

immuneCellIdx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
tumorCellIdx = [26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38]

# In[11]:

# helper function for finding data files
def get_sample_data(dir_name, sample_id):
    # Get the list of all files in directory tree at given path
    sample_data_files = []
    list_of_files = list()
    for (dir_path, dir_names, file_names) in os.walk(dir_name):
        list_of_files += [os.path.join(dir_path, file) for file in file_names]

    for elem in list_of_files:
        data_file = elem.split("/")[-1]
        if data_file.startswith(sample_id):
            sample_data_files.append(elem)

    return sample_data_files


# In[12]:

# register the voronoi based cell cell interactions
def register_voronoi_interactions(vor, df):
    interaction_count = np.zeros((len(ctdef.cellTypeMap), len(ctdef.cellTypeMap)), dtype="uint32")

    for interaction in vor.ridge_points:
        cTa = df["phenotype"][interaction[0]]
        cTb = df["phenotype"][interaction[1]]

        cTaIdx = ctdef.cellTypeMap[cTa]
        cTbIdx = ctdef.cellTypeMap[cTb]
        cTaPos = np.asarray([df["Location_Center_X"][interaction[0]], df["Location_Center_Y"][interaction[0]]])
        cTbPos = np.asarray([df["Location_Center_X"][interaction[1]], df["Location_Center_Y"][interaction[1]]])

        cellDist = np.linalg.norm(cTaPos - cTbPos)

        if cellDist < maxVoronoiDist:
            interaction_count[cTaIdx, cTbIdx] += 1

    return interaction_count


# In[13]:

# random voronoi interaction permutation
def register_rnd_voronoi_interactions(p):
    rnd_voronoi_interaction_count_table = np.zeros((len(ctdef.cellTypeMap), len(ctdef.cellTypeMap)), dtype="uint32")

    # report progress every 10th permutation
    if p % 10 == 0:
        print(f"vor {p}\r", end="")

    # random sample phenotypes
    rnd_phenotypes = df["phenotype"].sample(frac=1).reset_index(drop=True)

    for interaction in vor.ridge_points:
        cTa = rnd_phenotypes[interaction[0]]
        cTb = rnd_phenotypes[interaction[1]]

        cTaIdx = ctdef.cellTypeMap[cTa]
        cTbIdx = ctdef.cellTypeMap[cTb]
        cTaPos = np.asarray([df["Location_Center_X"][interaction[0]], df["Location_Center_Y"][interaction[0]]])
        cTbPos = np.asarray([df["Location_Center_X"][interaction[1]], df["Location_Center_Y"][interaction[1]]])

        # calculate cell-cell distance (center of mass to center of mass)
        cellDist = np.linalg.norm(cTaPos - cTbPos)

        # increase voronoi direct neighbor count if cellDist does not exceed maxVoronoiDist (default 200)
        if cellDist < maxVoronoiDist:
            rnd_voronoi_interaction_count_table[cTaIdx, cTbIdx] += 1

    return rnd_voronoi_interaction_count_table


# In[14]:

# ###   function to calculate the z-score
#
# for p-value calculation we use:
#
# Phipson, B. & Smyth, G. K. Permutation P-values should never be zero:
# calculating exact P-values when permutations are randomly drawn. Stat.
# Appl. Genet. Mol. Biol. 9, Article39 (2010).
#
def calculate_zscore_stat(sum_close_pairs, rnd_close_pairs_counts, num_permutations):
    z = (sum_close_pairs - np.nanmean(rnd_close_pairs_counts)) / np.nanstd(rnd_close_pairs_counts)
    p = stats.norm.sf(abs(z))

    rnd_greater = 0
    if not np.isnan(sum_close_pairs):
        rnd_greater = sum(i > sum_close_pairs for i in rnd_close_pairs_counts)
        rnd_smaller = sum(i < sum_close_pairs for i in rnd_close_pairs_counts)
        pr = (rnd_greater + 1) / (num_permutations + 1)
        pl = (rnd_smaller + 1) / (num_permutations + 1)

    return z, p, pr, pl


# In[15]:
# main
if __name__ == "__main__":

    # Argument parser
    parser = argparse.ArgumentParser(description="Run spatial cell-cell interaction analysis")
    parser.add_argument("--result_dir", required=True, default="", type=str, help="Results directory")
    parser.add_argument(
        "--data_dir",
        required=True,
        default="",
        type=str,
        help='Data directory with phenotype coordinate files: i.e. path to "Single_cell_data" directory',
    )
    parser.add_argument(
        "--run_samples", required=False, default=-1, type=int, help="Run only this many samples (default: all)"
    )
    parser.add_argument(
        "--use_cached_mc_data",
        required=False,
        default=False,
        type=bool,
        help="Use cached monte carlo data (avoids running random permutation again)",
    )
    parser.add_argument(
        "--num_cores", required=False, default=1, type=int, help="number of cores to use for parallel computing"
    )
    parser.add_argument(
        "--num_permutations", required=False, default=1000, type=int, help="Number of random permutations to run"
    )
    parser.add_argument(
        "--max_cell_distance",
        required=False,
        default=200,
        type=int,
        help="Maximum distance a cell pair can have to each other be called neighbor",
    )
    parser.add_argument(
        "--discrete_cm",
        required=False,
        default=False,
        type=bool,
        help="Use discrete color map with 15 categories in heatmap plots",
    )
    parser.add_argument(
        "--v_min",
        required=False,
        default=-130,
        type=int,
        help="Minimum z-score value to anchor the colormap (default -130)",
    )
    parser.add_argument(
        "--v_max",
        required=False,
        default=130,
        type=int,
        help="Maximum z-score value to anchor the colormap (default 130)",
    )
    parser.add_argument(
        "--num_categories",
        required=False,
        default=31,
        type=int,
        help="Number of categories for discrete colormap (default 31)",
    )
    parser.add_argument(
        "--legend",
        required=False,
        default=False,
        type=bool,
        help="Add phenotype legend (default false)",
    )
    parser.add_argument(
        "--run_sample",
        required=False,
        default="",
        type=str,
        help="Run specified sample (default not set)",
    )
    parser.add_argument(
        "--skip_phenotype",
        required=False,
        default="",
        type=str,
        nargs="*",
        help="Do not plot listed phenotypes",
    )
    parser.add_argument(
        "--tissue_type",
        required=False,
        default="",
        type=str,
        help="Analyze only cells in the given annotated tissue",
    )

    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)

    args = parser.parse_args()

    # Input output settings
    result_dir = args.result_dir
    data_dir = args.data_dir

    # stop after running this many samples (-1 runs all)
    nr_of_samples_to_run = args.run_samples

    # use cached random permutation data
    use_cached_mc_data = args.use_cached_mc_data

    # number of cores to use in parallel functions
    num_cores = args.num_cores

    # number of random permutations
    num_permutations = args.num_permutations

    # maximum distance a cell pair can have to be called neighbor
    maxVoronoiDist = args.max_cell_distance

    # use discrete cm
    discrete_cm = args.discrete_cm
    num_categories = args.num_categories

    # Min/Max z-score values to anchor the colormap
    vmin = args.v_min
    vmax = args.v_max

    # legend
    add_legend = args.legend

    # run only specified sample
    run_sample = args.run_sample

    # list of phenotypes to skip in plots
    skip_phenotype = args.skip_phenotype
    skipped_phenotypes = False

    # tissue type (e.g. TU for tumor)
    tissue_type = args.tissue_type

    #
    # Start processing data files
    #

    # image/sample counter
    i_counter = 0

    for s in sample_map.keys():

        if run_sample is not "" and sample_map[s] != run_sample:
            continue

        # stop after nr_of_samples_to_run (if requested)
        if i_counter == nr_of_samples_to_run and nr_of_samples_to_run > 0:
            print("Done " + str(i_counter) + " images.")
            break

        # reset count tables
        sample_voronoi_interaction_count_table = []
        sample_rnd_voronoi_interaction_count_table = []

        vor_z_score_table = np.zeros((len(ctdef.cellTypeMap), len(ctdef.cellTypeMap)), dtype="float64")
        vor_pr_value_table = np.zeros((len(ctdef.cellTypeMap), len(ctdef.cellTypeMap)), dtype="float64")
        vor_pl_value_table = np.zeros((len(ctdef.cellTypeMap), len(ctdef.cellTypeMap)), dtype="float64")

        sample_id = sample_map[s]
        # get images for sample
        sample_data_files = get_sample_data(data_dir, s)

        # iterate over sample ROI images
        for sample_data_file in sample_data_files:
            print(sample_map[s] + " - " + sample_data_file)
            df = pd.read_csv(sample_data_file, na_values="NA", low_memory=False)
            # df = df[df["phenotype"] != "Undefined"]

            # filter tissue type if given
            tissue_type_filename_part = ""
            if tissue_type is not "":
                df = df[df["tissue"] == tissue_type]
                tissue_type_filename_part = "_" + tissue_type

            df.reset_index(inplace=True)

            # get cell positions and calculate Voronoi diagram
            positions = np.transpose(np.array([df["Location_Center_X"].values, df["Location_Center_Y"].values]))
            vor = Voronoi(positions)

            # count cell-cell neighbors/interactions for current ROI image
            voronoi_interaction_count_table = register_voronoi_interactions(vor, df)

            # append counts of ROI to sample counts
            sample_voronoi_interaction_count_table.append(voronoi_interaction_count_table)

            # if we do not use cached monte carlo data from previous runs, calculate random permutations
            if use_cached_mc_data is False:
                nCells = len(vor.points)

                # set up process pool for parallel processing and run
                pool = Pool(num_cores)
                pool_result = pool.map(register_rnd_voronoi_interactions, [p for p in range(num_permutations)])
                pool.close()

                # store monte carlo results for current image ROI
                rnd_voronoi_interaction_count_table = np.asarray(pool_result)

                # append monte carlo results of current image ROI to sample monte carlo results
                sample_rnd_voronoi_interaction_count_table.append(rnd_voronoi_interaction_count_table)
        # done with all ROI images for current sample

        # if we do not use cached monte carlo data from previous runs, save as cache (pkl) file
        # else
        # read cached (pkl) monte carlo data
        if use_cached_mc_data is False:
            with open(
                result_dir + "/" + sample_id + tissue_type_filename_part + "_sample_rnd_vor_count_table.pkl", mode="wb"
            ) as pkl_f:
                pickle.dump(sample_rnd_voronoi_interaction_count_table, pkl_f)
        else:
            with open(
                result_dir + "/" + sample_id + tissue_type_filename_part + "_sample_rnd_vor_count_table.pkl", "rb"
            ) as pkl_f:
                sample_rnd_voronoi_interaction_count_table = pickle.load(pkl_f)

        # merge list of arrays (all sample ROI interaction counts) to n-dim array
        sample_voronoi_interaction_count_table = np.array(sample_voronoi_interaction_count_table)

        # sum up all ROI image arrays
        sample_voronoi_interaction_count_table = np.nansum(sample_voronoi_interaction_count_table, axis=0)

        tumor_immuneInteractions = np.sum(
            [sample_voronoi_interaction_count_table[x, immuneCellIdx] for x in tumorCellIdx]
        ) + np.sum([sample_voronoi_interaction_count_table[y, tumorCellIdx] for y in immuneCellIdx])
        immune_all_Interactions = np.sum(sample_voronoi_interaction_count_table[immuneCellIdx])

        mixingscore = tumor_immuneInteractions / immune_all_Interactions

        with open(result_dir + "/" + sample_id + tissue_type_filename_part + "_mixing_score.txt", "w") as f:
            f.write("{:s}\t{:.4f}\n".format("Mixingscore:", mixingscore))

        # merge list of monte carlo result arrays to n-dim array [#ROI images, #cell types, #cell types, #permutations]
        sample_rnd_voronoi_interaction_count_table = np.array(sample_rnd_voronoi_interaction_count_table)

        # concatenate single image permutation results for the sample, [#cell types, #cell types, (#images * #permutations)]
        sample_rnd_voronoi_interaction_count_table = np.nansum(sample_rnd_voronoi_interaction_count_table, axis=0)
        sample_rnd_voronoi_interaction_count_table = np.transpose(sample_rnd_voronoi_interaction_count_table, (1, 2, 0))

        # Calculate z-score to determine if pairs are spatially enriched, meaning if there is a
        # significant cell type interaction (positive z-score) or separation (negative z-score)
        for i in range(len(ctdef.cellTypeMap.keys())):
            for j in range(len(ctdef.cellTypeMap.keys())):
                z, p, pr, pl = calculate_zscore_stat(
                    sample_voronoi_interaction_count_table[i, j],
                    sample_rnd_voronoi_interaction_count_table[i, j],
                    num_permutations,
                )
                vor_z_score_table[i, j] = z
                vor_pr_value_table[i, j] = pr
                vor_pl_value_table[i, j] = pl

        # save z-score and p-value table as tsv
        z_score_df = pd.DataFrame(
            data=vor_z_score_table, columns=ctdef.cellTypeNameMap.values(), index=ctdef.cellTypeNameMap.values()
        )
        z_score_df.to_csv(
            result_dir + "/" + sample_id + tissue_type_filename_part + "_voronoi_neighbors_zscores.tsv",
            index=True,
            sep="\t",
            float_format="%.3f",
        )

        p_value_df = pd.DataFrame(
            data=vor_pr_value_table, columns=ctdef.cellTypeNameMap.values(), index=ctdef.cellTypeNameMap.values()
        )
        p_value_df.to_csv(
            result_dir + "/" + sample_id + tissue_type_filename_part + "_voronoi_neighbors_pvalues_rightSided.tsv",
            index=True,
            sep="\t",
            float_format="%.3f",
        )

        p_value_df = pd.DataFrame(
            data=vor_pl_value_table, columns=ctdef.cellTypeNameMap.values(), index=ctdef.cellTypeNameMap.values()
        )
        p_value_df.to_csv(
            result_dir + "/" + sample_id + tissue_type_filename_part + "_voronoi_neighbors_pvalues_leftSided.tsv",
            index=True,
            sep="\t",
            float_format="%.3f",
        )

        #
        # generate heatmaps
        #

        # make data frame from z-score matrix
        #

        # make mask to filter z-scores for significant enrichment (p < 0.01)
        mask_table = np.zeros((len(ctdef.cellTypeMap), len(ctdef.cellTypeMap)), dtype="float64")

        for i in range(len(ctdef.cellTypeMap.keys())):
            for j in range(len(ctdef.cellTypeMap.keys())):
                if np.abs(vor_pr_value_table[i, j]) < 0.01 or np.abs(vor_pl_value_table[i, j]) < 0.01:
                    mask_table[i, j] = 1

        # copy z-score table and filter with mask
        vor_z_score_table_copy = vor_z_score_table.copy()

        # remove phenotypes we do not want to show in plot
        if len(skip_phenotype) > 0:
            remove_phenotype_idx = []

            for phenotype in skip_phenotype:
                if phenotype in ctdef.cellTypeMap:
                    remove_phenotype_idx.append(ctdef.cellTypeMap[phenotype])

                    # only do this for the first image
                    if skipped_phenotypes is not True:
                        del ctdef.cellTypeIdxColorMap[ctdef.cellTypeNameColorMapIdx[ctdef.cellTypeNameMap[phenotype]]]
                        del ctdef.cellTypeColorMap[phenotype]
                        ctdef.cellTypeMapOrder_idx.remove(ctdef.cellTypeNameColorMapIdx[ctdef.cellTypeNameMap[phenotype]])

            # done skipping phenotypes
            skipped_phenotypes = True

            mask_table = np.delete(mask_table, remove_phenotype_idx, axis=0)
            mask_table = np.delete(mask_table, remove_phenotype_idx, axis=1)
            vor_z_score_table_copy = np.delete(vor_z_score_table_copy, remove_phenotype_idx, axis=0)
            vor_z_score_table_copy = np.delete(vor_z_score_table_copy, remove_phenotype_idx, axis=1)

        # and filter with mask
        vor_z_score_table_copy *= mask_table

        # deal with nan values for plotting
        # set nan to 0
        vor_z_score_table_copy[np.isnan(vor_z_score_table_copy)] = 0

        # deal with inf values for plotting
        # set inf first to nan and then to max(z-score) + 0.1
        vor_z_score_table_copy[np.isposinf(vor_z_score_table_copy)] = np.nan
        vor_z_score_table_copy[np.isnan(vor_z_score_table_copy)] = np.nanmax(vor_z_score_table_copy) + 0.1
        vor_z_score_table_copy[np.isneginf(vor_z_score_table_copy)] = np.nan
        vor_z_score_table_copy[np.isnan(vor_z_score_table_copy)] = np.nanmin(vor_z_score_table_copy) - 0.1

        # assign cell type names to rows and columns as used in manuscript (commented, we use numbers)
        # vor_z_score_df = pd.DataFrame(vor_z_score_table_copy, columns=ctdef.cellTypeNameMap.values(), index=ctdef.cellTypeNameMap.values())

        # assign cell type number to rows and columns as used in manuscript
        vor_z_score_df = pd.DataFrame(
            vor_z_score_table_copy, columns=ctdef.cellTypeIdxColorMap.keys(), index=ctdef.cellTypeIdxColorMap.keys()
        )

        # add legend
        if add_legend:
            # set cell type names and colors for legend
            patch_list = []
            for key, val in ctdef.cellTypeColorMap.items():
                data_key = mpatches.Patch(
                    color=val,
                    label=str(
                        ctdef.cellTypeNameColorMapIdx[ctdef.cellTypeNameMap[key]] + ": " + ctdef.cellTypeNameMap[key]
                    ),
                )
                patch_list.append(data_key)

        #
        # plot the clustered heatmap using seaborn clustermap: euclidean distance and
        #

        # color map

        # full
        cmap = "seismic"

        # discretized with 15 catergories
        if discrete_cm is True:
            cmap = sns.color_palette("seismic", num_categories)

        fig = sns.clustermap(
            vor_z_score_df,
            method="weighted",
            cmap=cmap,
            center=0.0,
            vmin=vmin,
            vmax=vmax,
            cbar_pos=(1.0, 0.15, 0.03, 0.5),
            cbar_kws={"shrink": 0.7},
            yticklabels=True,
            xticklabels=True,
            row_colors=ctdef.cellTypeIdxColorMap_s,
            col_colors=ctdef.cellTypeIdxColorMap_s,
            linewidths=0.5,
            figsize=(12, 12),
        )

        fig.ax_heatmap.tick_params(axis="both", which="both", length=0, pad=28)
        fig.ax_heatmap.yaxis.set_ticks_position("left")
        fig.ax_heatmap.xaxis.set_ticks_position("top")

        fig.ax_col_colors.set_yticks([])
        fig.ax_row_colors.set_xticks([])

        fig.ax_row_dendrogram.set_visible(False)
        fig.ax_col_dendrogram.set_visible(False)

        if tissue_type is not "":
            plot_title = sample_id + "(" + tissue_type + ")" + ": direct neighbors permutation z-scores"
        else:
            plot_title = sample_id + ": direct neighbors permutation z-scores"

        fig.ax_heatmap.set_title(plot_title, pad=20)

        # add legend
        if add_legend:
            fig.ax_heatmap.legend(
                handles=patch_list,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.05),
                fontsize="xx-small",
                ncol=5,
                labelspacing=0.25,
                handletextpad=0.5,
            )

        plt.tight_layout()

        fig.savefig(
            result_dir
            + "/"
            + sample_id
            + tissue_type_filename_part
            + "_voronoi_neighbors_frequency_zscores_clustered.png",
            dpi=300,
        )
        fig.savefig(
            result_dir
            + "/"
            + sample_id
            + tissue_type_filename_part
            + "_voronoi_neighbors_frequency_zscores_clustered.pdf",
            dpi=300,
        )

        #
        # plot the custom ordered heatmap using seaborn
        #

        # order z-score table by custom order
        custom_sorted_vor_z_score_df = vor_z_score_df[ctdef.cellTypeMapOrder_idx]
        custom_sorted_vor_z_score_df = custom_sorted_vor_z_score_df.reindex(ctdef.cellTypeMapOrder_idx)

        # color map

        # full
        cmap = "seismic"

        # discretized with 15 catergories
        if discrete_cm is True:
            cmap = sns.color_palette("seismic", num_categories)

        fig = sns.clustermap(
            custom_sorted_vor_z_score_df,
            row_cluster=False,
            col_cluster=False,
            cmap=cmap,
            center=0.0,
            vmin=vmin,  # -(np.max(pair_z_score_table_copy))*0.75,
            vmax=vmax,  # np.max(pair_z_score_table_copy)*0.75,
            cbar_pos=(1.0, 0.15, 0.03, 0.5),
            cbar_kws={"shrink": 0.7},
            yticklabels=True,
            xticklabels=True,
            row_colors=ctdef.cellTypeIdxColorMap_s,
            col_colors=ctdef.cellTypeIdxColorMap_s,
            linewidths=0.5,
            figsize=(12, 12),
        )

        fig.ax_heatmap.tick_params(axis="both", which="both", length=0, pad=28)
        fig.ax_heatmap.yaxis.set_ticks_position("left")
        fig.ax_heatmap.xaxis.set_ticks_position("top")

        fig.ax_col_colors.set_yticks([])
        fig.ax_row_colors.set_xticks([])

        fig.ax_row_dendrogram.set_visible(False)
        fig.ax_col_dendrogram.set_visible(False)

        fig.ax_heatmap.set_title(plot_title, pad=20)

        # add legend
        if add_legend:
            fig.ax_heatmap.legend(
                handles=patch_list,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.05),
                fontsize="xx-small",
                ncol=5,
                labelspacing=0.25,
                handletextpad=0.5,
            )

        plt.tight_layout()

        fig.savefig(
            result_dir
            + "/"
            + sample_id
            + tissue_type_filename_part
            + "_voronoi_neighbors_frequency_zscores_ordered.png",
            dpi=300,
        )
        fig.savefig(
            result_dir
            + "/"
            + sample_id
            + tissue_type_filename_part
            + "_voronoi_neighbors_frequency_zscores_ordered.pdf",
            dpi=300,
        )

        # Done with sample increase image/sample counter
        i_counter += 1
