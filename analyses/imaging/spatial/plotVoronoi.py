#!/usr/bin/env python
# coding: utf-8

"""
Plot Voronoi diagrams of phenotyped IMC images.

Requirements:
    * Python >= 3.7.0
    * scipy >= 1.4.1
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


# Import libraries
#

# In[1]:
# Import libraries
#

import os
import sys
import argparse
import numpy as np
import pandas as pd
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

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

# In[4]:
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


# In[5]:

if __name__ == "__main__":

    # Argument parser
    parser = argparse.ArgumentParser(description="Plot Voronoi diagrams with cell type specific colors")
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
        "--plot_size_x", required=False, default=10, type=float, help="Size of plot in inches horizontal"
    )
    parser.add_argument("--plot_size_y", required=False, default=10, type=float, help="Size of plot in inches vertical")
    parser.add_argument(
        "--image_size_x",
        required=False,
        default=1000,
        type=int,
        help="Size of original ROI image in pixels horizontal",
    )
    parser.add_argument(
        "--image_size_y",
        required=False,
        default=1000,
        type=int,
        help="Size of original ROI image in pixels vertical",
    )
    parser.add_argument(
        "--make_cropped_plot",
        required=False,
        default=False,
        action="store_true",
        help="Make a cropped region of the plot, it set default origin is 0,0 and copped size is 200x200, may be adjusted via --image_crop_origin x,y and image_crop_size x,y",
    )
    parser.add_argument(
        "--image_crop_origin",
        required=False,
        default=[0, 0],
        nargs=2,
        type=int,
        help="x y coordinates of origin for cropped image",
    )
    parser.add_argument(
        "--image_crop_size",
        required=False,
        default=[200, 200],
        nargs=2,
        type=int,
        help="x y size in pixels of cropped image",
    )
    parser.add_argument(
        "--ridge_line_width",
        required=False,
        default=0.25,
        type=float,
        help="line width of Voronoi ridges",
    )
    parser.add_argument(
        "--nucleus_marker_size",
        required=False,
        default=0.5,
        type=float,
        help="size of Voronoi point marker",
    )

    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)

    args = parser.parse_args()

    # Input output settings
    result_dir = args.result_dir
    data_dir = args.data_dir

    # plot in inches
    xmax = args.plot_size_x
    ymax = args.plot_size_y

    # image size in px
    xmaxImg = args.image_size_x
    ymaxImg = args.image_size_y

    # set voronoi ridge line width
    ridge_line_width = args.ridge_line_width

    # set voronoi point size
    nucleus_marker_size = args.nucleus_marker_size

    # lets see if we need to make a cropped plot
    if args.make_cropped_plot:
        image_crop_origin = args.image_crop_origin
        image_crop_size = args.image_crop_size

        if image_crop_origin[0] + image_crop_size[0] >= xmaxImg or image_crop_origin[1] + image_crop_size[1] >= ymaxImg:
            print("Cropped section must fit into ROI image")
            sys.exit(1)

        if image_crop_size[0] != image_crop_size[1]:
            if image_crop_size[0] > image_crop_size[1]:
                scale_plot_y = image_crop_size[1] / image_crop_size[0]
                ymax = ymax * scale_plot_y
            else:
                scale_plot_x = image_crop_size[0] / image_crop_size[1]
                xmax = xmax * scale_plot_x

        ridge_line_width = 1000 / max(image_crop_size) * ridge_line_width
        nucleus_marker_size = 1000 / max(image_crop_size) * nucleus_marker_size

    # stop after running this many samples (-1 runs all)
    nr_of_samples_to_run = args.run_samples

    # image/sample counter
    i_counter = 0
    for s in sample_map.keys():

        # stop after nr_of_samples_to_run (if requested)
        if i_counter == nr_of_samples_to_run and nr_of_samples_to_run > 0:
            print("Done " + str(i_counter) + " images.")
            break

        sample_id = sample_map[s]
        # get images for sample
        sample_data_files = get_sample_data(data_dir, s)

        # iterate over sample ROI images
        for sample_data_file in sample_data_files:
            print(sample_map[s] + " - " + sample_data_file)
            fname = os.path.basename(sample_data_file)
            tile_name = fname.replace("_with_phenotype.csv", "")

            df = pd.read_csv(sample_data_file, na_values="NA", low_memory=False)
            # df = df[df["phenotype"] != "Undefined"]
            # df.reset_index(inplace=True)

            # get cell positions and calculate Voronoi diagram
            positions = np.transpose(np.array([df["Location_Center_X"].values, df["Location_Center_Y"].values]))
            # add far distance points for edge filling in voronoi plot
            positions = np.append(positions, [[9999, 9999], [-9999, 9999], [9999, -9999], [-9999, -9999]], axis=0)

            vor = Voronoi(positions)
            # voronoi_plot_2d(vor, show_vertices=False, line_colors="white", line_width=2.5, line_alpha=0.9, point_size=5)

            # make Voronoi diagram figure
            plt.figure(figsize=(xmax, ymax))

            # define cell nucleus center of mass marker
            marker_style = dict(color="tab:blue", marker="o", markersize=nucleus_marker_size, fillstyle="full")

            # reset
            cell_color = {}
            cell_nr = 0
            cell_types_present = []

            # loop over point_regions and draw corresponding polygons from vertices
            # and fill with cell type specific color
            for p_reg in vor.point_region:
                if not -1 in vor.regions[p_reg]:
                    cell_type = df["phenotype"][cell_nr]
                    polygon = [vor.vertices[i] for i in vor.regions[p_reg]]
                    plt.fill(
                        *zip(*polygon), ctdef.cellTypeColorMap[cell_type], edgecolor="w", linewidth=ridge_line_width
                    )
                    cell_types_present.append(cell_type)
                cell_nr += 1

            # set cell type names and colors for legend
            patch_list = []
            for key, val in ctdef.cellTypeColorMap.items():
                if key in cell_types_present:
                    data_key = mpatches.Patch(color=val, label=ctdef.cellTypeNameMap[key])
                    patch_list.append(data_key)

            plt.title(sample_id + ": Voronoi diagram")  # title with fontsize 20
            plt.plot(positions[:-1, 0], positions[:-1, 1], "o", **marker_style)
            plt.xlim([0, xmaxImg]), plt.ylim([0, ymaxImg])

            if args.make_cropped_plot:
                # plot only cropped section
                plt.xlim([image_crop_origin[0], image_crop_origin[0] + image_crop_size[0]])
                plt.ylim([image_crop_origin[1], image_crop_origin[1] + image_crop_size[1]])

            # add legend
            plt.legend(
                handles=patch_list,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.05),
                fontsize="xx-small",
                ncol=5,
                labelspacing=0.25,
                handletextpad=0.5,
            )

            plt.tight_layout()

            if args.make_cropped_plot:
                plt.savefig(
                    result_dir
                    + "/"
                    + sample_id
                    + "_"
                    + tile_name
                    + "_voronoi_plot_cropped_"
                    + str(image_crop_origin[0])
                    + "_"
                    + str(image_crop_origin[1])
                    + "_"
                    + str(image_crop_size[0])
                    + "x"
                    + str(image_crop_size[1])
                    + ".png",
                    dpi=300,
                )
            else:
                plt.savefig(result_dir + "/" + sample_id + "_" + tile_name + "_voronoi_plot.png", dpi=300)
            plt.close()

        # Done with sample increase image/sample counter
        i_counter += 1
