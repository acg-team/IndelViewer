########################################################################################################################
#
# Licensed Materials - Property of Gholamhossein Jowkar
# Copyright (C) 2022-2025 by Gholamhossein Jowkar
#
########################################################################################################################
#
# This file is part of arpip_indel_viewer of ARPIP project: https://github.com/acg-team/IndelViewer

#                           ABOUT THE ARPIP PACKAGE
#                           =====================
# ARPIP: Ancestral Sequence Reconstruction with insertions and deletions under the Poisson Indel Process
# ARPIP is a joint maximum likelihood approach for phylogenetic ancestral sequence reconstruction, capable of modeling
# indels both biological and mathematically.
#

# arpip_indel_viewer is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# arpip_indel_viewer is a free software: you can redistribute it and/or modify it under the terms of the GNU Lesser
# General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# You should have received a copy of the GNU Lesser General Public
# License along with IndelViewer. If not, see <http://www.gnu.org/licenses/>.
#
########################################################################################################################

# Python 3
# main.py
# Created by: Gholam-Hossein Jowkar
# ACGT ZHAW
# Created date: September 2022
# Modified by:
# Modified date:

import logging

from viewer.helper import Helper
from viewer.indel_tree import DrawTree
from viewer.plotting_events import PlottingEvents
from viewer.helper import *
# from tree import DrawTree
import os, io, random, datetime, argparse, pathlib
import string
import numpy as np
import pandas as pd

# from Biopython
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO, Seq

# from panel
import panel as pn
import panel.widgets as pnw
from panel.widgets import Tqdm
pn.extension()

# from Bokeh
from bokeh import events
from bokeh.plotting import figure, show
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d, CustomJS, Slider, LabelSet
from bokeh.models import RangeTool, Button, CustomJS, Div, TapTool
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot, column, row

# form ETE3
from ete3 import Tree, PhyloTree, TreeStyle, NodeStyle, SeqGroup


class Initializer:
    parser = argparse.ArgumentParser(description='Indel Viewer')
    parser.add_argument('input_file', type=str, help='Path to the input file')
    parser.add_argument('--verbose', '-v', action='store_true', default=0, help='Verbosity level')
    # parser.add_argument('path', type=pathlib.Path)

    parsed_args = parser.parse_args()

    verbose = parsed_args.verbose

    with open(parsed_args.input_file) as file:
        dict_input_file = {}
        for line in file:
            key, value = line.rstrip().split('=', 1)
            dict_input_file[key] = value

        necessary_keys = ['fasta', 'relation', 'tree', 'indel']
        input_file_keys = list(dict_input_file.keys())
        # logging.info("The existing file extensions:%s", dict_input_file)
        # for key, value in dict_input_file.items():
        #     logging.INFO(f"{key}: {value}")

        if len(necessary_keys) == len(input_file_keys):
            if (dict_input_file['fasta'].split('.', -1)[-1] == 'fasta' and
                    dict_input_file['relation'].split('.', -1)[-1] == 'txt' and
                    dict_input_file['tree'].split('.', -1)[-1] == 'nwk' and
                    dict_input_file['indel'].split('.', -1)[-1] == 'txt'):
                files_path = dict_input_file
            else:
                logging.error("[Main:Initializer] The input text file is incorrect")
        else:
            logging.error("[Main:Initializer] There should be four files with different formats.")
            print("Input text file is incorrect, the sample file is running instead")
            logging.warning("[Main:Initializer]Input text file is incorrect, the sample file is running instead")
            path = './sample/'
            files_path = {
                'fasta': './sample/prank_phyml_OMAGroup_1003884.fa.best.arpipasr.fasta',
                'relation': './sample/prank_phyml_OMAGroup_1003884.fa.best.arpipnode_rel.txt',
                'tree': './sample/prank_phyml_OMAGroup_1003884.fa.best.arpiptree.nwk',
                'indel': './sample/prank_phyml_OMAGroup_1003884.fa.best.arpipindel.txt'
            }

        if verbose:
            print("Verbose mode is enabled")

    def __init__(self, files_path=files_path):

        self.aligned_data = AlignIO.read(f"{files_path['fasta']}", 'fasta')
        self.seq_len = self.aligned_data.get_alignment_length()

        lines = Helper.read_input_file(files_path['relation'])
        self.dic_relation = {}
        for line in lines[1:-1]:
            s = line.split('\t')
            self.dic_relation[s[0]] = s[1].strip()

        self.tree = Tree(f"{files_path['tree']}", format=1)

        self.list_names = []
        for node in self.tree.traverse("postorder"):
            # if node not in self.tree.iter_leaves():
            #     a = node.children[0]
            #     node.name = self.dic_relation.get(a.name)
            self.list_names.append(node.name)

        events = Helper.read_input_file(files_path['indel'])
        self.df_event, self.mat_insertion, self.mat_deletion, self.dic_name_index = Helper.get_indel_events(
            self.seq_len, events, self.list_names)


# The main panel
selected_site_number = 0

def main(args=None):

    # Configure logging settings
    logging.basicConfig(filename='indel_view.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    # Configure logging settings for both console and file handlers
    # logger = logging.getLogger('')
    # logger.setLevel(logging.DEBUG)

    ts = TreeStyle()
    initializer = Initializer()
    logging.info("[Main] The input files are read successfully.")

    DrawTree(ts, initializer, selected_site_number)
    tree_panel = pn.pane.image.PNG(initializer.tree.render("%%inline", tree_style=ts), height=300)
    tree_result = pn.pane.image.PNG(name='empty', height=600)

    def update_graph(event):
        global tree_panel, tree_result, ts
        ts = TreeStyle()
        DrawTree(ts, initializer, selected_site_number)
        tree_panel = pn.pane.image.PNG(initializer.tree.render("%%inline", tree_style=ts), height=300)
        tree_result = pn.pane.image.PNG(name='empty', height=600)

        app[0][0][0] = tree_panel

    def update_selected_site_number(number):
        global selected_site_number
        selected_site_number = int(number) - 1

    indel_view_btn = pn.widgets.Button(name='show the scenario', width=200, button_type='primary')
    # watch associates a button click with the update function
    indel_view_btn.param.watch(update_graph, 'clicks')

    tqdm = Tqdm(width=280)

    def generate_plotting_events(event):
        plottingTree = PhyloTree('./sample/prank_phyml_OMAGroup_1003884.fa.best.arpiptree.nwk', alg_format='fasta',
                                 format=1)
        plottingEvents = PlottingEvents(plottingTree, initializer, tqdm)
        plottingEvents.print_all_sites()

    indel_view_btn_plotting = pn.widgets.Button(name='Plotting All Events', width=200, button_type='primary')

    # watch associates a button click with the update function
    indel_view_btn_plotting.param.watch(generate_plotting_events, 'clicks')

    # viewing the alignment
    seq_pn = pn.pane.Bokeh(name='align', height=300)
    seq_pn.object = Helper.view_alignment(initializer.aligned_data, initializer.df_event, update_selected_site_number,
                                          molecule_type='AA', plot_width=600)

    # create widget tool to select the column
    site_slider = pn.widgets.IntSlider(name="value", start=1, end=initializer.seq_len, value=1, width=500)

    # bottom
    app = pn.Tabs(
        # ('Import Files', initializer.apps),

        ('Result View',
             pn.Column(
                 pn.Row(tree_panel, seq_pn),
                 pn.Column(
                     pn.Row(
                         pn.pane.Markdown('## ARPIP Indel Viewer'),
                     ),
                     pn.layout.Divider(margin=(-20, 0, 0, 0)),
                     pn.Row(
                         pn.Column(
                             pn.Row(
                                 "Show the scenario:"),
                             pn.Row(indel_view_btn)
                         ),

                         pn.Column(
                             pn.Row(
                                 "Print all the scenarios:"),
                             pn.Row(indel_view_btn_plotting),
                             pn.Row(tqdm)
                         ),
                     ),
                     # styles=dict(background='WhiteSmoke', width='320'),
                 )
             )
         )
    )

    app.show(title="Indel Viewer")  # html page
    # app  # inline


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()


