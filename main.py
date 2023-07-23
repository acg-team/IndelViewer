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
import os, io, random
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
from panel.widgets import Tqdm

class Initializer:

    path = './sample/'
    files_path = {
        'fasta': './sample/prank_phyml_OMAGroup_1003884.fa.best.arpipasr.fasta',
        'relation': './sample/prank_phyml_OMAGroup_1003884.fa.best.arpipnode_rel.txt',
        'tree': './sample/prank_phyml_OMAGroup_1003884.fa.best.arpiptree.nwk',
        'indel': './sample/prank_phyml_OMAGroup_1003884.fa.best.arpipindel.txt'
    }
    def __init__(self, files_path=files_path):

        self.file_input_fasta = pn.widgets.FileInput(accept='.fasta')
        self.file_input_arpipnode_rel_txt = pn.widgets.FileInput(accept='.txt')
        self.file_input_nwk = pn.widgets.FileInput(accept='.nwk')
        self.file_input_arpipindel_txt = pn.widgets.FileInput(accept='.txt')

        self.apps = pn.Column(
            pn.Column(pn.pane.Markdown('## Indel Viewer fasta Input'), self.file_input_fasta),
            pn.Column(pn.pane.Markdown('## Indel Viewer arpipindel_rel Input'), self.file_input_arpipnode_rel_txt),
            pn.Column(pn.pane.Markdown('## Indel Viewer nwk Input'), self.file_input_nwk),
            pn.Column(pn.pane.Markdown('## Indel Viewer arpipindel Input'), self.file_input_arpipindel_txt),
        )




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
            if node not in self.tree.iter_leaves():
                a = node.children[0]
                node.name = self.dic_relation.get(a.name)
            self.list_names.append(node.name)

        events = Helper.read_input_file(files_path['indel'])
        self.df_event, self.mat_insertion, self.mat_deletion, self.dic_name_index = Helper.get_indel_events(
            self.seq_len, events, self.list_names)


# The main panel
selected_site_number = 2546

def main(args=None):
    ts = TreeStyle()
    initializer = Initializer()

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

    # create widget tool to select the coloumn
    site_slider = pn.widgets.IntSlider(name="value", start=1, end=initializer.seq_len, value=1, width=500)
    file_input_fasta = pn.widgets.FileInput(accept='.fasta')
    # bottom
    app = pn.Tabs(
        # ('Import Files', initializer.apps),

        ('Work View',
         pn.Column(
             pn.Row(tree_panel, seq_pn),
             pn.Column(
                 pn.Row(
                     pn.pane.Markdown('## Indel Viewer'),
                 ),
                 pn.layout.Divider(margin=(-20, 0, 0, 0)),
                 pn.Row(
                     pn.Column(
                         pn.Row(
                             "Simple example."),
                         pn.Row(indel_view_btn)
                     ),

                     pn.Column(
                         pn.Row(
                             "Adjust the hyperparameters:"),
                         pn.Row(indel_view_btn_plotting),
                         pn.Row(tqdm)
                     ),
                 ),
                 styles=dict(background='WhiteSmoke', width='320'),
             )
         )
         )
    )

    app.show()  # html page
    # app  # inline


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()


