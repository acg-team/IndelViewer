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
# helper.py
# Created by: Gholam-hossein Jowkar
# ACGT ZHAW
# Created date: November 2022
# Modified by:
# Modified date:


import logging
import math
import os, io, random
import string
import numpy as np
import pandas as pd

# from Biopython
# from Bio.Seq import Seq
# from Bio.Align import MultipleSeqAlignment
# from Bio import AlignIO, SeqIO, Seq

import panel as pn
import panel.widgets as pnw

pn.extension()

# from Bokeh
from bokeh import events
from bokeh.plotting import figure, show
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d, Slider, LabelSet
from bokeh.models import RangeTool, Button, CustomJS, Div, TapTool, TextInput
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot, column, row
from bokeh.layouts import layout
from bokeh.events import DoubleTap, ButtonClick

# form ETE3
from ete3 import Tree, PhyloTree, TreeStyle, NodeStyle, SeqGroup


class Helper():

    @staticmethod
    def read_file(path, file):
        """
           Reading all line of a file
        """
        file_path = f"{path}/{file}"
        f = open(file_path, 'r')
        return f.readlines()

    @staticmethod
    def read_input_file(path):
        f = open(path, 'r')
        return f.readlines()

    @staticmethod
    def get_indel_events(seq_len, events, list_names):
        """
        param: events: event file in text format
        param: list_names: name of all nodes
        return: putting the event in 3 files: mat_in/del matrix of events separately
        mat_insertion/deletion of size(m*n) and events of size(2*m)
        """
        dic_event = {}
        mat_insertion = np.zeros([seq_len, len(list_names)])
        # mat_insertion
        mat_deletion = np.zeros([seq_len, len(list_names)])
        # mat_deletion
        dic_name_index = {}
        counter = 0
        # making a dictionary of name and index of tree's nodes
        for element in list_names:
            dic_name_index[element] = counter
            counter += 1
        for i, line in enumerate(events):
            try:
                if line.strip():
                    dic_event[i + 1] = line.strip()
                    s = line.split(';')
                    for sub in s:
                        if ':X' in sub:
                            strd = sub.split(':X')
                            for sub_str_id in strd[:-1]:
                                mat_deletion[i, dic_name_index[sub_str_id]] = 1
                        elif ':I' in sub:
                            stri = sub.split(':I')
                            for sub_str_id in stri[:-1]:
                                mat_insertion[i, dic_name_index[sub_str_id]] = 1
            except KeyError:
                if not (sub_str_id in list_names):
                    print("Error: Node with name [%s]  could not be found in the tree." % sub_str_id)
                else:
                    print("ERROR")

        df_event = pd.Series(dic_event)
        # df_event
        return df_event, mat_insertion, mat_deletion, dic_name_index

        @staticmethod
        def get_color(seqs, alphabet='DNA'):
            """
            :param seqs: A single MultipleSeqAlignment object.
            :param alphabet: The alphabet could be DNA or AA.
            :return: Dictionary of colors corresponding to the alphabet.
            Usage:
            >> Make colors for bases in sequence.
            """
            text = [i for s in list(seqs) for i in s]
            if alphabet == 'DNA':
                clrs = {'A': 'red', 'T': 'green', 'G': 'orange', 'C': 'blue', '-': 'white'}
            elif alphabet == 'AA':
                clrs = {'-': 'lightgray', 'A': 'royalblue', 'R': 'red', 'N': 'green',
                        'D': 'orchid', 'C': 'sandybrown', 'E': 'darkorchid', 'Q': 'lime',
                        'Z': 'snow', 'G': 'darkorange', 'H': 'turquoise', 'I': 'deepskyblue',
                        'L': 'lightskyblue', 'K': 'orangered', 'M': 'midnightblue', 'F': 'cornflowerblue',
                        'P': 'gold', 'S': 'darkgreen', 'T': 'limegreen', 'W': 'navy',
                        'Y': 'aqua', 'V': 'blue', 'B': 'snow'}
            colors = [clrs[i] for i in text]
            return colors

        @staticmethod
        def display_event(div, attributes=[]):
            """
            Function to build a suitable CustomJS to display the current event in the div model.
            """
            style = 'float: left; clear: left; font-size: 13px'
            return CustomJS(args=dict(div=div), code="""
                const attrs = %s;
                const args = [];
                for (let i = 0; i < attrs.length; i++) {
                    const val = JSON.stringify(cb_obj[attrs[i]], function(key, val) {
                        return val.toFixed ? Number(val.toFixed(2)) : val;
                    })
                    args.push(attrs[i] + '=' + val)
                }
                const line = "<span style=%r><b>" + cb_obj.event_name + "</b>(" + args.join(", ") + ")</span>\\n";
                const text = div.text.concat(line);
                const lines = text.split("\\n")
                if (lines.length > 35)
                    lines.shift();
                div.text = lines.join("\\n");
            """ % (attributes, style))

        @staticmethod
        def view_alignment(aln, indel_events, func, molecule_type='DNA', fontsize="5pt", plot_width=800):
            """
            :param aln: The alignment file using BioPython AlignIO.read()
            :param fontsize: The font size of the alphabet.
            :param molecule_type: The alphabet is 'DNA' or 'AA'.
            Usage:
            >> Bokeh sequence alignment view from https://dmnfarrell.github.io/
            """
            # make sequence and id list from the aln object
            seqs = [rec.seq for rec in (aln)]
            ids = [rec.id for rec in aln]
            text = [i for s in list(seqs) for i in s]
            colors = Helper.get_color(seqs, molecule_type)

            len_seq = len(seqs[0])
            nb_seq = len(seqs)
            width = 0.4

            x = np.arange(1, len_seq + 1)
            y = np.arange(0, nb_seq, 1)

            # creates a 2D grid of coords from the 1D arrays
            xx, yy = np.meshgrid(x, y)

            indel_events_nd_projection = np.tile(indel_events, (nb_seq, 1))  # n copy of indel event array

            # flattens the arrays
            gx = xx.ravel()
            gy = yy.flatten()
            g_indel = indel_events_nd_projection.flatten()

            # use recty for rect coords with an offset
            recty = gy + 0.6
            h = 1 / nb_seq

            # create the ColumnDataSource with all the arrays
            source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors, indel=g_indel))
            plot_height = nb_seq * 9 + 50
            x_range = Range1d(0, len_seq + 1, bounds='auto')

            if len_seq > 70:
                viewlen = 70
            else:
                viewlen = len_seq

            # view_range is for the close up view
            view_range = (0, viewlen)
            tools = "xpan, hover, doubletap, xwheel_zoom, reset, save"

            # select_tools = ['tap']

            # entire sequence view (no text, with zoom)
            p1 = figure(title=None, width=plot_width, height=50,
                        x_range=x_range, y_range=(0, nb_seq), tools="", toolbar_location=None,
                        min_border=0, background_fill_color="#efefef")  # toolbar_location='below'
            rects = Rect(x="x", y="recty", width=1, height=1, fill_color="colors",
                         line_color=None, fill_alpha=0.65)
            p1.add_glyph(source, rects)
            p1.yaxis.visible = False
            p1.grid.visible = False

            # sequence text view with ability to scroll along x axis
            p2 = figure(title=None, width=plot_width, height=plot_height,
                        x_range=view_range, y_range=ids, tools=tools,
                        min_border=0, toolbar_location='below')  # ,tools="xpan,reset", lod_factor=1)
            glyph = Text(x="x", y="y", text="text", text_align='center', text_baseline='bottom',
                         text_color="black", text_font="monospace", text_font_size=fontsize)
            rects = Rect(x="x", y="recty", width=1, height=1, fill_color="colors",
                         line_color=None, fill_alpha=0.7)
            p2.add_glyph(source, glyph)
            p2.add_glyph(source, rects)

            p2.grid.visible = False
            p2.xaxis.major_label_text_font_style = "bold"
            p2.yaxis.minor_tick_line_width = 0
            p2.yaxis.major_tick_line_width = 0

            # hover action
            p2.hover.tooltips = [
                ("site number", "@x"),
                ("event", "@indel")
            ]

            #     p2.js_on_event(events.Tap, display_event(p2, attributes=point_attributes))

            # connection between two figures
            range_tool = RangeTool(x_range=p2.x_range)
            range_tool.overlay.fill_color = "gold"
            range_tool.overlay.fill_alpha = 0.4
            #     range_tool.overlay.

            p1.ygrid.grid_line_color = None
            p1.add_tools(range_tool)
            p2.toolbar.active_multi = 'auto'
            p1.toolbar.active_multi = 'auto'

            # ------------- callback ----------------
            def callback(event):
                selected = source.selected.indices[0]
                selected_row_index = math.floor(selected / len_seq)

                result_t = []
                for k in range(0, nb_seq):
                    if selected_row_index < k:
                        result_t += [selected + (len_seq * (k - selected_row_index))]
                    elif selected_row_index > k:
                        result_t += [selected - (len_seq * (selected_row_index - k))]
                    result_t += [selected]

                source.selected.indices = result_t

            taptool = p2.select(type=TapTool)
            p2.on_event(DoubleTap, callback)

            # source.selected.on_event("indices", my_tap_handler)

            # events.SelectionGeometry
            # events.Tap,
            # events.DoubleTap

            # source.js_on_event(events.DoubleTap, CustomJS(args=dict(s=source), code="""
            # selected = s.selected.indices[0]
            # selected_row_index = math.floor(selected/len_seq)

            # result_t = []
            # for k in range(0,nb_seq):
            # 	if selected_row_index < k:
            # 		result_t += [selected + (len_seq * (k-selected_row_index))]
            # 	elif selected_row_index > k:
            # 		result_t += [selected - (len_seq * (selected_row_index-k))]
            # 	result_t += [selected]

            # s.selected.indices = result_t

            # """))

            code = """
            var selected = source.selected.indices[0],
            selected_row_index = Math.floor(selected/len_seq),
            result_t = [];
        
            for (let k = 0; k < nb_seq; k++) {
                if(selected_row_index < k)
                    result_t.push(selected + (len_seq * (k-selected_row_index)));
                else if(selected_row_index > k)
                    result_t.push(selected - (len_seq * (selected_row_index-k)));
                result_t.push(selected);
            }

            source.selected.indices = result_t;
            source.change.emit();
            """

            callback = CustomJS(args={'source': source, 'len_seq': len_seq, 'nb_seq': nb_seq}, code=code)
            # def call_back_column_number_tap(event):
            # 	func(text_input.value)
            # 	return callback

            p2.add_tools(TapTool(gesture='doubletap', callback=callback))

            text_input = TextInput(value="", max_width=250, title="Enter the sequence column number:")
            # text_input.js_on_change("value", CustomJS(code="""
            #     console.log('text_input: value=' + this.value, this.toString())
            # """))

            button = Button(label="Select Column", align="end", width=200, max_width=80, button_type="primary")
            button.js_on_click(CustomJS(
                args=dict(source=source, text_input=text_input, len_seq=len_seq, nb_seq=nb_seq, range_tool=range_tool),
                code="""
                var result_t = [],
                    selected_column = parseInt(text_input.value);
                for (let k = 0; k < nb_seq; k++) {
                    if(0 < k)
                        result_t.push(selected_column + (len_seq * k));
                    result_t.push(selected_column);
                }
                console.log((selected_column-35,selected_column+34))
                if(selected_column>36){
                    range_tool.x_range.start = selected_column-35;
                    range_tool.x_range.end = selected_column+34;
                }
            
                source.selected.indices = result_t;
            """))

            def call_back_column_number(event):
                func(text_input.value)

            def call_back_column_number_tap(event):
                func(source.selected.indices[0])

            # button.on_click(func(15))
            button.on_event(ButtonClick, call_back_column_number)

            p2.on_event(DoubleTap, call_back_column_number_tap)

            selecte_column = row(text_input, button, max_width=plot_width)  # , sizing_mode="stretch_width"

            p = gridplot([[p1], [p2], [selecte_column]], toolbar_location='right')

            return p

    @staticmethod
    def get_color(seqs, alphabet='DNA'):
        """
        :param seqs: A single MultipleSeqAlignment object.
        :param alphabet: The alphabet could be DNA or AA.
        :return: Dictionary of colors corresponding to the alphabet.
        Usage:
        >> Make colors for bases in sequence.
        """
        text = [i for s in list(seqs) for i in s]
        if alphabet == 'DNA':
            clrs = {'A': 'red', 'T': 'green', 'G': 'orange', 'C': 'blue', '-': 'white'}
        elif alphabet == 'AA':
            clrs = {'-': 'lightgray', 'A': 'royalblue', 'R': 'red', 'N': 'green',
                    'D': 'orchid', 'C': 'sandybrown', 'E': 'darkorchid', 'Q': 'lime',
                    'Z': 'snow', 'G': 'darkorange', 'H': 'turquoise', 'I': 'deepskyblue',
                    'L': 'lightskyblue', 'K': 'orangered', 'M': 'midnightblue', 'F': 'cornflowerblue',
                    'P': 'gold', 'S': 'darkgreen', 'T': 'limegreen', 'W': 'navy',
                    'Y': 'aqua', 'V': 'blue', 'B': 'snow'}
        colors = [clrs[i] for i in text]
        return colors

    @staticmethod
    def display_event(div, attributes=[]):
            """
            Function to build a suitable CustomJS to display the current event
            in the div model.
            """
            style = 'float: left; clear: left; font-size: 13px'
            return CustomJS(args=dict(div=div), code="""
                const attrs = %s;
                const args = [];
                for (let i = 0; i < attrs.length; i++) {
                    const val = JSON.stringify(cb_obj[attrs[i]], function(key, val) {
                        return val.toFixed ? Number(val.toFixed(2)) : val;
                    })
                    args.push(attrs[i] + '=' + val)
                }
                const line = "<span style=%r><b>" + cb_obj.event_name + "</b>(" + args.join(", ") + ")</span>\\n";
                const text = div.text.concat(line);
                const lines = text.split("\\n")
                if (lines.length > 35)
                    lines.shift();
                div.text = lines.join("\\n");
                """ % (attributes, style))

    @staticmethod
    def view_alignment(aln, indel_events, func, molecule_type='DNA', fontsize="5pt", plot_width=800):
        """
        :param aln: The alignment file using BioPython AlignIO.read()
        :param fontsize: The font size of the alphabet.
        :param molecule_type: The alphabet is 'DNA' or 'AA'.
        Usage:
        >> Bokeh sequence alignment view from https://dmnfarrell.github.io/
        """
        # make sequence and id list from the aln object
        seqs = [rec.seq for rec in (aln)]
        ids = [rec.id for rec in aln]
        text = [i for s in list(seqs) for i in s]
        colors = Helper.get_color(seqs, molecule_type)

        len_seq = len(seqs[0])
        nb_seq = len(seqs)
        width = 0.4

        x = np.arange(1, len_seq + 1)
        y = np.arange(0, nb_seq, 1)

        # creates a 2D grid of coords from the 1D arrays
        xx, yy = np.meshgrid(x, y)

        indel_events_nd_projection = np.tile(indel_events, (nb_seq, 1))  # n copy of indel event array

        # flattens the arrays
        gx = xx.ravel()
        gy = yy.flatten()
        g_indel = indel_events_nd_projection.flatten()

        # use recty for rect coords with an offset
        recty = gy + 0.6
        h = 1 / nb_seq

        # create the ColumnDataSource with all the arrays
        source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors, indel=g_indel))
        plot_height = nb_seq * 9 + 50
        x_range = Range1d(0, len_seq + 1, bounds='auto')

        if len_seq > 70:
            viewlen = 70
        else:
            viewlen = len_seq

        # view_range is for the close up view
        view_range = (0, viewlen)
        tools = "xpan, hover, doubletap, xwheel_zoom, reset, save"

        # select_tools = ['tap']

        # entire sequence view (no text, with zoom)
        p1 = figure(title=None, width=plot_width, height=50,
                    x_range=x_range, y_range=(0, nb_seq), tools="", toolbar_location=None,
                    min_border=0, background_fill_color="#efefef")  # toolbar_location='below'
        rects = Rect(x="x", y="recty", width=1, height=1, fill_color="colors",
                     line_color=None, fill_alpha=0.65)
        p1.add_glyph(source, rects)
        p1.yaxis.visible = False
        p1.grid.visible = False

        # sequence text view with ability to scroll along x axis
        p2 = figure(title=None, width=plot_width, height=plot_height,
                    x_range=view_range, y_range=ids, tools=tools,
                    min_border=0, toolbar_location='below')  # ,tools="xpan,reset", lod_factor=1)
        glyph = Text(x="x", y="y", text="text", text_align='center', text_baseline='bottom',
                     text_color="black", text_font="monospace", text_font_size=fontsize)
        rects = Rect(x="x", y="recty", width=1, height=1, fill_color="colors",
                     line_color=None, fill_alpha=0.7)
        p2.add_glyph(source, glyph)
        p2.add_glyph(source, rects)

        p2.grid.visible = False
        p2.xaxis.major_label_text_font_style = "bold"
        p2.yaxis.minor_tick_line_width = 0
        p2.yaxis.major_tick_line_width = 0

        # hover action
        p2.hover.tooltips = [
            ("site number", "@x"),
            ("event", "@indel")
        ]

        #     p2.js_on_event(events.Tap, display_event(p2, attributes=point_attributes))

        # connection between two figures
        range_tool = RangeTool(x_range=p2.x_range)
        range_tool.overlay.fill_color = "gold"
        range_tool.overlay.fill_alpha = 0.4
        #     range_tool.overlay.

        p1.ygrid.grid_line_color = None
        p1.add_tools(range_tool)
        p2.toolbar.active_multi = 'auto'
        p1.toolbar.active_multi = 'auto'

        # -------------- callback function ----------------
        def callback(event):
            selected = source.selected.indices[0]
            selected_row_index = math.floor(selected / len_seq)

            result_t = []
            for k in range(0, nb_seq):
                if selected_row_index < k:
                    result_t += [selected + (len_seq * (k - selected_row_index))]
                elif selected_row_index > k:
                    result_t += [selected - (len_seq * (selected_row_index - k))]
                result_t += [selected]

            source.selected.indices = result_t

        taptool = p2.select(type=TapTool)
        p2.on_event(DoubleTap, callback)

        # source.selected.on_event("indices", my_tap_handler)

        # events.SelectionGeometry
        # events.Tap,
        # events.DoubleTap

        # source.js_on_event(events.DoubleTap, CustomJS(args=dict(s=source), code="""
        # selected = s.selected.indices[0]
        # selected_row_index = math.floor(selected/len_seq)

        # result_t = []
        # for k in range(0,nb_seq):
        # 	if selected_row_index < k:
        # 		result_t += [selected + (len_seq * (k-selected_row_index))]
        # 	elif selected_row_index > k:
        # 		result_t += [selected - (len_seq * (selected_row_index-k))]
        # 	result_t += [selected]

        # s.selected.indices = result_t

        # """))

        code = """
        var selected = source.selected.indices[0],
            selected_row_index = Math.floor(selected/len_seq),
            result_t = [];

        for (let k = 0; k < nb_seq; k++) {
            if(selected_row_index < k)
                result_t.push(selected + (len_seq * (k-selected_row_index)));
            else if(selected_row_index > k)
                result_t.push(selected - (len_seq * (selected_row_index-k)));
            result_t.push(selected);
        }

        source.selected.indices = result_t;
        source.change.emit();
        """

        callback = CustomJS(args={'source': source, 'len_seq': len_seq, 'nb_seq': nb_seq}, code=code)
        # def call_back_column_number_tap(event):
        # 	func(text_input.value)
        # 	return callback

        p2.add_tools(TapTool(gesture='doubletap', callback=callback))

        text_input = TextInput(value="", max_width=250, title="Enter the sequence column number:")
        # text_input.js_on_change("value", CustomJS(code="""
        #     console.log('text_input: value=' + this.value, this.toString())
        # """))

        button = Button(label="Select Column", align="end", width=200, max_width=80, button_type="primary")
        button.js_on_click(CustomJS(
            args=dict(source=source, text_input=text_input, len_seq=len_seq, nb_seq=nb_seq, range_tool=range_tool),
            code="""
            var result_t = [],
                selected_column = parseInt(text_input.value);
            for (let k = 0; k < nb_seq; k++) {
                if(0 < k)
                    result_t.push(selected_column + (len_seq * k));
                result_t.push(selected_column);
            }
            console.log((selected_column-35,selected_column+34))
            if(selected_column>36){
                range_tool.x_range.start = selected_column-35;
                range_tool.x_range.end = selected_column+34;
            }

            source.selected.indices = result_t;
        """))

        def call_back_column_number(event):
            func(text_input.value)

        def call_back_column_number_tap(event):
            func(source.selected.indices[0])

        # button.on_click(func(15))
        button.on_event(ButtonClick, call_back_column_number)

        p2.on_event(DoubleTap, call_back_column_number_tap)

        selecte_column = row(text_input, button, max_width=plot_width)  # , sizing_mode="stretch_width"

        p = gridplot([[p1], [p2], [selecte_column]], toolbar_location='right')

        return p

    @staticmethod
    def view_tree(tr, tree_panel, fontsize='9pt', plot_width=800):
        """
        :param :tr The tree data structure from ETE3 library
        Usage:
        >> Viewing the phylogeny tree
        """
        ts = TreeStyle()
        ts.show_leaf_name = False
        ts.draw_aligned_faces_as_table
        #     circular_style.mode = "c" # draw tree in circular mode
        #     circular_style.scale = 20
        for n in tr.traverse():
            if n.is_leaf():
                n.add_feature("evoltype", "L")
            else:
                n.add_feature("evoltype", "S")

        tree_panel.object = tr.render("%%inline", tree_style=ts)
        #     tree_result.object = tr.render("%%inline",  tree_style=ts)
        return

    @staticmethod
    def view_indel(event):
        #result.object =
        tree_panel.object = get_image()

    @staticmethod
    def utils_io(event):
        pass

    @staticmethod
    def utils_locate_site(loc=1):
        pass

    @staticmethod
    def reorder_sequences(tree, aln):
        for rec in aln:
            print(rec.id)
