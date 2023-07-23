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
# plotting_events.py
# Plotting the events on the tree
# Created by: Gholam-hossein Jowkar
# ACGT ZHAW
# Created date: February 2023
# Modified by:
# Modified date:
# For more info see http://etetoolkit.org/docs/latest/tutorial/


# from ete3 import (PhyloTree, TreeStyle, NodeStyle, SeqGroup, TextFace, AttrFace, SequenceFace)
from ete3 import (NodeStyle, SeqGroup)
import os, sys, os.path
from IPython.display import Image
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from .indel_tree import IndelTree
import numpy as np


class PlottingEvents(IndelTree):
    def __init__(self, tr, initializer, tqdm=None):

        self.tr = tr
        self.initializer = initializer
        self.tqdm = tqdm
        ####### main code ########
        # The internal nodes'name is missing from this PhyloTree object, we add them manually
        for node in tr.traverse("postorder"):
            if node not in tr.iter_leaves():
                a = node.children[0]
                node.name = initializer.dic_relation.get(a.name)
                # print("The list of internal nodes:", node.name)
            initializer.list_names.append(node.name)
        # print("The list of all nodes:", initializer.list_names)

        align_fas = format(initializer.aligned_data, 'fasta')

        directory = "EventsPerSite/"
        self.path = os.getcwd()
        path_w_dir = os.path.join(self.path, directory)

        # define the access rights
        access_rights = 0o755

        # create new directory to store the
        if not (os.path.isdir(path_w_dir)):
            try:
                os.mkdir(path_w_dir)
                print("Directory '% s' is created" % directory)
            except OSError as error:
                print(error)
            else:
                print("Successfully created the directory %s " % path_w_dir)
        else:
            print("The directory is already exists.")

        os.chdir(path_w_dir)

    def print_selected_site(self, site_number):
        self.render_site(site_number, self.tr, self.initializer.aligned_data, self.initializer.mat_insertion,
                         self.initializer.mat_deletion, self.initializer.dic_name_index)

        print("Successfully change the directory back to %s " % self.path)
        os.chdir(self.path)

    def print_all_sites(self):
        self.render_all_sites(self.tr, self.initializer.aligned_data, self.initializer.mat_insertion,
                              self.initializer.mat_deletion, self.initializer.dic_name_index,
                              self.tqdm)  # all the sites in the msa file

        print("Successfully change the directory back to %s " % self.path)
        os.chdir(self.path)

    def link_to_align(self, tree, aln):
        """
        same function as for phyloTree
        :argument tree: the ete3 phylotree
        :argument aln the alignment file which should be linked with the tree

        """
        missing_leaves = []
        missing_intenals = []
        if type(aln) == SeqGroup:
            alg = aln
        else:
            alg = SeqGroup(aln, alg_format="fasta")  ###
        for n in tree.traverse():
            try:
                n.add_feature("sequence", alg.get_seq(n.name))
                # print(n.name)
            except KeyError:
                if n.is_leaf():
                    missing_leaves.append(n.name)
                else:
                    missing_intenals.append(n.name)

        if len(missing_leaves) > 0:
            print("Warning: [%d] terminal nodes could not be found in the alignment." % \
                  len(missing_leaves), file=sys.stderr)
        if len(missing_intenals) > 0:
            print("Warning: [%d] internal nodes could not be found in the alignment." % \
                  len(missing_leaves), file=sys.stderr)

    def render_site(self, site_number, tree, alignment, mat_insertion, mat_deletion, dic_name_index):
        """
        :param site_number: the number of the selected site
        :param tree: the ete3 phylotree
        :param alignment: the alignment file which should be linked with the tree
        """
        single_col_str = self.make_single_residue_msa(alignment, site_number)
        #     tr_style = apply_event_to_tree(site_number, tree)
        # print("The site number is", site_number)
        ins_event = mat_insertion[site_number, :]
        del_event = mat_deletion[site_number, :]
        child_lst = IndelTree.father_subtree_node(self, tree)
        self.reset_view(tree)
        IndelTree.event_on_tree(self, tree, ins_event, del_event, dic_name_index, child_lst)
        #     for n in tree.traverse():
        #         print("The node %s is inserted: %s" %(n.name, n.is_inserted))
        #         print("The node %s is deleted: %s" %(n.name, n.is_deleted))
        #         print("The node %s is exist: %s" %(n.name, n.is_exist))

        #     tss = TreeStyle()
        #     tss.optimal_scale_level = 'full'
        #     tss.extra_branch_line_type = 2
        #     tss.extra_branch_line_color= 'Black'
        # #     tss.layout_fn = indel_view_layout
        #     tss.show_leaf_name = False
        #     tss.draw_aligned_faces_as_table
        self.link_to_align(tree, single_col_str)
        self.print_view(tree)
        tree.render("tree_alignment_site_{}.png".format(site_number), w=1200, units='px')  # , tree_style = tss

    def render_all_sites(self, tree, alignment, mat_insertion, mat_deletion, dic_name_index, tqdm):
        """
        Saving the tree and msa columns for all the existing site in "png" format.
        :param tree: the ete3 phylotree
        :param alignment: the alignment file which should be linked with the tree
        """
        aln_len = alignment.get_alignment_length()
        for i in tqdm(range(0, aln_len)):
            self.render_site(i, tree, alignment, mat_insertion, mat_deletion, dic_name_index)

        print("All the sites were printed successfully.")

    def make_single_residue_msa(self, alignment, site_number):
        """
        Making text from a single column of msa with fasta format
        param: alignment: the alignment file which should be linked with the tree
        param: site_number: the selected site number which should be extracted from msa file
        """
        single_res = ""
        sites_count = 0
        # print(alignment[:, site_number])
        for record in alignment:
            single_res = single_res + ">" + record.id + "\n" + alignment[sites_count, site_number] + "\n"
            sites_count += 1
        return single_res

    def reset_view(self, tree):
        """
        Delete all the features and styles from last run
        since the style and feature are attached to the tree and
        per column these features and styles are different.
        """
        for n in tree.traverse():
            n.del_feature("is_inserted")
            n.del_feature("is_deleted")
            n.del_feature("is_exist")
            normal_stl = NodeStyle()
            normal_stl["fgcolor"] = "#000000"
            normal_stl["size"] = 0
            normal_stl["vt_line_color"] = "#000000"
            normal_stl["hz_line_color"] = "#000000"
            normal_stl["vt_line_width"] = 1
            normal_stl["hz_line_width"] = 1
            normal_stl["vt_line_type"] = 0  # 0:solid, 1:dashed, 2:dotted
            normal_stl["hz_line_type"] = 0
            n.set_style(normal_stl)

    def print_view(self, tree):
        """
        This function meant to apply events on the tree topology.
        The layout function could not change the style in this case,
        so I changed them manually.
        :param tree: the ete3 phylotree
        """
        for n in tree.traverse():
            if n.is_inserted:
                # set bold red branch to the insertion node
                # n.add_feature("evoltype", "D")
                ins_style = NodeStyle()
                ins_style["fgcolor"] = "#0f0f0f"
                ins_style["size"] = 0
                ins_style["vt_line_color"] = "#ff0000"
                ins_style["hz_line_color"] = "#ff0000"
                ins_style["vt_line_width"] = 6
                ins_style["hz_line_width"] = 6
                ins_style["vt_line_type"] = 0  # 0:solid, 1:dashed, 2:dotted
                ins_style["hz_line_type"] = 0
                n.set_style(ins_style)
                # faces.add_face_to_node(faces.AttrFace("name","Arial",11,"#ff0000",None), node, 0 )
            elif n.is_deleted:
                # n.add_feature("evoltype", "L")
                # set dashed black line to the deleted nodes
                del_style = NodeStyle()
                del_style["fgcolor"] = "#777777"
                del_style["shape"] = "circle"
                del_style["vt_line_color"] = "#777777"
                del_style["hz_line_color"] = "#777777"
                del_style["vt_line_width"] = 2
                del_style["hz_line_width"] = 2
                del_style["vt_line_type"] = 1  # 0:solid, 1:dashed, 2:dotted
                del_style["hz_line_type"] = 1
                n.set_style(del_style)
                # faces.add_face_to_node(faces.AttrFace("name","Arial",8,"#777777",None), node, 0 )
            elif (not n.is_exist) and (not n.is_deleted):
                # n.add_feature("evoltype", "L")
                still_brth_style = NodeStyle()
                still_brth_style["fgcolor"] = "#777777"
                still_brth_style["vt_line_color"] = "#777777"
                still_brth_style["hz_line_color"] = "#777777"
                still_brth_style["vt_line_width"] = 2
                still_brth_style["hz_line_width"] = 2
                still_brth_style["hz_line_type"] = 2
                still_brth_style["vt_line_type"] = 2
                n.set_style(still_brth_style)
                # faces.add_face_to_node( faces.AttrFace("name","Arial",8,"#777777",None), node, 0 )
            else:
                # n.add_feature("evoltype", "S")
                n.img_style["fgcolor"] = "#0f0f0f"
                n.img_style["vt_line_color"] = "#0f0f0f"
                n.img_style["hz_line_color"] = "#0f0f0f"
                n.img_style["vt_line_width"] = 2
                n.img_style["hz_line_width"] = 2
                n.img_style["hz_line_type"] = 0
                n.img_style["vt_line_type"] = 0
                # faces.add_face_to_node( faces.AttrFace("name","Arial",8,"#0f0f0f",None), node, 0 )




