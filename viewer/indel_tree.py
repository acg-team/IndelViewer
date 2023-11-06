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
# indel_tree.py
# Placing events on the tree data structure
# Created by: Gholam-hossein Jowkar
# ACGT ZHAW
# Created date: November 2022
# Modified by: Gholam-Hossein Jowkar
# Modified date: February 2023

__all__ = ['IndelTree']

import logging

import numpy as np
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace


class IndelTree:
    def event_on_tree(self, tree, ins_event, del_event, dic_name_index, child_lst):
        # extract event and map them to their real name
        deleted_nodes_id = np.where(del_event == 1)[0]
        deleted_nodes = [self.get_key_dic(ids, dic_name_index) for ids in deleted_nodes_id]
        # print("Deleted nodes", deleted_nodes)
        inserted_node_id = np.where(ins_event == 1)[0]
        inserted_node = [self.get_key_dic(ids, dic_name_index) for ids in inserted_node_id]
        # print("Insertion node", inserted_node)

        # if the insertion node is not the tree root, the node before insertion point should be considered as got deleted
        # print("insertion node id",inserted_node_id)

        is_exist_flag = 0

        subtree_del = []
        for n in tree.traverse(strategy='preorder'):
            n.add_feature("is_inserted", False)
            n.add_feature("is_deleted", False)
            n.add_feature("is_exist", False)  # still birth nodes

            if n.name in inserted_node:
                # print("The insertion point:",n.name)
                # insertion node
                n.is_inserted = True
                n.is_exist = True
                is_exist_flag = 0
                root_name = n.name
                # If the insertion is not the root, all ancestor node should
                # be flaged as deleted or other flags for the style section.
                # if inserted_node !='root':
                # print("Warning: the insertion point is not the root of the tree")
            elif n.name in deleted_nodes:
                # deleted nodes
                n.is_deleted = True
                if not n.is_leaf():
                    subtree_del = child_lst[n.name]
            elif n.name in subtree_del:
                # node extant of the deleted nodes
                n.is_deleted = True
                n.is_exist = True
            elif is_exist_flag and not (n.is_deleted):
                n.is_exist = True
                # nodes that exist but is not member of is_inserted but might be is_deleted
            elif n.name in child_lst[inserted_node[0]]:
                n.is_exist = True

    def get_key_dic(self, val, dic):
        for key, value in dic.items():
            if val == value:
                return key
        pass

    def father_subtree_node(self, tree):
        lst_child = {}
        for node in tree.traverse():
            father_node = node
            lst_child[father_node.name] = [node.name for node in father_node.iter_descendants()]
        return lst_child


class DrawTree(IndelTree):
    def __init__(self, ts, initializer, selected_site_number):
        # selected_site_number = 1188 #2546
        ins_event = initializer.mat_insertion[selected_site_number, :]
        del_event = initializer.mat_deletion[selected_site_number, :]
        child_lst = IndelTree.father_subtree_node(self, initializer.tree)
        result = IndelTree.get_key_dic(self, 1, initializer.dic_name_index)
        IndelTree.event_on_tree(self, initializer.tree, ins_event, del_event, initializer.dic_name_index, child_lst)
        logging.debug("[DrawTree:__init__] %s", child_lst)
        # for n in tree.traverse():
        #     print("The node %s is inserted: %s" %(n.name, n.is_inserted))
        #     print("The node %s is deleted: %s" %(n.name, n.is_deleted))
        #     print("The node %s is exist: %s" %(n.name, n.is_exist))

        # ts = TreeStyle()
        # ts.optimal_scale_level = 'full'
        ts.extra_branch_line_type = 2
        ts.extra_branch_line_color = 'Black'
        ts.layout_fn = self.indel_view_layout
        ts.show_leaf_name = False
        # tree.render("%%inline",  tree_style=ts)

    def indel_view_layout(self, node):
        if node.is_inserted:
            # set bold red branch to the insertion node
            ins_style = NodeStyle()
            ins_style["fgcolor"] = "#0f0f0f"
            ins_style["size"] = 0
            ins_style["vt_line_color"] = "#ff0000"
            ins_style["hz_line_color"] = "#ff0000"
            ins_style["vt_line_width"] = 6
            ins_style["hz_line_width"] = 6
            ins_style["vt_line_type"] = 0  # 0:solid, 1:dashed, 2:dotted
            ins_style["hz_line_type"] = 0
            node.set_style(ins_style)
            faces.add_face_to_node(faces.AttrFace("name", "Arial", 11, "#ff0000", None), node, 0)
        elif node.is_deleted:
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
            node.set_style(del_style)
            faces.add_face_to_node(faces.AttrFace("name", "Arial", 8, "#777777", None), node, 0)
        elif (not node.is_exist) and (not node.is_deleted):
            still_brth_style = NodeStyle()
            still_brth_style["fgcolor"] = "#777777"
            still_brth_style["vt_line_color"] = "#777777"
            still_brth_style["hz_line_color"] = "#777777"
            still_brth_style["vt_line_width"] = 2
            still_brth_style["hz_line_width"] = 2
            still_brth_style["hz_line_type"] = 2
            still_brth_style["vt_line_type"] = 2
            node.set_style(still_brth_style)
            faces.add_face_to_node(faces.AttrFace("name", "Arial", 8, "#777777", None), node, 0)
        else:
            node.img_style["fgcolor"] = "#0f0f0f"
            node.img_style["vt_line_color"] = "#0f0f0f"
            node.img_style["hz_line_color"] = "#0f0f0f"
            node.img_style["vt_line_width"] = 2
            node.img_style["hz_line_width"] = 2
            node.img_style["hz_line_type"] = 0
            node.img_style["vt_line_type"] = 0
            faces.add_face_to_node(faces.AttrFace("name", "Arial", 8, "#0f0f0f", None), node, 0)
