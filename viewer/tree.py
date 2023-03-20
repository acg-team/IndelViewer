########################################################################################################################
# Licensed Materials - Property of Gholamhossein Jowkar
#
#
# Copyright (C) 2022-2025 by Gholamhossein Jowkar
#
# This work was supported by the Swiss National Science Foundation (SNF) grants $31003A\_176316$ to Dr. M. Anisimova.
# The funding body did not play any role in the design of the study and collection, analysis, and interpretation of
# data and in writing the code.
#
########################################################################################################################
#
# This file is part of IndelViewer of ARPIP project: https://github.com/acg-team/IndelViewer

#                               ABOUT THE ARPIP PACKAGE
#                               =====================
# ARPIP: Ancestral Sequence Reconstruction with insertions and deletions under the Poisson Indel Process
# ARPIP is a joint maximum likelihood approach for phylogenetic ancestral sequence reconstruction, capable of modeling
# indels both biological and mathematically.
#

# IndelViewer is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# IndelViewer is a free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General
# Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# You should have received a copy of the GNU Lesser General Public
# License along with IndeViewer. If not, see <http://www.gnu.org/licenses/>.
#
########################################################################################################################

# Python 3
# Processing Tree file
# Created by:  Gholamhossein Jowkar <jowk@zhaw.ch>
# ACGT ZHAW
# Created date: November 2022
# Modified by: Gholamhossein Jowkar
# Modified date:
__all__ = ['Tree']

import sys, os
import numpy as np
# form ETE3
from ete3 import Tree, PhyloTree, TreeStyle, NodeStyle, SeqGroup, faces, AttrFace, CircleFace


class Tree:
    def event_on_tree(self, tree, ins_event, del_event, dict_name_index, lst_child):
        """
        This function place the indel events on the tree.
        :param tree: The ete3 Tree
        :param ins/del_event: The numpy array (column) corresponding to indel event with values of zero and one for
        presence and absence of events.
        :param dict_name_index: A Dictionary of 'value:index'.
        :param lst_child: The child name for the subtree we're working on in format of 'father: child'.
        :return: Nothing, but placing events on the tree.
        """
        # extract event and map them to their real name
        deleted_nodes_id = np.where(del_event == 1)[0]
        deleted_nodes = [self.get_key_dic(ids, dict_name_index) for ids in deleted_nodes_id]
        # print("Deleted nodes", deleted_nodes)
        inserted_node_id = np.where(ins_event == 1)[0]
        inserted_node = [self.get_key_dic(ids, dict_name_index) for ids in inserted_node_id]
        # print("Insertion node", inserted_node)

        # if the insertion node is not the root, the node before insertion point should be considered as got deleted
        print("insertion node id", inserted_node_id)

        is_exist_flag = 0

        subtree_del = []
        for n in tree.traverse(strategy='preorder'):
            n.add_feature("is_inserted", False)
            n.add_feature("is_deleted", False)
            n.add_feature("is_exist", False)  # still birth nodes

            if n.name in inserted_node:
                print("The insertion point:", n.name)
                # insertion node
                n.is_inserted = True
                n.is_exist = True
                is_exist_flag = 0
                root_name = n.name
                # If the insertion is not the root, all ancestor node should
                # be flaged as deleted or other flags for the style section.
                if inserted_node != 'root':
                    print("Warning: the insertion point is not the root of the tree")
            elif n.name in deleted_nodes:
                # deleted nodes
                n.is_deleted = True
                if not n.is_leaf():
                    subtree_del = lst_child[n.name]
            elif n.name in subtree_del:
                # node extant of the deleted nodes
                n.is_deleted = True
                n.is_exist = True
            elif is_exist_flag and not (n.is_deleted):
                n.is_exist = True
                # nodes that exist but is not a member of is_inserted but might be is_deleted.
            elif n.name in lst_child[inserted_node[0]]:
                n.is_exist = True

    def get_key_of_dict(self, val, dict):
        """
       Getting the key value of any arbitrary dictionary.
       """
        for key, value in dict.items():
            if val == value:
                return key

    def get_father_subtree_of_node(self, tree):
        """
        For each subtree rooted at all the nodes, extracting a dictionary of child's name.
        :param tree: The ete3 tree.
        :return: list of child a subtree
        """
        lst_child = {}
        for node in tree.traverse():
            father_node = node
            lst_child[father_node.name] = [node.name for node in father_node.iter_descendants()]
        return lst_child



