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
# utils.py
# Created by: Gholam-hossein Jowkar
# ACGT ZHAW
# Created date: July 2023
# Modified by:
# Modified date:


import os, sys

import numpy as np
import pandas as pd

from ete3 import Tree

from Bio import SeqIO, AlignIO, Seq, SeqRecord


class TreeUtils:

    dict_distances_to_root = {}
    total_tree_length = 0
    def __init__(self, tree, dict_relation):
        self.tree = tree
        self.dict_relation = dict_relation

    def get_all_distances_to_root(self):
        """
        Get the distance of each node to the root of the tree
        :return:
        """
        for node in self.tree.traverse('postorder'):
            if node not in self.tree.iter_leaves():
                a = node.children[0]
                node.name = self.dict_relation.get(a.name)
            self.dict_distances_to_root[node.name] = self.tree.get_distance(node.name)
        return self.dict_distances_to_root

    @staticmethod
    def get_relationship_nodes(relation_lines):
        """
        Extract the relationship between nodes from the relation file
        :param relation_lines: Lines containing the relationship between nodes from ARPIP output
        :return: dictionary of nodes and their parents (kid:father)
        """
        dict_node_relationship = {}
        for line in relation_lines[1:-1]:
            s = line.split('\t')
            dict_node_relationship[s[0]] = s[1].strip()
        return dict_node_relationship

    def compute_total_tree_length(self):
        """
        Compute the total tree length
        :return None
        """
        for node in self.tree.traverse('postorder'):
            self.total_tree_length += node.dist

    def print_relationship_nodes(self):
        """
        Print the relationship between nodes
        :return None
        """
        if self.dict_relation:
            print("Relation file (kid:father)", self.dict_relation)
        else:
            print("No relation file is provided, please run the get_relationship_nodes method first.")



class AlignmentUtils:
    dict_alignment_lengths = {}
    def __init__(self, alignment):
        self.alignment = alignment

    def get_all_alignment_lengths(self, with_gap=False):
        """
        Compute the length of each sequence in the alignment
        :param with_gap: to consider the gap characters in the sequence or not
        :return: dictionary of sequence ids and their lengths
        """
        for record in self.alignment:
            if with_gap:
                self.dict_alignment_lengths[record.id] = len(record.seq)
            else:
                self.dict_alignment_lengths[record.id] = len(record.seq.ungap('-'))
        return self.dict_alignment_lengths

    def split_multiple_ancestral_alignment(self, asr_pref="V"):
        """
        Split the alignment file into msa and maa using the prefix provided by user.
        :param asr_pref: The prefix which the default value is "V". The prefix is used to identify the ASR sequences.
        :return: msa file and maa file separately.
        """
        msa = []
        asr = []
        for rec in self.alignment:
            try:
                if asr_pref in rec.id or "root" in rec.id:
                    asr.append(rec)
                else:
                    msa.append(rec)
            except KeyError as ker:
                print("ERROR: In the alignment file %s", str(ker))
        MSA = AlignIO.MultipleSeqAlignment(msa)
        ASR = AlignIO.MultipleSeqAlignment(asr)
        return MSA, ASR

    def get_gappy_cols(self):
        """
        Extract gappy columns from the alignment
        :param msa: alignment file in AlionIO format (see BioPython for more)
        return: list of unique column numbers containing gap characters in the alignment.
        """
        gappy_col = []
        # row counter
        r_counter = 0
        for record in self.alignment:
            # col counter
            c_counter = 0
            for i in record:
                if i == '-':
                    gappy_col.append(c_counter)
                c_counter += 1
            r_counter += 1
        return list(set(gappy_col))

    # @njit it is possible to run this function with numba
    @staticmethod
    def get_miss_match(np_inferred_alignment, np_ref_alignment):
        """
        The total number of inferred characters similarly, excluding gaps, divided by the total number of characters
        :param np_inferred_alignment:
        :param np_ref_alignment:
        :return:
        """
        lst_mis_match_lin = []
        total_sub_lin = np.count_nonzero((np_ref_alignment != 0) & (np_ref_alignment != -1), axis=1)
        mis_match_lin = 0
        mis_match = 0
        total_sub = np.count_nonzero((np_ref_alignment != 0) & (np_ref_alignment != -1))
        for i in range(np_inferred_alignment.shape[0]):
            mis_match_lin = 0
            for j in range(np_inferred_alignment.shape[1]):
                if np_ref_alignment[i, j] != 0 and np_ref_alignment[i, j] != -1 and np_inferred_alignment[i, j] != np_ref_alignment[i, j]:
                    mis_match += 1
                    mis_match_lin += 1
            lst_mis_match_lin.append(mis_match_lin)
        result = [x/y for x, y in zip(lst_mis_match_lin, total_sub_lin)]
        return total_sub, mis_match, 1 - (mis_match/total_sub)

    @staticmethod
    def convert_fasta_to_phylip(fasta_aln_path, phylip_aln_path):
        """
        Convert fasta alignment file to phylip format
        :param fasta_aln_path: path and name of to fasta alignment file
        :param phylip_aln_path: path and name of phylip alignment file
        :return: None
        """
        AlignIO.convert(fasta_aln_path, "fasta", phylip_aln_path, "phylip-relaxed")

    @staticmethod
    def convert_phylip_to_fasta(phylip_aln_path, fasta_aln_path):
        """
        Convert phylip alignment file to fasta format
        :param phylip_aln_path:
        :param fasta_aln_path:
        :return: None
        """
        AlignIO.convert(phylip_aln_path, "phylip-relaxed", fasta_aln_path, "fasta")