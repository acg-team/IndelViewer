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
# version.py
# Created by: Gholam-Hossein Jowkar
# ACGT ZHAW
# Created date: November 2022
# Modified by:
# Modified date:

import importlib.metadata


try:
	__version__ = importlib.metadata.version('arpip_indel_viewer')
except importlib.metadata.PackageNotFoundError:
	__version__ = None

