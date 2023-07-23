import argparse
import pathlib

from main import Initializer
from viewer.helper import Helper
from viewer.indel_tree import DrawTree
from viewer.plotting_events import PlottingEvents

from tqdm import tqdm

from ete3 import PhyloTree


def cmd_add(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--single-result', dest='singleResult', type=lambda x: int(x) if int(x) >= 1 else parser.error(
        '--single-result N must be zero or positive'), metavar='N', help='Show the scenario of entered column number')
    parser.add_argument('inputFile', type=pathlib.Path)
    parsed_args = parser.parse_args(args)

    with open(parsed_args.inputFile) as file:
        resources = {}
        for line in file:
            key, value = line.rstrip().split('=', 1)
            resources[key] = value
        # Or as a dictionary comprehension:
        # resources2 = {key:value for line in file for key, value in [line.rstrip().split('=', 1)]}

        necessary_keys = ['fasta', 'relation', 'tree', 'indel']
        input_file_keys = resources.keys()
        print(resources)

        if len(necessary_keys) == len(input_file_keys):
            if len([x for x in necessary_keys if x in input_file_keys]) == 4:

                if (resources['fasta'].split('.', -1)[-1] == 'fasta' and
                        resources['relation'].split('.', -1)[-1] == 'txt' and
                        resources['tree'].split('.', -1)[-1] == 'nwk' and
                        resources['indel'].split('.', -1)[-1] == 'txt'):

                    initializer = Initializer(resources)


                else:
                    print("Input text file is incorrect")
            else:
                print("Input text file is incorrect")
        else:
            print("Input text file is incorrect")

        if parsed_args.singleResult:
            plottingTree = PhyloTree(resources['tree'], alg_format='fasta', format=1)
            plottingEvents = PlottingEvents(plottingTree, initializer)
            plottingEvents.print_selected_site(parsed_args.singleResult)
        else:
            plottingTree = PhyloTree(resources['tree'], alg_format='fasta', format=1)
            plottingEvents = PlottingEvents(plottingTree, initializer, tqdm)
            plottingEvents.print_all_sites()


if __name__ == '__main__':
    cmd_add()


