![Tux, the Indelviewer mascot](./docs/icons/indelViewer%20logo.jpg)

# Indel Points Visualization

This repo contains an implementation of the `arpip indel viewer` which help the user to visualize the indel points in the phylogenetic tree. 

## Getting Started

### Installation
To install the package you can do any of the following:

Simply download the repository and run the following command in the root directory:

### Usage
Using a phylogenetic tree, the user can visualize the indel points in the tree. The user can also visualize the indel points in the tree using the command line.

for run the application use this command:

```console
python ./main.py ./input.txt
```

for print a single site tree with indel points use this command:
```console
python ./_cli.py --single-result 102  ./input.txt
```

for print all site tree with indel points use this command:

```console
python ./_cli.py  ./input.txt             
```

for help use this command:

```console
python ./_cli.py —help
```

## Articles
* For more information on IndelPoints algorithm, see Jowkar et al.'s 2022 paper
which can be located [here](https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syac050/6648472).


## Contributing
If you want to contribute to this repo simply submit a pull request!


## Author
Gholam-Hossein Jowkar [E-mail](jowk@zhaw.ch)
