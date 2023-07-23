# Indel Points Visualization

This repo contains an implementation of the `arpip indel viewer` which help the user to visualize the indel points in the phylogenetic tree. 

## Contributing
If you want to contribute to this repo simply submit a pull request!

## Getting Started

### Installation
To install the package you can do any of the following:

- Run the command 
```bash
pip install arpip_indel_viewer
```

### Usage
Using a phylogenetic tree, the user can visualize the indel points in the tree. The user can also visualize the indel points in the tree using the command line.

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


## Author
Gholom-Hossein Jowkar [E-mail](jowk@zhaw.ch)
