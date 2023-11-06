![Tux, the Indelviewer mascot](./docs/icons/indelViewer%20logo.jpg)

# ARPIP's Indel Visualizer

This repository contains the implementation of the `arpip indel viewer` which help the user to visualize the indel points
inferred by ARPIP in the phylogenetic tree. 

## Getting Started
Before installation make sure that you have the Python version 3.8 or higher installed on your system. 
[Need help?](https://realpython.com/installing-python/)

### Installation

To install the package you can simply download the repository and run the following command in the root directory.

Install the dependencies using this command:

```console
pip3 install -r requirements.txt
```

### Usage
The user can run the visualizer using the command line. For running the application use this command:

```console
python ./main.py ./input.txt
```

To print a single site tree with indel points use this command:

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
