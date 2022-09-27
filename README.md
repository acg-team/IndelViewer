# Indel Points Visualization

This repo contains an implementation of the `arpip indel viewer` which... 
## Contributing
If you want to contribute to this repo simply submit a pull request!

## Getting Started

### Installation
To install the package you can do any of the following:

- Run the command `pip install ...`

### Using RobustRandomCutForests
Using a RobustRandomCutForest to classify potential anomalies in your data is simple. Assuming you already have a vector of data stored in `X` you would run the following:

```python
from arpip.indelpoints.visualizer import  
tree = ....
```


The function `` will return an 

```python
# Given an array of points....
for point in points:
    
```

## Testing
All tests are written using `pytest`. Simply `pip install pytest` to be able to run tests. All tests are located under the `tests` folder. Any new tests are always welcome!

## Articles
* For more information on IndelPoints algorithm, see Jowkar et al.'s 2022 paper
which can be located [here](https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syac050/6648472).


## Contact
<jowk@zhaw.ch>
