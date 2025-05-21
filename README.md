# LRG Distance/Pathfinding Code Documentation

## Table of contents
1. [Introduction](#introduction)
2. [Installation](#installation)
    1. [Requirements](#requirements)
    2. [Walkthrough](#walkthrough)
3. [Usage](#usage)
    1. [Initialization](#initialization)
    2. [Data Storage and Access](#datastorage)

## Introduction <a name="introduction"></a>

This package contains the C code used to run efficient pathfinding (dijkstra's method) calculations over the longrungrowth data-grid to find the optimal path and travel cost between any two locations. The implementation is seperated into two parts covering the calculation of compatible bilateral cost matrix for use in the model simulation portion (link to py-sim!) as well as the tools to recreate the optimization of the parameters for the cost function against a set of training data with real-world grid-to-grid travel times/costs.

The package contains all of the C code for the above mentioned functions but also serves as an integration point with Python, allowing for the easy implementation of these functions within a python script or larger data processing pipeline while taking care of most of the tricky technical details behind the scenes.

## Installation <a name="installation"></a>

### Requirements <a name="requirements"></a>

- git
- Python 3.11 or greater
- A python virtual environment manager 
  - i.e. pyenv w/ pyenv-virtualenv or conda
- gcc compiler
- GNU Make
- nlopt optimizer

### Walkthrough <a name="walkthrough"></a>

1. Clone the repo into the desired location
2. Install the required python packages (ideally within a virtual environment):
  - numpy 
3. From the project root directory run `make all` in the terminal to compile the C++/CUDA source code into shared libraries (.so files) in the lib folder
4. If not already present make a folder called `data` in project root directory, this is where input and output data will be placed for the model
5. To call the model from the associated python module either create a python script in the `src` directory and run it from the project root directory with `python src/<script_name.py>` or load model interface from a REPL or interactive shell started from the project root directory (in which case the core module should be imported with `from src.dist.core import *`)

## Usage <a name="usage"></a>

In the abstract, the pathfinding system operates on three intersecting data structures, a grid of cells forming the extent of geographic data that is input to the model (this will likely be a grid covering the entire world due to data formatting), a grid of cells that forms a subset of the data grid that forms the geographic boundaries that the pathfinding model will be operating within (this will generally be a much larger scale than the data grid, i.e. many cells of the data grid will be contained within one cell of the model grid, even if the *total* model grid is a smaller than the *total* data grid), and finally a list of actual nodes composed of a subset of the model grid cells that will form the nodes of the model graph. The graph is then formed through finding connections and calculating the edge weight of travelling between the nodes of the model.

### Initialization <a name="initialization"></a>

Within a script, the core model implementation can be imported with
```python
from dist.core import *
```

Then, it is necessary to create a `Model` object with
```python
<model_name> = Model(data_grid_rows: int, data_grid_columns: int, 
                     model_grid_rows: int, model_grid_columns: int,
                     model_grid_ul_y: int, model_grid_ul_x: int,
                     model_grid_lr_y: int, model_grid_lr_x: int,
                     model_grid_size: int, model_node_count: int)
```
Where the dimensions (number of rows and columns) of both the full grid of data (most likely a grid covering the entire world) and the subset of that grid that forms the boundaries of the relevant model. Additionally, the size of the model grid (i.e. the number of data grid cells in each dimension that form one model grid cell) and the x, y index (referring to columns and rows respectively, starting from zero, in the upper left corner of the grid) of the upper left and lower right corners of the model grid within the data grid. 

Furthermore, the number of actual nodes that will form the model graph (assuming not all cells in the model grid will be considered as full nodes) must be specified. The actual details of these nodes must then be added from an external data input and loaded into the model, this first step only allocates the model arrays and forms the basic API structure.

### Data Storage and Access <a name="datastorage"></a>

The core `Model` object allocates and makes available the following fields of memory, accessible as numpy arrays. From the perspective of the end python user, the data is arranged in the following structure with three object members in the form of python dictionaries containing arrays:

1. `edge_weights` -- contains information on the current edge weights for the model graph
    1. `'neighbor_id'`: for each `node_id` index, contains a list of 8 ids of neighboring nodes in the graph 
    2. `'cost'`: a parallel array containing the edge weight 'cost' for each neighboring pair of `node_id`s in the graph
2. `node_data` -- contains information on the nodes (and non-node cells) of the model grid 
    1. `'node_id'`: contains either the `node_id` of the cell if it is a node or -1 if the cell is not a proper node 
    2. `'count'`: for each cell contains the number of proper node cells within the immediate adjacent cells of the grid, including diagonals (max 8)
    3. `'expected'`: for each cell contains the 'expected' number of proper node neighbors if all cells of the grid were nodes. Notably this distinguishes cells on the edges of the model grid as they have less 'natural' neighbors
3. `grid_data` -- contains information on the data grid
    1. `'elevation'`: for each cell contains the average elevation in meters above sea levels (note this may be negative)
    2. `'terrain'`: for each cell indicates whether the cell is considered 'land' (1) or 'sea' (0)
    3. `'river'`: for each cell indicates whether or not the cell contains a major river, 1 for yes, 0 for no
    4. `'cell_weight'` for each cell contains the relative weighting of that cell within the larger grouping of each model grid





