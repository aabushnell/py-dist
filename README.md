# LRG Distance/Pathfinding Code Documentation

## Table of contents
1. [Introduction](#introduction)
2. [Installation](#installation)
    1. [Requirements](#requirements)
    2. [Walkthrough](#walkthrough)
3. [Usage](#usage)
    1. [Initialization](#initialization)

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
5. To call the model from the associated python module either create a python script in the `src` directory and run it from the project root directory with `python src/<script_name.py>` or load model interface from a REPL or interactive shell started from the project root directory (in which case the core module should be imported with `from src.sim.core import *`)

## Usage <a name="usage"></a>

### Initialization <a name="initialization"></a>
