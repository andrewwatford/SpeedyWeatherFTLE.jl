# SpeedyWeatherFTLE

[![Build Status](https://github.com/andrewwatford/SpeedyWeatherFTLE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/andrewwatford/SpeedyWeatherFTLE.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![docs](https://img.shields.io/badge/documentation-latest_release-blue.svg)](https://andrewwatford.github.io/SpeedyWeatherFTLE.jl/)

## Setting up the project for development
To set up the project for development for the first time, first clone this repository. Then, from the repository directory, open `julia` and run:
```julia
]instantiate
]activate .
```
`instantiate` creates a `Manifest.toml` file, which creates the environment that is needed for all future use of the package - you usually only need to instantiate once. `activate .` activates the environment in the present directory, that is, the current project.

## Running the test suite located in `./test`
Run the following code:
```julia
]test
```

## TODOs:
- Publish this as a package so that it can be installed directly?
- Make some actual documentation
- Visualization functions:
    - Globe function (with slider?)
    - Slider through time
- From Ambre's code:
    - currently only supports 1-layer simulations
    - add ability to prescribe velocity field
    - Suppress output when advecting
    - Would be nice to get output directly rather than write to then read from file. Cannot find how to do this in doc
- Write nice docstrings for everything
- General refactoring
