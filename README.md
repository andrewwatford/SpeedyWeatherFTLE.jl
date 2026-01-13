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
    - Make sure everyone is added as co-authors / maintainers in the .toml file
- Make some adaptive documentation
    - Make some documentation that runs the greet function
    - Make some documentation that plots something in SpeedyWeather
- Figure out what's going on with that PR
- Flesh out this TODO list :)
