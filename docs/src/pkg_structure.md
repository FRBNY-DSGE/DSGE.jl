# [`DSGE.jl`](https://github.com/FRBNY-DSGE/DSGE.jl)
## `doc/`: Code and model documentation.
### `examples/`: example scripts
## `save/`: Sample input files; default input/output directories.
## `src/`
### `DSGE.jl`: The main module file.
### `abstractdsgemodel.jl`: Defines the `AbstractModel` type.
### `abstractvarmodel.jl`: Defines the `AbstractVARModel`, `AbstractDSGEVARModel`, and `AbstractDSGEVECMModel` types.
### `defaults.jl`: Default settings for models.
### `statespace_types.jl`: Defines types for computing the state-space representation of models.
### `statespace_functions.jl`: Defines functions for computing the state-space representation of models.
### `data/`: Manipulating and updating input dataset.
### `solve/`: Solving the model; includes `gensys.jl` code.
### `estimate/`: Optimization, posterior sampling, and other functionality.
### `forecast/`: Forecasts, smoothing, shock decompositions, and impulse response functions.
### `decomp/`: Decompose changes in forecasts into three reasons: new data, data revisions, and changes in the calibration.
### `analysis/`: Moment tables of estimated parameters, computation of forecast means and bands.
### `altpolicy/`: Infrastructure for forecasting under alternative monetary policy rules.
### `scenarios/`: Forecasting alternative scenarios.
### `plot/`: Plot estimation results, forecasts, etc.
### `packet/`: Automatically generate documents with results from forecasts and estimations.
### `models/`
#### `representative/`: Representative agent models.
##### `m990/`: Contains code to define and initialize version 990 of the New York Fed DSGE model.
###### `m990.jl`: Constructs a `Model990` object.
###### `eqcond.jl`: Constructs `Model990` equilibrium condition matrices
###### `measurement.jl`: Constructs `Model990` measurement equation matrices.
###### `pseudo_measurement.jl`: Constructs `Model990` pseudo-measurement equation matrices.
###### `subspecs.jl`: Code for model sub-specifications is defined here. See [Editing or Extending a Model](@ref editing-extending-model) for details on constructing model sub-specifications.
###### `augment_states.jl`: Code for augmenting the state space system after model solution.
##### [`[m991/]`]: Code for new models should be kept in directories at this level in the directory tree
#### `heterogeneous/`: Heterogeneous agent models
#### `poolmodel/`: PoolModel type for model averaging.
#### `var/`: DSGE-VAR and DSGE-VECM models.
## `test/`: Module test suite.
