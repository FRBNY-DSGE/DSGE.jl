# [`DSGE.jl`](https://github.com/FRBNY-DSGE/DSGE.jl)
## `doc/`: Code and model documentation.
## `save/`: Sample input files; default input/output directories.
## `src/`
### `DSGE.jl`: The main module file.
### `abstractdsgemodel.jl`: Defines the `AbstractModel` type.
### `parameters.jl`: Implements the `AbstractParameter` type and its
       subtypes.
### `settings.jl`: Implements the `Setting` type.
### `defaults.jl`: Specifies default `Setting`s.
### `distributions_ext.jl`: DSGE-specific extensions of `Distributions` functionality.
### `data/`: Manipulating and updating input dataset.
### `solve/`: Solving the model; includes `gensys.jl` code.
### `estimate/`: Optimization, posterior sampling, and other functionality.
### [`xxx/`]: Other model functionality, such as forecasts, impulse response functions, and shock decompositions.
### `models/`
#### `m990/`: Contains code to define and initialize version 990 of the FRBNY DSGE model.
##### `m990.jl`: Constructs a `Model990` object.
##### `eqcond.jl`: Constructs `Model990` equilibrium condition matrices
##### `measurement.jl`: Constructs `Model990` measurement equation matrices.
##### `subspecs.jl`: Code for model sub-specifications is defined here. See [Editing or Extending a Model](@ref editing-extending-model) for details on constructing model sub-specifications.
##### `augment_states.jl`: Code for augmenting the state space system after model solution.
#### [`[m991/]`]: Code for new models should be kept in directories at this level in the directory tree
## `test/`: Module test suite.
