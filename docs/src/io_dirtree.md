# `<dataroot>/`: Root data directory.
## `data/`:  Macroeconomic input data series.
## `cond/`: Conditional data, i.e.
    ["nowcast"](https://en.wikipedia.org/wiki/Nowcasting_%28economics%29).
## `user/`: User-created or sample model input files.

# `<saveroot>/`: Root save directory.
## `output_data/`
### `m990/`: Input/output files for the `Model990` type. A model of type `SPEC` will create its own save directory `SPEC/` at this  level in the directory tree.
### `ss0/`: Subdirectory for subspec 0. A model of a different subspec will have similar directories at this level of the tree.
#### `estimate/`
##### `figures/`: Plots and other figures
##### `tables/`: LaTeX tables
##### `raw/`: Raw output data from estimation step
##### `work/`: Derived data files created using `raw/` files as input
#### [`xxx/`]: Other model outputs, such as forecasts, impulse response functions, and shock decompositions; subdirectory structure mirrors that of `estimate`.
