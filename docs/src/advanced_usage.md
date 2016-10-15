# Advanced Usage

```@meta
CurrentModule = DSGE
```

## Package Directory Structure

The package directory structure follows Julia module conventions. Directories in square brackets indicate future additions. *Note that this directory tree is not linked, although it appears to be.*

```@contents
Pages = ["pkg_structure.md"]
Depth = 5
```

## Working with Settings

There are many computational settings that affect how the code runs without affecting the
mathematical definition of the model. 

Below, we describe several important settings for package usage.

For more details on implementation and usage of settings, see [Model Settings](@ref).

See [defaults.jl](https://github.com/FRBNY-DSGE/DSGE.jl/blob/master/src/defaults.jl) for the complete description of default settings.

#### General

- `dataroot`: The root directory for
  model input data.
- `saveroot`: The root directory for model output.
- `use_parallel_workers`: Use available parallel workers in computaitons.
- `data_vintage`: Data vintage identifier, formatted
  `yymmdd`. By default, `data_vintage` is set to today's date. It is (currently) the only
  setting printed to output filenames by default.

#### Dates
- `date_presample_start`: Start date of pre-sample.
- `date_mainsample_start`: Start date of main sample.
- `date_zlbregime_start`: Start date of zero lower bound regime.
- `date_mainsample_end`: End date of main sample.
- `date_forecast_start`: Start date of forecast period.
- `date_forecast_end`: End date of forecast period.

#### Anticipated Shocks
- `n_anticipated_shocks`: Number of anticipated policy shocks.
- `n_anticipated_shocks_padding`: Padding for anticipated shocks.

#### Estimation
- `reoptimize`: Whether to reoptimize the posterior mode. If `true`
    (the default), `estimate()` begins reoptimizing from the model
    object's parameter vector. See [Optimizing or Reoptimizing](@ref
    estimation-reoptimizing) for more details.
- `calculate_hessian`: Whether to compute the Hessian. If `true` (the
    default), `estimate()` calculates the Hessian at the posterior mode.

#### Metropolis-Hastings
- `n_mh_simulations`: Number of draws from the posterior distribution per block.
- `n_mh_blocks`: Number of blocks to run Metropolis-Hastings.
- `n_mh_burn`: Number of blocks to discard as burn-in for Metropolis-Hastings.
- `mh_thin`: Metropolis-Hastings thinning step.

### Accessing Settings
The function `get_setting(m::AbstractModel, s::Symbol)` returns the value of the setting `s`
in `m.settings`. Some settings also have explicit getter methods that take only the model
object `m` as an argument. Note that not all are exported.

*I/O*:

  * `saveroot(m)`,
  * `dataroot(m)`,
  * `data_vintage(m)`,

*Parallelization*:

  * `use_parallel_workers(m)`

*Estimation*:

  * `reoptimize(m)`,
  * `calculate_hessian(m)`,

*Metropolis-Hastings*:

  * `n_mh_blocks(m)`,
  * `n_mh_simulations(m)`,
  * `n_mh_burn(m)`,
  * `mh_thin(m)`

### Overwriting Default Settings

To overwrite default settings added during model construction, a user must define a new
`Setting` object and update the corresponding entry in the model's `settings` dictionary
using the `<=` syntax. If the `print`, `code`, and `description` fields of the new `Setting`
object are not provided, the fields of the existing setting will be maintained. If new
values for `print`, `code`, and `description` are specified, and if these new values are
distinct from the defaults for those fields, the fields of the existing setting will be
updated.

For example, overwriting `use_parallel_workers` should look like this:
```julia
m = Model990()
m <= Setting(:use_parallel_workers, true)
```

## [Editing or Extending a Model](@id editing-extending-model)

Users may want to extend or edit `Model990` in a number of different ways.  The most common
changes are listed below, in decreasing order of complexity:

1. Add new parameters
2. Modify equilibrium conditions or measurement equations
3. Change the values of various parameter fields (i.e. initial `value`, `prior`,
   `transform`, etc.)
4. Change the values of various computational settings (i.e. `reoptimize`,
   `n_mh_blocks`)

Points 1 and 2 often go together (adding a new parameter guarantees a change in equilibrium
conditions), and are such fundamental changes that they increment the model specification
number and require the definition of a new subtype of `AbstractModel` (for instance,
`Model991`).  See [Model specification](@ref model-specification-mspec) for more details.

Any changes to the initialization of preexisting parameters are defined as a new model
*sub-specification*, or *subspec*. While less significant than a change to the model's
equilibrium conditions, changing the values of some parameter fields (especially priors) can
have economic significance over and above settings we use for computational purposes.
**Parameter definitions should not be modified in the model object's constructor.** First,
incrementing the model's sub-specification number when parameters are changed improves
model-level (as opposed to code-level) version control. Second, it avoids potential output
filename collisions, preventing the user from overwriting output from previous estimations
with the original parameters. The protocol for defining new sub-specifications is described
in [Model sub-specifications](@ref model-sub-specifications-msubspec).

Overriding default settings is described in the [Model Settings](@ref) section.

### [Model specification (`m.spec`)](@id model-specification-mspec)

A particular model, which corresponds to a subtype of `AbstractModel`, is defined as a set
of parameters, equilibrium conditions (defined by the `eqcond` function) and measurement
equations (defined by the `measurement` function).  Therefore, the addition of new
parameters, states, or observables, or any changes to the equilibrium conditions or
measurement equations necessitate the creation of a new subtype of `AbstractModel.`

To create a new model object, we recommend doing the following:

1. Duplicate the `m990` directory within the [models](https://github.com/FRBNY-DSGE/DSGE.jl/tree/master/src/models) directory. Name the new
   directory `mXXX.jl`, where `XXX` is your chosen model specification number or string.
   Rename `m990.jl` in this directory to `mXXX.jl`.
2. In the `mXXX/` directory, change all references to `Model990` to `ModelXXX`.
3. Edit the `m990.jl`, `eqcond.jl`, and `measurement.jl` files as you see fit.  If adding
   new states, equilibrium conditions, shocks, or observables, be sure to add them to the
   appropriate list in `init_model_indices`.
4. Open the module file, `src/DSGE.jl`. Add `ModelXXX` to the list of functions to export,
   and include each of the files in `src/model/mXXX`.

### [Model sub-specifications (`m.subspec`)](@id model-sub-specifications-msubspec)

`Model990` sub-specifications are initialized by overwriting initial parameter definitions
before the model object is fully constructed. This happens via a call to `init_subspec` in
the `Model990` constructor. (Clearly, an identical protocol should be followed for new model
types as well.)

To create a new sub-specification (e.g., subspec 1) of `Model990`, edit the file
`src/models/subspecs.jl` as follows (note that this example is not actually
sub-specification `1` of `Model990`. In the source code, our sub-specification `5` is
provided as additional example.):

**Step 1.** Define a new function, `ss1`, that takes an object of type `Model990` (not
   `AbstractModel`!) as an argument. In this function, construct new parameter objects and
   overwrite existing model parameters using the `<=` syntax. For example,

```julia
function ss1(m::Model990)
    m <= parameter(:ι_w, 0.000, (0.0, .9999), (0.0,0.9999), DSGE.Untransformed(), Normal(0.0,1.0), fixed=false,
                   description="ι_w: Some description.",
                   tex_label="\\iota_w")
    m <= parameter(:ι_p, 0.0, fixed=true,
                   description= "ι_p: Some description"
                   tex_label="\\iota_p")
end
```

**Step 2.** Add an `elseif` condition to `init_subspec`:

```julia
    ...
    elseif subspec(m) == "ss1"
        return ss1(m)
    ...
```
To construct an instance of `Model990`, `ss1`, call the constructor
for `Model990` with `ss1` as an argument. For example,

```julia
m = Model990("ss1")
```