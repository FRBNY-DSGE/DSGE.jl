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

- `saveroot::String`: The root directory for model output.
- `use_parallel_workers::Bool`: Use available parallel workers in computations.
- `n_anticipated_shocks`: Number of anticipated policy shocks.

#### Data and I/O

- `dataroot::String`: The root directory for model input data.
- `data_vintage::String`: Data vintage, formatted `yymmdd`. By default,
  `data_vintage` is set to today's date. It is (currently) the only setting
  printed to output filenames by default.
- `dataset_id::Int`: Dataset identifier. There should be a unique dataset ID for
  each set of observables.
- `cond_vintage::String`: Conditional data vintage, formatted `yymmdd`.
- `cond_id::Int`: Conditional dataset identifier. There should be a unique
  conditional dataset ID for each set of input, raw data mnemonics (not
  observables!).
- `cond_semi_names::Vector{Symbol}` and `cond_full_names::Vector{Symbol}`: names
  of observables for which we want to use semi- and full conditional data. All
  other observables are `NaN`ed out in the conditional data periods.
- `population_mnemonic::Nullable{Symbol}`: population series mnemonic in form
  `Nullable(:<mnemonic>__<source>)` (for example, `Nullable(:CNP16OV__FRED)`),
  or `Nullable{Symbol}()` if the model doesn't use population data

#### Dates

- `date_presample_start::Date`: Start date of pre-sample.
- `date_mainsample_start::Date`: Start date of main sample.
- `date_zlb_start::Date`: Start date of zero lower bound regime.
- `date_forecast_start::Date`: Start date of forecast period (or the period
  after the last period for which we have GDP data).
- `date_conditional_end::Date`: Last date for which we have conditional
  data. This is typically the same as `date_forecast_start` when we condition on
  nowcasts and current quarter financial data.

#### Estimation

- `reoptimize::Bool`: Whether to reoptimize the posterior mode. If `true` (the
    default), `estimate` begins reoptimizing from the model object's parameter
    vector.  See [Optimizing or Reoptimizing](@ref estimation-reoptimizing) for
    more details.
- `calculate_hessian::Bool`: Whether to compute the Hessian. If `true` (the
    default), `estimate` calculates the Hessian at the posterior mode.
- `n_mh_simulations::Int`: Number of draws from the posterior distribution per
  block.
- `n_mh_blocks::Int`: Number of blocks to run Metropolis-Hastings.
- `n_mh_burn::Int`: Number of blocks to discard as burn-in for Metropolis-Hastings.
- `mh_thin::Int`: Metropolis-Hastings thinning step.

#### Forecasting

- `forecast_jstep::Int`: Forecast thinning step.
- `forecast_block_size::Int`: Number of draws in each forecast block *before*
  thinning by `forecast_jstep`.
- `forecast_input_file_overrides::Dict{Symbol, String}`: Maps `input_type`(s) to
  the file name containing input draws for that type of forecast. See
  [Forecasting](@ref).
- `forecast_horizons::Int`: Number of periods to forecast.
- `impulse_response_horizons::Int`: Number of periods for which to calculate IRFs.

#### Alternative Policy

- `alternative_policy::AltPolicy`: See [Alternative Policies](@ref).

### Accessing Settings

The function `get_setting(m::AbstractModel, s::Symbol)` returns the value of the setting `s`
in `m.settings`. Some settings also have explicit getter methods that take only the model
object `m` as an argument. Note that not all are exported.

### Overwriting Default Settings

To overwrite default settings added during model construction, a user must
create a `Dict{Symbol, Setting}` and pass that into the model constructor as the
keyword argument `custom_settings`. If the `print`, `code`, and `description`
fields of the new `Setting` object are not provided, the fields of the existing
setting will be maintained. If new values for `print`, `code`, and `description`
are specified, and if these new values are distinct from the defaults for those
fields, the fields of the existing setting will be updated.

For example, overwriting `use_parallel_workers` should look like this:
```julia
custom_settings = Dict{Symbol, Setting}(
    :use_parallel_workers => Setting(:use_parallel_workers, true))
m = Model990(custom_settings = custom_settings)
```

By default, passing in `custom_settings` overwrites the entries in the model
object's `settings` field. However, with the additional keyword argument
`testing = true`, it will overwrite the entries in `test_settings`:

```julia
m = Model990(custom_settings = custom_settings, testing = true)
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
