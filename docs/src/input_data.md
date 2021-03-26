# [Input Data](@id input-data-step)

```@meta
CurrentModule = DSGE
```

Given all of the hard work put into specifying the model, one should be able to maintain
the input data painlessly. To that extent, DSGE.jl provides facilities to download
appropriate vintages of data series from FRED (Federal Reserve Economic Data).

Note that a sample input dataset for use with model `m990` is provided; see
[New York Fed Model 990 Data](@ref frbny-data) for more details. To update this
sample dataset for use with model `m990`, see [Update sample input data](@ref).

## Setup

To take advantage of the ability to automatically download data series from FRED via the
FredData.jl package, set up your FRED API access by following the directions
[here](https://github.com/micahjsmith/FredData.jl/blob/master/README.md).

## Loading data

At the most basic, loading data looks like this:

```julia
m = Model990()
df = load_data(m)
```

By default, `load_data` will look on the disk first to see if an appropriate vintage of data
is already present. If data on disk are not present, or if the data are invalid for any
reason, a fresh vintage will be downloaded from FRED and merged with the other data sources
specified. See [`load_data`](@ref) for more details.

The resulting DataFrame `df` contains all the required data series for this model, fully
transformed. The first row is given by the Setting `date_presample_start` and the last row
is given by `date_mainsample_end`. The first `n_presample_periods` rows of `df` are the
presample.

Driver functions including `estimate` accept this `df` as an argument and convert it into a
`Matrix` suitable for computations using `df_to_matrix`, which sorts the data, ensures the
full sample is present, discards the date column, and sorts the observable columns according
to the `observables` field of the model object.

## Non-FRED data sources

Some data series may not be available from FRED or one may simply wish to use a different
data source, for whatever reason. The data sources and series are specified in the
`input_series` field of an `Observable` object (see [ModelConstructors.jl](https://github.com/FRBNY-DSGE/ModelConstructors.jl)).
For each data source that is *not* `:fred`, a
well-formed CSV of the form `<source>_<yymmdd>.csv` is expected in the directory indicated
by `inpath(m, "raw")`.  For example, the following might be the contents of a data source
for two series `:series1` and `:series2`:

```
date,series1,series2
1959-06-30,1.0,NaN
1959-09-30,1.1,0.5
# etc.
```

Note that quarters are represented by the date of the *last* day of the quarter and missing
values are specified by `NaN`.

### Example

Let's consider an example dataset comprised of 10 macro series sourced from FRED and one
survey-based series sourced from, say, the Philadelphia Fed's [Survey of Professional
Forecasters](http://www.philadelphiafed.org/research-and-data/real-time-center/survey-of-professional-forecasters/historical-data/inflation-forecasts) via Haver Analytics. The `Observable` for that data series might look like this:

```
Observable(:obs_longcpi, [:ASAXC10__SPF], annualtoquarter, quartertoannual,
           "Median 10Y CPI Expectations", "Median 10Y CPI Expectations")
```

If the data vintage specified for the model is `151127` (Nov. 27, 2015), then the following
files are expected in `inpath(m, "raw")`:

```
spf_151127.csv
fred_151127.csv
```

The FRED series will be downloaded and the `fred_151127.csv` file will be automatically
generated, but the `spf_151127.csv` file must be manually compiled as shown above:

```
date,ASACX10
1991-12-31,4.0
# etc.
```

Now, suppose that we set the data vintage to `151222`, to incorporate the BEA's third
estimate of GDP. The `fred_151222.csv` file will be downloaded, but there are no updates to
the SPF dataset during this period. Regardless, the file `spf_151222.csv` must be present to
match the data vintage. The solution in this case is to manually copy and rename the older
SPF dataset. Although this is not an elegant approach, it is consistent with the concept of a
vintage as the data available at a certain point in time --- in this example, it just so
happens that the SPF data available on Nov. 27 and Dec. 22 are the same.

## Incorporate population forecasts

Many variables enter the model in per-capita terms. To that extent, we use data on
population levels to adjust aggregate variables into per-capita variables. Furthermore, we
apply the [Hodrick-Prescott filter](https://en.wikipedia.org/wiki/Hodrick%E2%80%93Prescott_filter)
("H-P filter") to the population levels to smooth cyclical components.

The user will ultimately want to produce forecasts of key variables such as GDP and then
represent these forecasts in standard terms. That is, one wants to report GDP forecasts in
aggregate terms, which is standard, rather than per-capita terms. To do this, we either
extrapolate from the last periods of population growth in the data, or use external
population forecasts.

Note that if external population forecasts are provided, non-forecast procedures, such as
model estimation, are also affected because the H-P filter smoothes back from the latest
observation.

To incorporate population forecasts,

1. Set the model setting `use_population_forecast` to `true`.
2. Provide a file `population_forecast_<yymmdd>.csv` to `inpath(m, "raw")`. Population
   forecasts should be in levels, and represent the same series as given by the
   `population_mnemonic` setting (defaults to `:CNP16OV`, or "Civilian Noninstitutional
   Population, Thousands"). If your population forecast is in growth rates, convert it to
   levels yourself. The first row of data should correspond to the last period of
   the main sample, such that growth rates can be computed. As many additional rows of
   forecasts as desired can be provided.

The file should look like this:
```
date,POPULATION
2015-12-31,250000
2016-03-31,251000
# etc.
```

## Dataset creation implementation details

Let's quickly walk through the steps DSGE.jl takes to create a suitable dataset.

First, a user provides a detailed specification of the data series and transformations used
for their model.

- the user specifies `m.observables`; the keys of this dictionary name
    the series to be used in estimating the model.

- the user specifies `m.observable_mappings`; the keys of this dictionary name observed variables, and the values correspond to the observable object, which contains information about the forward and reverse transforms as well as the input data series from which the observable is constructed.

- For a given observable, an input series, e.g.
    `m.observable_mappings[:obs_gdp].input_series`, is an array of mnemonics to be
    accessed from the data source listed after the mnemonic (separated by the double
    underscore). Note that these mnemonics do not correspond to observables one-to-one,
    but rather are usually series in *levels* that will be further transformed.

- There are also both forward and reverse transforms for a given observable,
    e.g. `m.observable_mappings[:obs_gdp].fwd_transform` and
    `m.observable_mappings[:obs_gdp].rev_transform`. The forward transform operates on a
    single argument, `levels`, which is a DataFrame of the data in levels returned by the
    function `load_data_levels`. The reverse transform operates on a forward transformed
    series (which is in model units) transforming it into human-readable units, such
    as one quarter percent changes or per-capita adjustments. Both transforms return a
    DataArray for a single series. These functions could do nothing, or they could
    perform a more complex transformation. See
    [Data Transforms and Utilities](@ref) for more information about series-specific
    transformations.

- the user adjusts data-related settings, such as `data_vintage`, `data_id`,
    `dataroot`, `date_presample_start`, `date_zlb_start`, `date_forecast_start`,
    and `use_population_forecast`. See [Working with Settings](@ref working-with-settings) for details.

Second, DSGE.jl attempts to construct the dataset given this setup through a call to
`load_data`. See [`load_data`](@ref) for more details.

- Intermediate data in levels are loaded. See [`load_data_levels`](@ref) for more details.
- Transformations are applied to the data in levels. See [`transform_data`](@ref) for more details.
- The data are saved to disk. See [`save_data`](@ref) for more details.

## Conditional data
The user can easily add conditional data for any observables. By "conditional data",
we mean that, in reality, some data has not become available yet, but we believe that
a certain number is a decent guess, so we want to forecast conditional on
our guessed data. For example, suppose we are in 2019:Q4, in which case we have not
observed 2019:Q4 GDP growth yet. However, we might have some idea of the number, so
we want our forecasts to be conditional on that guess.

To load such data, the user needs to include a "cond" folder within the input data folder,
i.e. this folder `joinpath(get_setting(m, :input_data), "cond")` should exist. Within this folder,
the user can create a csv file taking the form `cond_cdid=<xx>_cdvt=<yymmdd>.csv`.
The user should then make sure that the model object being used has the following settings

- `cond_id::Int64`: the conditional data's equivalent of `data_id` and will be inserted
  after the `cdid`. Note that the ID must be less than 100.

- `cond_vintage::String`: the conditional data's equivalent of `data_vintage` and will
  be inserted after the `cdvt`.

The contents of `cond_cdid=<xx>_cdvt=<yymmdd>.csv` should have columns for each raw data series
that is then used to construct a given conditional observable. The first column should be `date`
for the quarters of the conditional horizon, and the following columns should be for the raw data series.
For example, to obtain real GDP growth, we need to have a population forecast file with both CNP16OV and CE16OV,
the forecasted value of nominal GDP (under pnemonic GDP),
and the forecasted value of the GDP deflator (under pnemonic GDPDEF) since these series are all required to compute `obs_gdp`,
which is per-capita real GDP growth. For core inflation, we just need the index level for core PCE (under pnemonic PCEPILFE).

*Note that the csv should have only conditional horizon data*. If you have data for any historical quarters, then
the `DataFrame` with both historical and conditional data will not be created in REPL correctly. For example, if I am forecasting
2019:Q4 with a conditional forecast of 2019:Q4 values, then the data conditional csv should have only values
for 2019:Q4 (and onward). No values for 2019:Q3 or before should be in the conditional data csv.

Finally, to specify which variables should have conditional observations, make sure to set

- `cond_full_names::Vector{Symbol}`: variables when running a "full" conditional forecast.
  For Model 1002, this means averages of the current quarter's daily financial data
  as well as nowcasts of real GDP growth and core PCE inflation.

- `cond_semi_names::Vector{Symbol}`: variables when running a "semi" conditional forecast.
  For Model 1002, this means averages of the current quarter's daily financial data.

See the [default](https://github.com/FRBNY-DSGE/DSGE.jl/src/defaults.jl) settings
for an example of how these `cond_full_names` and `cond_semi_names` are initialized.

## Common pitfalls

Given the complexity of the data download, you may find that the dataset generated by
`load_data` is not exactly as you expect.
It is a good idea to compare the `observables.jl` file for your model with the [one
used by `Model1002`](https://github.com/FRBNY-DSGE/DSGE.jl/src/models/representative/m1002/observables.jl),
which uses all the features provided by the package for handling data.
Be certain that any significant differences are intentional. Here are also some common pitfalls to look out for:

- Ensure that the `data_vintage` and `cond_vintage` model settings are as you expect. (Try checking
    `data_vintage(m)` and `cond_vintage(m)`.)
- Ensure that the `data_id` and `cond_id` model settings are correct for the given model.
- Ensure that the `date_forecast_start` model setting is as you expect, and that is not
    logically incompatible with `data_vintage`.
- Ensure that the `date_conditional_end` model setting is as you expect, and that is not
    logically incompatible with `cond_vintage`.
- Double check the transformations specified in the `data_transforms` field of the model
    object.
- Ensure that the keys of the `observables` and `data_transforms` fields of the model object
    match.
- Check the input files for [Non-FRED data sources](@ref). They should be
    in the directory indicated by `inpath(m, "raw")`, be named appropriately given the
    vintage of data expected, and be formatted appropriately. One may have to copy and
    rename files of non-FRED data sources to match the specified vintage, even if the
    contents of the files would be identical.
- Look for any immediate issues in the final dataset saved
    (`data_dsid=<xx>_vint=<yymmdd>.csv`). If a data series in this file is all
    `NaN` values, then likely a non-FRED data source was not provided correctly.
- Ensure that the column names of the data CSV match the keys of the `observables` field of
    the model object.
- You may receive a warning that an input data file "does not contain the entire date range
    specified". This means that observations are not provided for some periods in which the
    model requires data. This is perfectly okay if your data series starts after
    `date_presample_start`.
- If you successfully created a data set but it is missing observations that you want to add,
  you may need to recreate the data set. By default, `load_data` checks if a data set
  with the correct vintage already exists. If it does, then `load_data` loads the saved data
  rather than recreate a data set from scratch. However, if the saved data set is missing
  observations, then you want to recreate it by calling `load_data(m; try_disk = false)`.
- If you have a column that is completely empty (all missing/NaN data), but you still want
  to load the data, then use the keyword `check_empty_columns = false`.

If you experience any problems using FredData.jl, ensure your API key is provided correctly
and that there are no issues with your firewall, etc. Any issues with FredData.jl proper
should be reported on that project's page.


## Update sample input data

A sample dataset is provided for the 2015 Nov 27 vintage. To update this dataset:

**Step 1**. See [Setup](@ref) to setup automatic data pulls using FredData.jl.

**Step 2**. Specify the exact data vintage desired:

```julia
julia>  m <= Setting(:data_vintage, "yymmdd")
```

**Step 3**. Create data files for the non-FRED data sources. For model `m990`, the required data files include
   `spf_<yymmdd>.csv` (with column `ASACX10`), `longrate_<yymmdd>.csv` (with
   column `FYCCZA`), and `fernald_<yymmdd>.csv` (with columns `TFPJQ` and
   `TFPKQ`). To include data on expected interest rates, the file
   `ois_<yymmdd>.csv` is also required. To include data on population forecasts,
   the file `population_forecst_<yymmdd>.csv` is also required (see
   [Incorporate population forecasts](@ref). See
   [New York Fed Model Input Data](@ref frbny-data) for details on the series
   used and links to data sources.

**Step 4**. Run `load_data(m)`; series from FRED will be downloaded and merged with the series from
   non-FRED data sources that you have already created. See [Common pitfalls](@ref) for some potential issues.


## Data Transforms and Utilities

```@autodocs
Modules = [DSGE]
Pages   = ["load_data.jl", "fred_data.jl", "transform_data.jl", "transformations.jl", "src/data/util.jl"]
Order   = [:function, :type]
```
