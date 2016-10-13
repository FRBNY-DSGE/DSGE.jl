# Input Data

Given all of the hard work put into specifying the model, one should be able to maintain
the input data painlessly. To that extent, *DSGE.jl* provides facilities to download
appropriate vintages of data series from FRED (Federal Reserve Economic Data).

Note that a sample input dataset for use with model `m990` is provided; see [Sample input
data](#sample-input-data) for more details. To update this sample dataset for use with
model `m990`, see [Update sample input data](#update-sample-input-data).

## Setup

To take advantage of the ability to automatically download data series from FRED via the
*FredData.jl* package, set up your FRED API access by following the directions
[here](https://github.com/micahjsmith/FredData.jl/blob/master/README.md).

## Loading data

At the most basic, loading data looks like this:

```
m = Model990()
df = load_data(m)
```

By default, `load_data` will look on the disk first to see if an appropriate vintage of data
is already present. If data on disk are not present, or if the data are invalid for any
reason, a fresh vintage will be downloaded from FRED and merged with the other data sources
specified. See `?load_data` for more details.

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
`data_series` field of the model object. For each data source that is *not* `:fred`, a
well-formed CSV of the form `<source>_<yymmdd>.csv` is expected in the directory indicated
by `inpath(m, "data")`.  For example, the following might be the contents of a data source
for two series `:series1` and `:series2`:

```
date,series1,series2
1959-06-30,1.0,NaN
1959-09-30,1.1,0.5
etc.
```

Note that quarters are represented by the date of the *last* day of the quarter and missing
values are specified by `NaN`.

### Example

Let's consider an example dataset comprised of 10 macro series sourced from FRED and one
survey-based series sourced from, say, the Philadelphia Fed's [Survey of Professional
Forecasters](https://www.philadelphiafed.org/research-and-data/real-time-center/survey-of-professional-forecasters/historical-data/inflation-forecasts)
via Haver Analytics:

```
julia> m.data_series
Dict{Symbol,Array{Symbol,1}} with 2 entries:
 :spf   => [:ASACX10]
 :fred  => [:GDP, :PCE, ...] #etc
```

If the data vintage specified for the model is `151127` (Nov. 27, 2015), then the following
files are expected in `inpath(m, "data")`:

```
spf_151127.csv
fred_151127.csv
```

The FRED series will be downloaded and the `fred_151127.csv` file will be automatically
generated, but the `spf_151127.csv` file must be manually compiled as shown above:

```
date,ASACX10
1991-12-31,4.0
etc.
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
2. Provide a file `population_forecast_<yymmdd>.csv` to `inpath(m, "data")`. Population
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
etc.
```

## Dataset creation implementation details

Let's quickly walk through the steps *DSGE.jl* takes to create a suitable dataset.

First, a user provides a detailed specification of the data series and transformations used
for their model.
- the user specifies `m.observables`; the keys of this dictionary name the series to be used
    in estimating the model.
- the user specifies `m.data_series`; the keys of this dictionary name data sources, and the
    values of this dictionary are lists of mnemonics to be accessed from that data source.
    Note that these mnemonics do not correspond to observables one-to-one, but rather are
    usually series in *levels* that will be further transformed.
- the user specifies `m.data_transforms`; the keys of this dictionary name the series to be 
    constructed and match the keys of `m.observables` exactly; the values of this dictionary
    are functions that operate on a single argument (`levels`) which is a DataFrame of the
    series specified in `m.data_series`. These functions return a DataArray for a single
    series. These functions could do nothing (e.g. return `levels[:, :SERIES1]`) or
    perform a more complex transformation, such as converting to one quarter percent changes
    or adjusting into per-capita terms.
- the user adjusts data-related settings, such as `data_vintage`, `dataroot`,
    `date_presample_start`, `date_mainsample_end`, and `date_zlbregime_start`, and
    `use_population_forecast`.

Second, *DSGE.jl* attempts to construct the dataset given this setup through a call to
`load_data`. See `?load_data` for more details.
- Intermediate data in levels are loaded. See `?load_data_levels` for more details.
- Transformations are applied to the data in levels. See `?transform_data` for more details.
- The data are saved to disk. See `?save_data` for more details.
    
## Common pitfalls

Given the complexity of the data download, you may find that the dataset generated by
`load_data` is not exactly as you expect. Here are some common pitfalls to look out for:
- Ensure that the `data_vintage` model setting is as you expect. (Try checking
    `data_vintage(m)`.)
- Ensure that the `date_mainsample_end` model setting is as you expect, and that is not
    logically incompatible with `data_vintage`.
- Ensure that the `data_series` field of the model object is set as expected.
- Double check the transformations specified in the `data_transforms` field of the model
    object.
- Ensure that the keys of the `observables` and `data_transforms` fields of the model object
    match.
- Check the input files for [non-FRED data sources](#non-fred-data-sources). They should be
    in the directory indicated by `inpath(m, "data")`, be named appropriately given the
    vintage of data expected, and be formatted appropriately. One may have to copy and
    rename files of non-FRED data sources to match the specified vintage, even if the
    contents of the files would be identical.
- Look for any immediate issues in the final dataset saved (`data_<yymmdd>.csv`). If a data
    series in this file is all `NaN` values, then likely a non-FRED data source was not
    provided correctly.
- Ensure that the column names of the data CSV match the keys of the `observables` field of
    the model object.
- You may receive a warning that an input data file "does not contain the entire date range
    specified". This means that observations are not provided for some periods in which the
    model requires data. This is perfectly okay if your data series starts after
    `date_presample_start`.

If you experience any problems using *FredData.jl*, ensure your API key is provided correctly
and that there are no issues with your firewall, etc. Any issues with *FredData.jl* proper
should be reported on that project's page.

## Sample input data

For more details on the sample input data provided -- which is used to estimate the provided
model `m990`, please see [Data](doc/Data.md).

For more details on using market interest rate expectations to treat the zero lower bound,
see [Anticipated Policy Shocks](doc/AnticipatedPolicyShocks.md). In particular, note that
our model, as used to compute the forecasts referenced in Liberty Street Economics posts,
is trained on data that includes six quarters of interest rate expectations. The user is
responsible for procuring interest rate expectations and appending it to the provided sample
data set, as discussed in the linked documentation here.

## Update sample input data

A sample dataset is provided for the 2015 Nov 27 vintage. To update this dataset:

1. See [above](#setup) to setup automatic data pulls using *FredData.jl*.
2. Specify the exact data vintage desired:

    ```
    julia> m <= Setting(:data_vintage, "yymmdd")
    ```
3. Create data files for the non-FRED data sources (specified in `m.data_series`). For model
   `m990`, the required data files include `spf_<yymmdd>.csv` (with column `ASACX10`),
   `longrate_<yymmdd>.csv` (with column `FYCCZA`), and `fernald_<yymmdd>.csv` (with columns
   `TFPJQ` and `TFPKQ`). To include data on expected interest rates, the file
   `ois_<yymmdd>.csv` is also required. To include [data on population
   forecasts](#incorporate-population-forecasts), the file `population_forecst_<yymmdd>.csv`
   is also required. See [Data](doc/Data.md) for details on the series used and links to
   data sources.
4. Run `load_data(m)`; series from FRED will be downloaded and merged with the series from
   non-FRED data sources that you have already created. See [Common
   pitfalls](#common-pitfalls) for some potential issues.
