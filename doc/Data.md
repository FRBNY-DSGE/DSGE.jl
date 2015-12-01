# Data Description

The FRBNY DSGE Model takes an HDF5 file containing a matrix of data as
input. The columns of this matrix contain transformations of the
following series (the number corresponds to the column of data matrix):

1.  Output growth (Bureau of Economic Analysis)
2.  Hours worked (Bureau of Labor Statistics)
3.  Real wage growth (Bureau of Labor Statistics)
4.  Inflation (GDP Deflator) (Bureau of Economic Analysis)
5.  Inflation (Core PCE) (Bureau of Economic Analysis)
6.  FFR (Federal Reserve Board)
7.  Consumption growth (Bureau of Economic Analysis)
8.  Investment growth (Bureau of Economic Analysis)
10. Spread (Baa) (Moody’s)
11. 10y Inflation Expectations (Survey of Professional Forecasters)
12. Total Factor Productivity (Federal reserve Bank of San Francisco)

The following series were used to transform some series into per capita terms:

1. Civilian Noninstitutional Population 16 Years and Over (Bureau of Labor
Statistics)

The Total Factor Productivity series is made available by the Federal Reserve
Bank of San Francisco, and can be found
[here](http://www.frbsf.org/economic-research/total-factor-productivity-tfp/). 

The Survey of Professional Forecasters series is made available by the Federal
Reserve Bank of Philadelphia, and can be found
[here](https://www.philadelphiafed.org/research-and-data/real-time-center/survey-of-professional-forecasters/historical-data/inflation-forecasts).

All other series are retrieved from
[FRED](https://research.stlouisfed.org/fred2), Federal Reserve Bank of St.
Louis.

For more details on the series, including mnemonics and transformations used,
please see Appendix A.I of
```
Del Negro, Marco, Marc P. Giannoni, and Frank Schorfheide. "Inflation in the
Great Recession and New Keynesian Models." American Economic Journal:
Macroeconomics 7.1 (2015): 168-196.
```
