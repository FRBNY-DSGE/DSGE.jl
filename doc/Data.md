# Data Description

The FRBNY DSGE Model takes an HDF5 file containing a matrix of data as
input. The columns of this matrix contain transformations of the
following series (the number corresponds to the column of data matrix):

1)  GDP (from the Bureau of Economic Analysis)
2)  Hours per capita (from the Bureau of Labor Statistics)
3)  Real wage growth (from the Bureau of Labor Statistics)
4)  GDP deflator (from the Bureau of Economic Analysis)
5)  Core PCE inflation (from the Bureau of Economic Analysis)
6)  Daily Nominal Federal Funds Rate (from the Federal Reserve Board)
7)  Personal consumption expenditure growth (from the Bureau of Economic Analysis)
8)  Investment growth (from the Bureau of Economic Analysis)
9)  Unemployment (from the Bureau of Labor Statistics)
10) Baa spreads (from Moody’s)
11) Long term expected CPI, (from the Survey of Professional Forecasters)
12) Total factor productivity (from the Federal reserve Bank of San Francisco)

The Total Factor Productivity series is made available by the Federal
Reserve Bank of San Francisco, and can be found
[here](http://www.frbsf.org/economic-research/total-factor-productivity-tfp/). The
Survey of Professional Forecasters series is made available by the
Federal Reserve Bank of Philadelphia, and can be found
[here](https://www.philadelphiafed.org/research-and-data/real-time-center/survey-of-professional-forecasters/historical-data/inflation-forecasts).

All other series are obtained via the Federal Reserve Economic Data (FRED) database.