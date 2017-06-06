# [New York Fed DSGE Model 990 Data](@id frbny-data)

```@meta
CurrentModule = DSGE
```

## Data Series

The New York Fed DSGE Model takes an CSV file containing a matrix of data as
input. The columns of this file contain transformations of the following series
(the number corresponds to the column of data matrix):

1.  Output Growth (Bureau of Economic Analysis)
2.  Hours Worked (Bureau of Labor Statistics)
3.  Real Wage Growth (Bureau of Labor Statistics)
4.  Inflation (GDP Deflator) (Bureau of Economic Analysis)
5.  Inflation (Core PCE) (Bureau of Economic Analysis)
6.  Federal Funds Rate (Board of Governors of the Federal Reserve System)
7.  Consumption Growth (Bureau of Economic Analysis)
8.  Investment Growth (Bureau of Economic Analysis)
9.  Spread (Baa) (Board of Governors of the Federal Reserve System)
10. 10-year Inflation Expectations (Federal Reserve Bank of Philadelphia)
11. 10-year Interest Rate (Board of Governors of the Federal Reserve System)
12. Total Factor Productivity (Federal Reserve Bank of San Francisco)

The following series are used to transform some series into per capita terms:

1. Civilian Noninstitutional Population 16 Years and Over (Bureau of Labor Statistics)

Most data series used to construct the above are retrieved from
[FRED](https://research.stlouisfed.org/fred2) (Federal Reserve Bank of St.
Louis). Other data sources include:

- The *Total Factor Productivity* series components are made available by the Federal Reserve
    Bank of San Francisco, and can be found
    [here](http://www.frbsf.org/economic-research/total-factor-productivity-tfp/) (series
    `alpha` and `dtfp` from the linked spreadsheet). Alternatively, they can be
    found as series `TFPJQ@USECON` (`alpha`) and `TFPKQ@USECON` (`dtfp`) via Haver Analytics. For more details on the series, see

```
Fernald, John. "A Quarterly, Utilization-Adjusted Series on Total Factor Productivity."
Federal Reserve Bank of San Francisco Working Paper 19 (2012): 20912.
```

- The *10-year Inflation Expectations* series from the *Survey of Professional
    Forecasters* is made available by the Federal Reserve Bank of Philadelphia, and
    can be found
    [here](http://www.philadelphiafed.org/research-and-data/real-time-center/survey-of-professional-forecasters/historical-data/inflation-forecasts/)
    (series `INFCPI10YR` from the linked spreadsheet). Alternatively, it can be found as series
    `ASACX10@SURVEYS` via Haver Analytics.
- The *10-year Treasury Yield* (zero-coupon, continuously compounded) series is made
    available by the Board of Governors of the Federal Reserve System, and can be found
    [here](http://www.federalreserve.gov/pubs/feds/2006/200628/200628abs.html) (series `SVENY10`
    from the linked spreadsheet). Alternatively, it can be found as series `FYCCZA@DAILY` via Haver
    Analytics. For more details on the series, see

```
Gurkaynak, Refet S., Brian Sack, and Jonathan H. Wright. "The U.S. Treasury Yield Curve:
1961 to the Present." Journal of Monetary Economics 54.8 (2007): 2291-2304.
```

For additional details on the series, including mnemonics and transformations used,
please see Appendix A.I of

```
Del Negro, Marco, Marc P. Giannoni, and Frank Schorfheide. "Inflation in the
Great Recession and New Keynesian Models." American Economic Journal:
Macroeconomics 7.1 (2015): 168-196.
```

## Interest Rate Expectations Data

In our model (as used to compute the
forecasts referenced in Liberty Street Economics posts), we treat the zero lower bound by adding
*anticipated policy shocks* and data on the market-implied Federal
Funds rate path. We do this by giving the model the market-implied
Federal Funds rate path for the next `n_anticipated_shocks` quarters
and forcing the model's interest rate path to hit those values in
those quarters. Afterwards, the path is unconstrained. The model is
trained on data that includes six quarters of interest rate expectations. The user is
responsible for procuring interest rate expectations and appending it to the provided sample
data set, as discussed in this documentation.


### Implementation

If you are able to access data on the market-implied FFR path (or another form of interest
rate expectations), you can augment the sample dataset or your own dataset to enable the
anticipated policy shocks feature. We use internal data from the Federal Reserve Board on
the implied Federal Funds Rate derived from OIS quotes. (One could also use interest rate
expectations from Blue Chip Financial Forecasts or Survey of Professional Forecasters.)

**Step 1**. Choose a value for `n_anticipated_shocks` (we suggest `6`):

```julia
m <= Setting(:n_anticipated_shocks, 6, true, "nant", "Number of ant. pol. shocks")
```

**Step 2**. Add implied FFR data to the `data` matrix:

**2a.** Append `n_anticipated_shocks` columns of `NaN` values to the end of the
       `data` matrix.

**2b.** Construct a matrix of data, say `ImpliedFFR`, on anticipated policy
       shocks. Define

```pseudocode
For t from first quarter ZLB binds to last quarter ZLB binds
   For h from 1 quarter ahead to n_anticipated_shocks quarters ahead
       ImpliedFFR[t,h] := FFR at horizon h quarters ahead implied as of quarter t.
   End
End
```

**2c.** Fill in the `data` matrix with the `ImpliedFFR` matrix. The first
   row of the `ImpliedFFR` matrix should go in the row of the `data` matrix in
   which the ZLB first bound and the last row of the `ImpliedFFR` matrix should
   go in the row of the `data` matrix in which the ZLB last bound.

**Step 3**. With your updated input `data` matrix, the code will add the appropriate
  number of states, shocks, equilibrium condition equations, and measurement
  equations.

### Discussion

The implementation of anticipated policy shocks may not be immediately clear.
Consider the following made-up `data` matrix:

| t      | GDP  | FFR | Inf | ... | Spread | ImpFFR\_1 | ... | ImpFFR_H   |
| ------ | :--: | :-: | :-: | :-: | :----: | :-----:   | :-: | :--------: |
| 1960Q1 | 2.5  | 5.0 | 2.5 | ... | 1.5    | NaN       | ... | NaN        |
| 1960Q2 | 2.2  | 5.2 | 1.5 | ... | 1.3    | NaN       | ... | NaN        |
| ...    | ...  | ... | ... | ... | ...    | ...       | ... | ...        |
| 2008Q3 | 1.1  | 2.2 | 1.0 | ... | 1.5    | NaN       | ... | NaN        |
| 2008Q4 | -4.5 | 2.0 | 2.0 | ... | 1.3    | 1.0       | ... | 1.5        |
| ...    | ...  | ... | ... | ... | ...    | ...       | ... | ...        |
| 2013Q1 | 2.2  | 0.2 | 1.7 | ... | 1.7    | 0.2       | ... | 1.5        |
| 2013Q2 | 2.3  | 0.2 | 1.8 | ... | 1.6    | 0.2       | ... | 1.4        |

Interpret this as follows:

- For periods before 2008Q4, there was no forward guidance or ZLB to enforce,
  and we have no implied FFR values to enter.
- In 2008Q4, actual FFR (made-up) was 2.2. Market prices *implied* that markets
  expected an interest rate of 1.0 in 2009Q1 -- 1 period from now --
  and 1.5 `n_anticipated_shocks` periods from 2008Q4.
- In 2013Q2, actual FFR (made-up) was 0.2. Markets expected FFR to remain at 0.2
  in 2013Q3, ..., and expected FFR of 1.4 `n_anticipated_shocks` periods from
  2013Q2.

## References

For a more comprehensive treatment of anticipated policy shocks, see NY Fed
Staff Report
[The FRBNY DSGE Model](https://www.newyorkfed.org/medialibrary/media/research/staff_reports/sr647.pdf)
- page 12, for how the anticipated policy shocks are incorporated into the
  monetary policy rule,
- page 16, for how the anticipated policy shocks entered the log-linear
  equilibrium conditions,
- page 18, for how the anticipated policy shocks and data on market expectations
  enter the measurement equation,
- page 23, for how the anticipated policy shocks propagate through the model.

For more in depth discussion of anticipated policy shocks/forward guidance and
the impact on the macroeconomy, see NY Fed Working Paper
[The Forward Guidance Puzzle](https://www.newyorkfed.org/medialibrary/media/research/staff_reports/sr574.pdf)
by Marco Del Negro, Marc Giannoni, and Christina Patterson.

Thanks to [Matthew Cocci](https://github.com/MattCocci) for the *Discussion*.

## Disclaimer

The sample input data provided with *DSGE.jl* is made available for purposes of
demonstrating the function of the model only. By using the data you acknowledge
and agree to the following terms and conditions. If you do not agree to these
terms and conditions, do not use the data provided with the model.

Some data provided with *DSGE.jl* may be copyrighted by its owner, and
permission to use such copyrighted materials other than to demonstrate the
functioning of *DSGE.jl* for your personal use must be obtained from the owner.
The Federal Reserve Bank of New York cannot provide permission to use the data
other than as permitted in this agreement.

You may not use the name of the Federal Reserve Bank of New York to endorse or
promote products derived from the use of the data provided with *DSGE.jl*, nor
for any other commercial purpose.

The Federal Reserve Bank of New York does not guarantee the completeness or
accuracy of the data and does not provide updates or corrections to the data
provided with the model. By downloading and using the data, you acknowledge
and agree that your use of the data is at your own risk and that none of the
parties involved in creating, producing or delivering *DSGE.jl* is liable
for any loss, injury, claim, liability or damage of any kind resulting in any
way from: (a) any errors in or omissions from the data; (b) your use of the
data or any conclusions you draw from it, regardless of whether you received
any assistance from the Federal Reserve Bank of New York or its employees with
regard to the data; or (c) the files containing the data or use of the website
from which the data files were downloaded, including anything caused by any
viruses, bugs or malfunctions.

ALL DATA AND MATERIALS ARE PROVIDED ON AN "AS IS", "AS AVAILABLE" BASIS WITHOUT
WARRANTY OF ANY KIND. THE FEDERAL RESERVE BANK OF NEW YORK EXPRESSLY DISCLAIMS
ALL WARRANTIES EXPRESS AND IMPLIED, INCLUDING WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE, QUALITY AND NON-INFRINGEMENT. NEITHER
THE FEDERAL RESERVE BANK OF NEW YORK NOR ANY EMPLOYEE OR AFFILIATE SHALL BE
LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES OF ANY KIND
WHATSOEVER (INCLUDING, WITHOUT LIMITATION, DAMAGES FOR LOSS OF PROFITS,
BUSINESS INTERRUPTION, LOSS OF INFORMATION, OR ATTORNEYS' FEES) IN ANY WAY DUE
TO, RESULTING FROM OR ARISING IN CONNECTION WITH THE USE OR PERFORMANCE OF, OR
INABILITY TO USE DATA OR MATERIALS, WHETHER IN AN ACTION OF CONTRACT,
NEGLIGENCE OR OTHER TORTIOUS ACTION, AND REGARDLESS OF THE NEGLIGENCE OF THE
BANK OR ANY EMPLOYEE OR AFFILIATE, EVEN IF THE FEDERAL RESERVE BANK OF NEW YORK
HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

Reference to any specific commercial product, process or service does not
constitute or imply its endorsement, recommendation or favoring by the Federal
Reserve Bank of New York.

Company and product names mentioned in connection with the data remain the
trademark and property of their respective owners.

