# Data

## Description

The FRBNY DSGE Model takes an CSV file containing a matrix of data as
input. The columns of this file contain transformations of the
following series (the number corresponds to the column of data matrix):

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

1. The *Total Factor Productivity* series components are made available by the Federal Reserve
    Bank of San Francisco, and can be found
    [here](http://www.frbsf.org/economic-research/total-factor-productivity-tfp/) (series
    `alpha` and `dtfp` from the linked spreadsheet). Alternatively, they can be
    found as series `TFPJQ@USECON` (`alpha`) and `TFPKQ@USECON` (`dtfp`) via Haver Analytics. For more details on the series, see
    ```
    Fernald, John. "A quarterly, utilization-adjusted series on total factor productivity."
    Federal reserve bank of San Francisco working paper 19 (2012): 20912.
    ```

2. The *10-year Inflation Expectations* series from the *Survey of Professional
    Forecasters* is made available by the Federal Reserve Bank of Philadelphia, and
    can be found
    [here](https://www.philadelphiafed.org/research-and-data/real-time-center/survey-of-professional-forecasters/historical-data/inflation-forecasts)
    (series `INFCPI10YR` from the linked spreadsheet). Alternatively, it can be found as series
    `ASACX10@SURVEYS` via Haver Analytics.

3. The *10-year Treasury Yield* (zero-coupon, continuously compounded) series is made
    available by the Board of Governors of the Federal Reserve System, and can be found
    [here](http://www.federalreserve.gov/pubs/feds/2006/200628/200628abs.html) (series `SVENY10`
    from the linked spreadsheet). Alternatively, it can be found as series `FYCCZA@DAILY` via Haver
    Analytics. For more details on the series, see
    ```
    Gurkaynak, Refet S., Brian Sack, and Jonathan H. Wright. "The U.S. Treasury Yield Curve:
    1961 to the Present." Journal of monetary Economics 54.8 (2007): 2291-2304.
    ```

For additional details on the series, including mnemonics and transformations used,
please see Appendix A.I of
```
Del Negro, Marco, Marc P. Giannoni, and Frank Schorfheide. "Inflation in the
Great Recession and New Keynesian Models." American Economic Journal:
Macroeconomics 7.1 (2015): 168-196.
```

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
