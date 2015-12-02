# Anticipated Policy Shocks

## Description

We treat the zero lower bound in our model by adding *anticipated policy
shocks* and data on the market-implied Federal Funds rate path. We do this by
giving the model the market-implied Federal Funds rate path for the next
`n_anticipated_shocks` quarters and forcing the model's interest rate path
to hit those values in those quarters. Afterwards, the path is unconstrained.

## Implementation

If you are able to access data on the market-implied FFR path, you can augment
the sample dataset or your own dataset to enable the anticipated policy shocks
feature. We use internal data from the Federal Reserve Board on the implied
Federal Funds Rate derived from OIS quotes.

1. Choose a value for `n_anticipated_shocks` (we suggest `6`):

    ```julia
    m <= Setting(:n_anticipated_shocks, 6, true, "nant", "Number of ant. pol. shocks)
    ```

2. Add implied FFR data to the `data` matrix:
    1. Append `n_anticipated_shocks` columns of `NaN` values to the end of the
       `data` matrix.
    2. Construct a matrix of data, say `ImpliedFFR`, on anticipated policy
       shocks. Define

         ```pseudocode
         For t from first quarter ZLB binds to last quarter ZLB binds
            For h from 1 quarter ahead to n_anticipated_shocks quarters ahead
                ImpliedFFR[t,h] := FFR at horizon h quarters ahead implied as of quarter t.
            End
         End
         ```

    3. Fill in the `data` matrix with the `ImpliedFFR` matrix. The first
       row of the `ImpliedFFR` matrix should go in the row of the `data` matrix in
       which the ZLB first bound and the last row of the `ImpliedFFR` matrix should
       go in the row of the `data` matrix in which the ZLB last bound.
3. With your updated input `data` matrix, the code will add the appropriate
  number of states, shocks, equilibrium condition equations, and measurement
  equations.

## Discussion

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
