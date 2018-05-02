using DSGE
using Base.Test
using DataFrames: DataFrame
using DataArrays: @data

# Previous and next quarter arithmetic
q = Date(1913,12,13)
DSGE.prev_quarter()
DSGE.next_quarter()
@test DSGE.prev_quarter(q) == Date(1913,09,30)
@test DSGE.next_quarter(q) == Date(1914,03,31)

start_date = Date(2000,01,01)
end_date   = Date(2010,01,01)
DSGE.get_quarter_ends(start_date,end_date)
DSGE.subtract_quarters(end_date, start_date)

df1 = DataFrame(a = [1])
df2 = DataFrame()
df1, df2 = DSGE.reconcile_column_names(df1, df2)
@test names(df1) == names(df2)
@test size(df2, 1) == 0

quartertodate("11q3")
quartertodate("1997q4")
quartertodate("1985-Q1")
@test_throws ParseError quartertodate("2005q9")
@test_throws ParseError quartertodate("12345")

df = DataFrame(date = ["1913-12-23", "1992-11-14", "2002-01-01", "2014-12-19"],
               x = @data([1, 2, NA, 4]))
DSGE.format_dates!(:date, df)

#DSGE.na2nan!(df)

# Test hpfilter, ensuring missing data are handled correctly
nanfront = [NaN; NaN]
nanback  = [NaN]
y  = [0.; 1.; 0.; 1.]
yn = [nanfront; y; nanback]
yf, yt   = hpfilter(y, 1600)
yfn, ytn = hpfilter(yn, 1600)
@test isequal(nanfront, yfn[1:length(nanfront)])
@test isequal(nanfront, ytn[1:length(nanfront)])
@test isequal(nanback, yfn[(end-length(nanback)+1):end])
@test isequal(nanback, ytn[(end-length(nanback)+1):end])

nothing
