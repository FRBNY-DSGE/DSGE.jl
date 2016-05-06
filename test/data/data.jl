using Base.Test
using DataFrames: DataFrame, @data

DSGE.last_quarter_end()

start_date = Date(2000,01,01)
end_date = Date(2010,01,01)
DSGE.get_quarter_ends(start_date,end_date)
DSGE.subtract_quarters(end_date, start_date)

DSGE.stringstodates(["1913-12-23", "1992-11-14", "2002-01-01", "2014-12-19"])

quartertodate("11q3")
quartertodate("1997q4")
quartertodate("1985-Q1")
@test_throws ParseError quartertodate("2005q9")
@test_throws ParseError quartertodate("12345")

df = DataFrame(date = ["1913-12-23", "1992-11-14", "2002-01-01", "2014-12-19"],
               x = @data([1, 2, NA, 4]))
DSGE.format_dates!(:date, df)

#DSGE.na2nan!(df)

nothing
