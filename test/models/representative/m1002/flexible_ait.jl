# Test that regime-switch in 2020-Q3 to flexible AIT rule runs properly

using Test, ModelConstructors, DSGE, Dates, FileIO, Random

generate_fulldist_forecast_data = false
if VERSION < v"1.5"
    ver = "111"
else
    ver = "150"
end


m = Model1002("ss30")
