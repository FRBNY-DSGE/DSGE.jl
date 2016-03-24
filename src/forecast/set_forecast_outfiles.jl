function set_forecast_outfiles(m, cond, output_type, forecast_settings)
    output_file_names = [
    "hist", # cond
    "hist_s", # cond
    "forecast", # cond
    "forecast_s", # cond
    "states",
    "datashocks", # shocks over uncond period
    "datashocks_ns", # non-standardized
    "shocks", # cond # shocks over period of conditional data
    "shocks_ns", # non-standardized
    "counter",
    "coutner_s",
    "shockdec",
    "shockdec_s",
    "ytrend",
    "ytrend_s",
    "dettrend", # cond
    "dettrend_s"] # cond
end
