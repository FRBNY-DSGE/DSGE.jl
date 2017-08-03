using HDF5, Plots

path = dirname(@__FILE__)

nmh10m4000 = h5read("$path/../../test/reference/lik_nmh=10_npart=4000.h5", "tpf")
nmh5m4000 = h5read("$path/../../test/reference/lik_nmh=5_npart=4000.h5","tpf")
nmh2m4000 = h5read("$path/../../test/reference/lik_nmh=2_npart=4000.h5","tpf")
kal = h5read("$path/../../test/reference/lik_nmh=2_npart=4000.h5", "kal")

plotly()
plot(kal, linewidth=2, label="Kalman")
plot!(nmh2m4000, label="TPF, N_{MH}=2, M=4000")
plot!(nmh5m4000, label="TPF, N_{MH}=5, M=4000")
plot!(nmh10m4000,label="TPF, N_{MH}=10, M=4000")
plot!(legend=:bottomright, xlabel="time", ylabel="log likelihood")
savefig("/data/dsge_data_dir/dsgejl/interns2017/plots/tuning_parameters.html")
gui()