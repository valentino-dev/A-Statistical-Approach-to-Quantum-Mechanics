using DelimitedFiles, DataFrames, CSV, Plots

# path = joinpath("data", "data_fig5_check.csv\n");
# data = readdlm(path, ";", header=false);

data = Matrix(CSV.read("data/data_fig5_check.csv", DataFrame, header=false))
println(size(data[10,:]))
println(typeof(data[10,:]))
plot = histogram(vec(data))
savefig(plot, "plot.png")