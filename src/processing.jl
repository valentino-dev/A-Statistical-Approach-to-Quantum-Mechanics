using DelimitedFiles, DataFrames, CSV, Plots
a = 0.5
m0 = 0.5
μ² = 2.0
λ = 0.0

data = Matrix(CSV.read("data/data_fig5_check.csv", DataFrame, header=false))
potential(x) = μ²*x .^2*0.5 .+ λ*x .^4
action(x_0, x_1) = m0 .*(x_1.-x_0) .^2 ./a .+potential(x_0) .*a
println(size(data)[1])
tot_action = zeros(size(data)[1])
println("got here")
for i in eachindex(tot_action)
    tot_action[i] = sum(action(data[i], circshift(data[i, :], 1)))
end
# print("tot_action", tot_action)
println(size(data))
# plots = histogram(vec(data), bins=range(-200, 200, length=22))
plots = plot([1, 2], [3, 4])
savefig(plots, "plot.png")