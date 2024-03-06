using DelimitedFiles, DataFrames, CSV, Plots
using Statistics
a = 0.5
m0 = 0.5
μ² = 2.0
λ = 0.0

data = Matrix(CSV.read("data/data_fig5.csv", DataFrame, header=false))

if false
    potential(x) = μ²*x .^2*0.5 .+ λ*x .^4
    action(x_0, x_1) = m0 .*(x_1.-x_0) .^2 ./a .+potential(x_0) .*a
    tot_action = sum(action(data, circshift(data, (0, -1))), dims=2)[1:50, :]
    plots = plot(eachindex(tot_action), tot_action)
    savefig(plots, "plots/ActionDecay.pdf")
end

if false
    nbins = 201
    radius = 4
    plots = histogram(vec(data), bins=range(-radius, radius, length=nbins+1), normalize=:pdf)
    savefig(plots, "plots/ProbabilityDensity.pdf")
end

if true
    Δτ = 0.5
    τ = 1:Δτ:2.5
    t = τ/Δτ
    correlation = ones(size(t, 1), size(data, 1))
    for i=eachindex(t)
        dataProd = data.*circshift(data, (0, t[i]))
        # println(dataProd)
        mean!(correlation[i, :], dataProd)
    end
    bined_correlation = ones(size(t, 1))
    println(correlation)
    mean!(bined_correlation, correlation)
    
    plots = plot(t, bined_correlation)
    savefig(plots, "plots/CorrelationFunction.pdf")
end