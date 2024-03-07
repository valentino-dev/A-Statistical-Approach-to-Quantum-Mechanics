using DelimitedFiles, DataFrames, CSV, Plots, Statistics
a = 0.5
m0 = 0.5
μ² = 2.0
λ = 0.0

println("Loading data..")
data = Matrix(CSV.read("data/data_fig5.csv", DataFrame, header=false))
println("Done loading.")


if false
    println("Ploting Action Decay..")
    sd_data = data[30:40, :]
    potential(x) = μ²*x .^2*0.5 .+ λ*x .^4
    action(x_0, x_1) = m0 .*(x_1.-x_0) .^2 ./a .+potential(x_0) .*a
    tot_action = sum(action(sd_data, circshift(sd_data, (0, -1))), dims=2)
    plots = plot(eachindex(tot_action), tot_action)
    savefig(plots, "plots/ActionDecay.pdf")
    println("Done ploting.")
end

if false
    println("Plotting Probability Density")
    nbins = 201
    radius = 4
    plots = histogram(vec(data), bins=range(-radius, radius, length=nbins+1), normalize=:pdf)
    savefig(plots, "plots/ProbabilityDensity.pdf")
    println("Done ploting.")
end

if true
    println("Plotting correrlation function..")
    Δτ = 0.5
    τ = 1:Δτ:2.5
    t = τ/Δτ
    corr_data = data[36:end, :]
    correlation = zeros(size(t, 1), size(corr_data, 1))
    for i=eachindex(t)
        singleton = ones(size(corr_data, 1))
        mean!(singleton, corr_data.*circshift(corr_data, (0, t[i])))
        correlation[i, :] = singleton[:]
    end
    bined_correlation = ones(size(t, 1))
    mean!(bined_correlation, correlation)
    
    plots = scatter(t*Δτ, bined_correlation, yaxis=:log)
    savefig(plots, "plots/CorrelationFunction.pdf")
    println(τ, bined_correlation)
    println("Done ploting.")
end