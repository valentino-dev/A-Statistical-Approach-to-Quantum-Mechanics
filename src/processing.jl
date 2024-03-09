using DelimitedFiles, DataFrames, CSV, Plots, Statistics,PGFPlotsX
using LaTeXStrings 
a = 0.5
m0 = 0.5
μ² = 2.0
λ = 0.0

pgfplotsx()

println("Loading data..")
data = Matrix(CSV.read("data/data_fig5_10k.csv", DataFrame, header=false))
println("Done loading.")



if true
    println("Ploting Action Decay..")
    start = 30
    sd_data = data[start:start+10, :]
    potential(x) = μ²*x .^2*0.5 .+ λ*x .^4
    action(x_0, x_1) = m0 .*(x_1.-x_0) .^2 ./a .+potential(x_0) .*a
    tot_action = sum(action(sd_data, circshift(sd_data, (0, -1))), dims=2)
    
    plot(eachindex(tot_action).+start, tot_action, label="Simulated  data")

    plot!(minorgrid=true, grid=true, gridwidth=1, gridalpha=0.4, gridstyle=:dot)
    title!("Action decay")
    xlabel!(L"t")
    ylabel!(L"S(t)")

    savefig("plots/ActionDecay.pdf")
    println("Done ploting.")
end

if false
    println("Plotting Probability Density")
    nbins = 201
    radius = 4
    plots = histogram(vec(data), bins=range(-radius, radius, length=nbins+1), normalize=:pdf, label="Simulated data")

    plot!(minorgrid=true, grid=true, gridwidth=1, gridalpha=0.4, gridstyle=:dot)
    title!("Probability density")
    xlabel!(L"x")
    ylabel!(L"P(x)")

    savefig(plots, "plots/ProbabilityDensity.pdf")
    println("Done ploting.")
end

if false
    println("Plotting correrlation function..")
    Dtau = 0.5
    tau = 1:Dtau:5
    t = tau/Dtau
    corr_data = data[50:1:end, :]
    correlation = zeros(size(t, 1), size(corr_data, 1))
    for i=eachindex(t)
        singleton = ones(size(corr_data, 1))
        mean!(singleton, corr_data.*circshift(corr_data, (0, t[i])))
        correlation[i, :] = singleton[:]
    end
    bined_correlation = ones(size(t, 1))
    mean!(bined_correlation, correlation)
    println(size)
    scatter(t*Dtau, bined_correlation, yaxis=:log, label="Simulated Data", marker=:x)

    plot!(minorgrid=true, grid=true, gridwidth=1, gridalpha=0.4, gridstyle=:dot)
    title!("Correlation function")
    xlabel!(L"\tau")
    ylabel!(L"\langle x(0)x(\tau)\rangle")

    savefig("plots/CorrelationFunction.pdf")
    println(size(tau))
    println(size(bined_correlation))
    println("Done ploting.")
end

if true
    println("Plotting correrlation function error..")
    Dtau = 0.5
    tau = 1:Dtau:5
    t = tau/Dtau
    corr_data = data[50:1:end, :]
    correlation = zeros(size(t, 1), size(corr_data, 1))
    for i=eachindex(t)
        singleton = ones(size(corr_data, 1))
        mean!(singleton, corr_data.*circshift(corr_data, (0, t[i])))
        correlation[i, :] = singleton[:]
    end
    bined_correlation = ones(size(t, 1))
    mean!(bined_correlation, correlation)
    println(size)
    scatter(t*Dtau, bined_correlation, yaxis=:log, label="Simulated Data", marker=:x)

    plot!(minorgrid=true, grid=true, gridwidth=1, gridalpha=0.4, gridstyle=:dot)
    title!("Correlation function error")
    xlabel!(L"\tau")
    ylabel!(L"\langle x(0)x(\tau)\rangle")

    savefig("plots/CorrelationFunctionError.pdf")
    println(size(tau))
    println(size(bined_correlation))
    println("Done ploting.")
end

function Gamma(t, _alpha, _beta, a_est, data)
    return mean((data[1:(size(data, 1)-t+1), _alpha].-a_est[_alpha]).*(data[t:end, _beta].-a_est[_beta]))
end

if false
    eff_data = data[1:10000, :]
    println("Plotting Gamma Method..")
    a_est = ones(size(eff_data, 2))
    mean!(a_est, eff_data')
    t = 1:size(eff_data, 1)
    gamma = zeros(size(t))
    _alpha = 1
    for i=t
        gamma[i] = Gamma(i, _alpha, _alpha, a_est, eff_data)
    end

    plot(eachindex(gamma[1:100]), gamma[1:100], label="Simulated data")
    plot!(minorgrid=true, grid=true, gridwidth=1, gridalpha=0.4, gridstyle=:dot)
    title!(L"\textrm{Gamma-Method: }\alpha=\beta=1")
    xlabel!(L"t")
    ylabel!(L"\Gamma_{\alpha\beta}(t)")
    savefig("plots/Gamma_Method_zoom.pdf")
end
if false
    eff_data = data[1:1000, :]
    println("Plotting Gamma Method..")
    a_est = ones(size(eff_data, 2))
    mean!(a_est, eff_data')
    relaxation_time = zeros(size(eff_data, 2))
    t = 1:size(eff_data, 1)
    gamma = zeros(size(t))


    for k=eachindex(relaxation_time)
        for i=t
            gamma[i] = Gamma(i, k, k, a_est, eff_data)
        end
        for i=range(1, size(gamma, 1)-1)
            if (gamma[i]>0 && gamma[i+1]<0)
                relaxation_time[k] = i
                break
            end
        end
    end
    rt_mean = mean(relaxation_time)
    histogram(relaxation_time, label="Simulated data")
    plot!(minorgrid=true, grid=true, gridwidth=1, gridalpha=0.4, gridstyle=:dot)
    title!(L"\textrm{Statistical independecy: mean}=23.74")
    xlabel!(L"t")
    ylabel!(L"\textrm{Frequency when }\Gamma(t)=0")
    savefig("plots/StatistcalIndependeny.pdf")
    println("Done plotting.")
end