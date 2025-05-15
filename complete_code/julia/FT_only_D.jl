using Plots, Roots

# 定数
const k_B = 8.617333e-11  # Boltzmann定数 (MeV/K)

function f_lambda(λ, epsilons, Δ, n,kT)
    sum = 0.0
    for ek in epsilons
        sum += (1 -((1 - (2*fermi_dirac(ek,Δ,λ,kT)))* (ek - λ) / sqrt((ek - λ)^2 + Δ^2)))
    end
    return sum - n
end

function gap_equation(λ,epsilons,G,kT,Δ)
    sum = 0.0
    for ek in epsilons
        sum += 0.5*((1 - (2*fermi_dirac(ek,Δ,λ,kT)))/sqrt((ek - λ)^2 + Δ^2))
    end
    return sum - 1.0 / G 
end

function fermi_dirac(ek,Δ,λ,kT)
    if kT==0
        return 0
    else
        Ek=sqrt((ek-λ)^2 +Δ^2)
        return 1 / (1 + exp(Ek/kT))
    end
end

function main()
    #Δ=
    G = 0.116 #MeV
    lambda_T0 = -9.519 #MeV
    AAA=0.0
    n = 16
    epsilons = [
        -9.81, -9.81, -9.81, -9.81, -9.81, -9.81,
        -9.284, -9.284, -9.284, -9.284, -9.284, -9.284, -9.284, -9.284,
        -7.852, -7.852,
        -7.568, -7.568, -7.568, -7.568,
        -6.55, -6.55, -6.55, -6.55, -6.55, -6.55, -6.55, -6.55, -6.55, -6.55, -6.55, -6.55, 
    ]
    kT_range = 0.0:0.00001:0.75
    Δ_vals=Float64[]
    notified = false
    for kT in kT_range

        try
            Δ_solution = find_zero(
                Δ -> gap_equation(lambda_T0, epsilons, G,kT, Δ),
                (0.0, 2.0)  # bracket による探し方
            )
            push!(Δ_vals, Δ_solution)
        catch e
            if !notified
                AAA=kT
                println(kT)
                notified = true
            end
            push!(Δ_vals, 0.0)
        end
    end
    # 取得した Δ vs kT を Plot
    Δ0=Δ_vals[1]
    println(Δ0)
    println(AAA/k_B)
    plt=plot(
        kT_range/Δ0,
        Δ_vals/Δ0,
        xlabel = "kT/Δ0",
        ylabel = "Δ/Δ0",
        title  = "116Sn Δ/Δ0 vs kT/Δ0",
        legend = false,
        size = (600,600)
    )
end

main()