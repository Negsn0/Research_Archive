using Plots
using Roots

# Binding Energy
function Binding(Z::Int, N::Int)
    A = Z + N
    av = 15.753
    as = 17.804
    aA = 23.69
    ac = 0.710
    a5 = 12.0
    
    val_v = av * A
    val_s = as * A^(2.0 / 3.0)
    val_a = aA * ((A - 2 * Z)^2) / A
    val_c = (ac * Z^2) / (A^(1.0 / 3.0))
    del = (N % 2 == 0) ? a5 / sqrt(A) : 0.0
    
    return val_v - val_s - val_a - val_c + del
end

# Functions
function f_lambda(λ, epsilons, Δ, n)
    sum = 0.0
    for ek in epsilons
        sum += (1 - (ek - λ) / sqrt((ek - λ)^2 + Δ^2))
    end
    return sum - n
end

function cal_G(epsilons, Δ, λ)
    sum = 0.0
    for ek in epsilons
        sum += 1.0 / sqrt((ek - λ)^2 + Δ^2)
    end
    return 2.0 / sum
end

function recal_D(epsilons, Δ, λ, G)
    sum = 0.0
    for ek in epsilons
        sum += Δ / sqrt((ek - λ)^2 + Δ^2)
    end
    return 0.5 * G * sum
end

function recal_N(λ, epsilons, Δ)
    sum = 0.0
    for ek in epsilons
        sum += 1 - ((ek - λ) / sqrt((ek - λ)^2 + Δ^2))
    end
    return sum
end

function get_parameters(N::Int)
    # N=100,102,...,132 に対応する5つの基準値を行ごとに格納
    table = [
        [-10.737, -9.231, -8.633, -7.941, -6.949],  # N=100
        [-10.605, -9.239, -8.522, -7.889, -6.888],  # N=102
        [-10.478, -9.246, -8.415, -7.838, -6.83 ],  # N=104
        [-10.365, -9.253, -8.313, -7.789, -6.776],  # N=106
        [-10.238, -9.26 , -8.214, -7.742, -6.725],  # N=108
        [-10.125, -9.266, -8.118, -7.696, -6.677],  # N=110
        [-10.016, -9.272, -8.026, -7.652, -6.632],  # N=112
        [-9.911 , -9.278, -7.938, -7.609, -6.59 ],  # N=114
        [-9.81  , -9.284, -7.852, -7.568, -6.55 ],  # N=116
        [-9.713 , -9.289, -7.77 , -7.528, -6.512],  # N=118
        [-9.619 , -9.294, -7.691, -7.49 , -6.477],  # N=120
        [-9.528 , -9.298, -7.614, -7.452, -6.444],  # N=122
        [-9.44  , -9.303, -7.54 , -7.416, -6.412],  # N=124
        [-9.356 , -9.307, -7.469, -7.382, -6.383],  # N=126
        [-9.311 , -9.274, -7.4  , -7.348, -6.355],  # N=128
        [-9.315 , -9.195, -7.333, -7.315, -6.329],  # N=130
        [-9.318 , -9.119, -7.284, -7.269, -6.304]   # N=132
    ]

    # (N-100) ÷ 2 の結果をテーブルの行インデックスに
    idx = div(N - 100, 2) + 1
    eps_vals = table[idx]

    # 各区画ごとに fill で繰り返し配列を作り、vcat で結合
    epsilons = vcat(
        fill(eps_vals[1], 6),
        fill(eps_vals[2], 8),
        fill(eps_vals[3], 2),
        fill(eps_vals[4], 4),
        fill(eps_vals[5], 12)
    )

    # 返り値は (N-100) と epsilons のタプル
    return (N - 100), epsilons
end

# Main Program
function main(A::Int)
    Z = 50
    N = A - Z
    Δ = 0.5 * (2 * Binding(Z, N) - Binding(Z, N - 1) - Binding(Z, N + 1))
    #println("Calculated Δ: ", Δ)

    n,epsilons=get_parameters(A)

    λ_sol = find_zero(
        λ -> f_lambda(λ, epsilons, Δ, n),
        (minimum(epsilons)-3.0, maximum(epsilons)+3.0)
    )
    #println("Calculated λ: ", λ_sol)

    G = cal_G(epsilons, Δ, λ_sol)
    #println("Calculated G: ", G)
    #println("Target n: ", n)
    #=
    open("116Sn_parameter.txt","w") do file
        println(file,"Calculated Δ: ", Δ)
        println(file,"Calculated λ: ", λ_sol)
        println(file,"Calculated G: ", G)
        println(file,"Recalculated Δ: ", reΔ)
        println(file,"Recalculated N: ", re_N)
        println(file,"Target n: ", n)
        println(file,"Defference N:",differN)
    end
    if abs(re_N - n) < 1e-6
        println("Particle number condition satisfied!")
    else
        println("Discrepancy in particle number condition.")
    end
    =#
    

    return λ_sol
end

x=Float64[]
y=Float64[]
#f(z)= 12/sqrt(z)

for A in 102:2:130
    push!(x,A)
    r=main(A)
    push!(y,r)
    #println("$(A) => $(r),")
end


#ΔvsA

scatter(x,y,
        xlabel="nucleon number A",ylabel="λ (MeV)",legend=false,label="Calculation")
#plot!(x,f,label="Δ=12/√A")
savefig("lambda_vs_A.pdf")


#GvsA
# g(x)=22/x
# # scatter(x,y,ylims=(0.0,0.3))
# plot!(x,g,label="g=22/A")
# savefig("G_vs_A.pdf")
# println