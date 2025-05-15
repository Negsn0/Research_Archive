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
        [-10.876531601057266,-9.231054204559157,-8.801834344939756,-7.943138868718755,-6.948651483728945],  # N=100
        [-10.767494559580605,-9.238926343666868,-8.712820105692124,-7.8944179520266005,-6.8875121218095785],  # N=102
        [-10.664552106723697,-9.246376564418853,-8.629605617184374,-7.849643391064106,-6.830035411864607],  # N=104
        [-10.567306261475395,-9.253430984323195,-8.55180659193121,-7.808540978373863,-6.775982380899364],  # N=106
        [-10.47539048707391,-9.260113826698749,-8.479068740639846,-7.770857984179109,-6.725132989881736],  # N=108
        [-10.388466739588686,-9.26644758386822,-8.411065002312263,-7.736361207167576,-6.677284343588283],  # N=110
        [-10.306222838159016,-9.272453164001888,-8.34749307166367,-7.704835228757716,-6.632249097438617],  # N=112
        [-10.22837011700927,-9.278150023479341,-8.288073187467596,-7.676080846753962,-6.589854036635472],  # N=114
        [-10.154641324892388,-9.283556286397221,-8.23254615043161,-7.649913667488811,-6.549938806388283],  # N=116
        [-10.084788742297022,-9.2886888526454,-8.180671543443275,-7.626162838273004,-6.512354774922464],  # N=118
        [-10.018582490733394,-9.293563495797592,-8.132226130630258,-7.604669904307855,-6.476964013456936],  # N=120
        [-9.955809011802316,-9.29819495190931,-8.087002414752163,-7.58528777621315,-6.443638379441069],  # N=122
        [-9.896269696649775,-9.302597000185877,-8.044807335075511,-7.567879796048906,-6.412258691142284],  # N=124
        [-9.839779648887664,-9.30678253636723,-8.005461090136862,-7.552318891191568,-6.382713983211936],  # N=126
        [-9.786166566194119,-9.31076363957871,-7.968796071744052,-7.538486806712029,-6.354900834177016],  # N=128
        [-9.735269727639519,-9.31455163330972,-7.934655898236407,-7.526273408014326,-6.32872275793865],  # N=130
        [-9.686939075367633,-9.318157141107704,-7.902894536474321,-7.515576046461511,-6.304089652334923]   # N=132
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
    
    open("116Sn_parameter1.txt","w") do file
        println(file,"Calculated Δ: ", Δ)
        println(file,"Calculated λ: ", λ_sol)
        println(file,"Calculated G: ", G)
    end

    return nothing
end


main(116)


# x=Float64[]
# y=Float64[]
# #f(z)= 12/sqrt(z)

# for A in 102:2:130
#     push!(x,A)
#     r=main(A)
#     push!(y,r)
#     #println("$(A) => $(r),")
# end


#ΔvsA

# scatter(x,y,
#         xlabel="nucleon number A",ylabel="λ (MeV)",legend=false,label="Calculation")
# #plot!(x,f,label="Δ=12/√A")
# savefig("lambda_vs_A.pdf")


#GvsA
# g(x)=22/x
# # scatter(x,y,ylims=(0.0,0.3))
# plot!(x,g,label="g=22/A")
# savefig("G_vs_A.pdf")
# println