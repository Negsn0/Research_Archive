#熱力学的な量を再計算するために用いるコード
using Plots, Roots, Statistics
using NLsolve  # 2変数(λ, Δ)を同時に解くために使用

# --- 定数 ---
const k_B = 8.617333e-11  # Boltzmann定数 (MeV/K)
# --- フェルミ分布関数 ---
function fermi_dirac(ek, Δ, λ, kT)
    # 温度が極めて低い場合の閾値（必要に応じて調整してください）
    kT_threshold = 1e-10
    if kT < kT_threshold
        return 0.0
    else
        Ek = sqrt((ek - λ)^2 + Δ^2)
        return 1.0 / (1.0 + exp(Ek / kT))
    end
end

# --- 粒子数固定のための方程式 ---
function f_lambda(λ, epsilons, Δ, n, kT)
    sum_val = 0.0
    for ek in epsilons
        Ek = sqrt((ek - λ)^2 + Δ^2)
            ratio = (ek - λ) / Ek
        # フェルミ・ディラック分布の値
        f_val = fermi_dirac(ek, Δ, λ, kT)
        # 各状態の寄与： 1 - (1 - 2f_val) * (ek-λ)/Eₖ
        sum_val += 1.0 - (1.0 - 2.0 * f_val) * ratio
    end
    return sum_val - n
end
# --- ギャップ方程式 ---
function gap_equation(λ, epsilons, G, kT, Δ)
    sum_val = 0.0
    for ek in epsilons
        Ek = sqrt((ek - λ)^2 + Δ^2)
        term = 0.5 * ((1.0 - (2.0 * fermi_dirac(ek, Δ, λ, kT))) / Ek)
        sum_val += term
    end
    return sum_val - 1.0 / G
end
# --- 2変数(λ, Δ)を同時に解くためのシステム ---
"""
    system_equations!(F, x, epsilons, n, G, kT)

    x[1] = λ
    x[2] = Δ

    F[1] = f_lambda(λ, epsilons, Δ, n, kT)
    F[2] = gap_equation(λ, epsilons, G, kT, Δ)
"""
function system_equations!(F, x, epsilons, n, G, kT)
    λ = x[1]
    Δ = x[2]
    F[1] = f_lambda(λ, epsilons, Δ, n, kT)
    F[2] = gap_equation(λ, epsilons, G, kT, Δ)
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

chem_pt = Dict(
    102 => -11.913765216101263,
    104 => -10.92796941426944,
    106 => -10.302637575488095,
    108 => -9.79850424368181,
    110 => -9.36880138797778,
    112 => -8.967695718459005,
    114 => -8.552386157716848,
    116 => -8.115050478710664,
    118 => -7.682093164254721,
    120 => -7.270809583982142,
    122 => -6.887453679322329,
    124 => -6.529781623549784,
    126 => -6.171020399074939,
    128 => -5.745012991641002,
    130 => -5.022406554420971
)

G_dict = Dict(
    102 => 0.26983794136816835,
    104 => 0.1991503783688934,
    106 => 0.16814905755002513,
    108 => 0.14943144839921188,
    110 => 0.13680183990072334,
    112 => 0.12767776052896732,
    114 => 0.12097959949668663,
    116 => 0.11614984979358076,
    118 => 0.11283525855960676,
    120 => 0.1107590171328641,
    122 => 0.10969483769997944,
    124 => 0.10936039319714089,
    126 => 0.10928457736919238,
    128 => 0.1087962209347509,
    130 => 0.10742872699345221
)

function entropy(epsilons,Δ,λ,kT)
    sum = 0.0
    for ek in epsilons
        f=fermi_dirac(ek, Δ, λ, kT)
        sum += f*log(f) + (1-f)*log(1-f)
    end
    return - sum
end

function vk2(eps,Δ,λ)
    return 0.5*(1 - (eps-λ)/sqrt((eps-λ)^2 + Δ^2))
end

function vk2_T(eps, Δ, λ, T)
    Ek = sqrt((eps-λ)^2 + Δ^2)
    return 0.5 * (1 - ((eps-λ)/Ek) * tanh(Ek/(2*T)))
end

function Hamiltonian_T(epsilons, Δ, λ, G, T)
    sum = 0.0
    for ek in epsilons
        sum += vk2_T(ek, Δ, λ, T) * (ek-λ)
    end
    return sum - Δ^2/G
end

function Hamiltonian(epsilons, Δ, λ, G)
    sum=0.0
    for ek in epsilons
        sum +=vk2(ek,Δ,λ)*(ek-λ )
    end
    return sum - Δ^2/G
end

function insert_breaks(x::Vector, y::Vector, threshold::Real)
    newx = Float64[]
    newy = Float64[]
    
    # 最初の点を追加
    push!(newx, x[1])
    push!(newy, y[1])
    
    # 隣接する点間の差をチェック
    for i in 1:length(x)-1
        if abs(y[i+1] - y[i]) > threshold
            # 差がしきい値を超えたら、NaN を挿入して線を切る
            push!(newx, NaN)
            push!(newy, NaN)
        end
        push!(newx, x[i+1])
        push!(newy, y[i+1])
    end
    return newx, newy
end

function Free_E(H,λ,n,kT,S)
    return H - λ*n -kT*S
end

function main(N::Int)
    # --- パラメータ設定 ---
    G = 0.11924213367895861      # (MeV)
    #G = get(G_dict,N,nothing)
    lambda_init = get(chem_pt,N,nothing) # (MeV) : 化学ポテンシャルの初期値
    d_init = 1.0
    n,epsilons = get_parameters(N)
    
    # 解を格納する配列
    λ_vals = Float64[]   # 化学ポテンシャル
    Δ_vals = Float64[]   # ギャップ
    S_vals = Float64[]   # エントロピー
    F_vals = Float64[]
    HT_vals = Float64[]
    # 「0以外の解が見つからなかった」タイミングを記録
    no_solution_kT = NaN
    # 収束判定で「これは 0 とみなす」際のしきい値 (適宜調整)
    delta_zero_threshold = 5e-3
    kT_vals = 0.0:0.00001:1.0

    # 相転移点に到達したかどうかのフラグ
    phase_transition_reached = false

    # --- ループして (λ, Δ) を同時に求める ---
    for this_kT in kT_vals
        if !phase_transition_reached
            try
                # (λ, Δ) を同時に解く
                sol = nlsolve(
                    (F, x) -> system_equations!(F, x, epsilons, n, G, this_kT),
                    [lambda_init, d_init],  # 初期値
                    xtol = 1e-10, ftol = 1e-10, iterations = 1000
                )
                
                λ_sol = sol.zero[1]
                Δ_sol = sol.zero[2]

                # もしギャップが閾値未満なら、明示的にΔ = 0.0に指定し、相転移フラグを立てる
                if abs(Δ_sol) < delta_zero_threshold
                    Δ_sol = 0.0
                    phase_transition_reached = true
                    no_solution_kT = this_kT
                end
                
                S = entropy(epsilons, Δ_sol, λ_sol, this_kT)
                Ht= Hamiltonian_T(epsilons,Δ_sol, λ_sol, G,this_kT)
                F = Free_E(Ht,λ_sol,n,this_kT,S)
                push!(F_vals,F)
                push!(HT_vals,Ht)
                push!(λ_vals, λ_sol)
                push!(Δ_vals, Δ_sol)
                push!(S_vals, S)
                lambda_init = λ_sol  # 次の初期値として更新
            
            catch e
                # 例外処理：解が見つからなかった場合は、Δ = 0.0 とし、相転移フラグを立てる
                FL(λ) = f_lambda(λ, epsilons, 0.0, n, this_kT)
                l_min = lambda_init - 1.0
                l_max = lambda_init + 1.0
                λ_sol = find_zero(FL, (l_min, l_max))
                Δ_sol = 0.0
                S = entropy(epsilons, Δ_sol, λ_sol, this_kT)
                Ht= Hamiltonian_T(epsilons,Δ_sol, λ_sol, G,this_kT)
                F = Free_E(Ht,λ_sol,n,this_kT,S)
                push!(F_vals,F)
                push!(HT_vals,Ht)
                push!(λ_vals, λ_sol)
                push!(Δ_vals, Δ_sol)
                push!(S_vals, S)
                lambda_init = λ_sol
                phase_transition_reached = true
            end
        else
            # すでに相転移点（Δが閾値以下）に達している場合、
            # 以降は常に Δ = 0.0 として、例外処理（またはΔ=0専用の枝）を実行する
            F(λ) = f_lambda(λ, epsilons, 0.0, n, this_kT)
            l_min = lambda_init - 1.0
            l_max = lambda_init + 1.0
            λ_sol = find_zero(F, (l_min, l_max))
            Δ_sol = 0.0
            S = entropy(epsilons, Δ_sol, λ_sol, this_kT)
            Ht= Hamiltonian_T(epsilons,Δ_sol, λ_sol, G,this_kT)
            F = Free_E(Ht,λ_sol,n,this_kT,S)
            push!(F_vals,F)
            push!(HT_vals,Ht)
            push!(λ_vals, λ_sol)
            push!(Δ_vals, Δ_sol)
            push!(S_vals, S)
            lambda_init = λ_sol
        end
    end

    # 最初に得られたギャップを参照ギャップ(Δ0)とする
    #   (もし最初が 0 に限りなく近かったら別途対処が必要)
    Δ0 = Δ_vals[1]  
    
    # --- 数値微分による比熱の計算 ---
    # kT_vals は Range 型なので、配列に変換してから微分を計算
    T_array = collect(kT_vals)
    # 隣接する内部エネルギーの差分 dH と温度差 dT を計算
    dH = diff(HT_vals)
    dT = diff(T_array)
    # 比熱 C = dH/dT の計算（単位は H の単位あたりの MeV/MeV）
    C_vals = dH ./ dT
    # 差分は区間の中間点に対応するので、中間温度を計算
    T_mid = (T_array[1:end-1] .+ T_array[2:end]) ./ 2

    # --- 負の比熱値を除去するフィルタ ---
    # ここでは、C_vals が正の値となる要素のみ、その対応する T_mid の値も抽出する
    positive_indices = findall(c -> c > 0, C_vals)
    T_mid_positive = T_mid[positive_indices]
    C_vals_positive = C_vals[positive_indices]
    
    threshold = 0.01
    T_mod,C_mod = insert_breaks(T_mid_positive,C_vals_positive,threshold)
    d0=Δ_vals[2]
    #return kT_vals,Δ_vals,HT_vals,T_mod,C_mod,d0,no_solution_kT
    return S_vals,F_vals,kT_vals,d0
end

S,F,kT,d0=main(116)
plot(kT/d0,F,ylabel="F(MeV)",xlabel="kT/Δ₀",legend=false)
savefig("F_kT.pdf")
# T_vals=Float64[]
# for A in 102:2:130
#     #a,b,c,d,e.f=main(A)
#     a,d=main(A)
#     push!(T_vals,a/d)
# end
# A=102:2:130

# plot(xlabel="nucleon number A",ylabel="kTc/Δ₀",xlims=(101,131),ylims=(0.0,1.0))
# scatter!(A,T_vals,label="Calculation")
# hline!([0.55],label="kTc/Δ₀=0.55")
# savefig("kTc_vs_A.pdf")
#=p1=plot(xlabel="kT/Δ₀",ylabel="Δ/Δ₀",xlims=(0.0,1.0),xticks=(0.0:0.1:1.0),
        ylims=(0.0,1.0))

p2=plot(xlabel="kT/Δ₀",ylabel="<H'> (MeV)",xlims=(0.0,1.0),xticks=(0.0:0.1:1.0))

p3=plot(xlabel="kT/Δ₀",ylabel="C",xlims=(0.0,1.0),xticks=(0.0:0.1:1.0))

#p1はD/D0 vs kT/D0,p2はE vs kT/D0,p3はC vs kT/D0
for A in 102:2:110
    kT,d,HT,T,C,d0=main(A)
    plot!(p1,kT/d0,d/d0,label="$(A)Sn")
    plot!(p2,kT/d0,HT/d0,label="$(A)Sn")
    plot!(p3,T/d0,C,label="$(A)Sn")
end

savefig(p1,"102d.pdf")
savefig(p2,"102H.pdf")
savefig(p3,"102C.pdf")

p1=plot(xlabel="kT/Δ₀",ylabel="Δ/Δ₀",xlims=(0.0,1.0),xticks=(0.0:0.1:1.0),
        ylims=(0.0,1.0))

p2=plot(xlabel="kT/Δ₀",ylabel="<H'> (MeV)",xlims=(0.0,1.0),xticks=(0.0:0.1:1.0))

p3=plot(xlabel="kT/Δ₀",ylabel="C",xlims=(0.0,1.0),xticks=(0.0:0.1:1.0))

#p1はD/D0 vs kT/D0,p2はE vs kT/D0,p3はC vs kT/D0
for A in 112:2:120
    kT,d,HT,T,C,d0=main(A)
    plot!(p1,kT/d0,d/d0,label="$(A)Sn")
    plot!(p2,kT/d0,HT/d0,label="$(A)Sn")
    plot!(p3,T/d0,C,label="$(A)Sn")
end

savefig(p1,"112d.pdf")
savefig(p2,"112H.pdf")
savefig(p3,"112C.pdf")

p1=plot(xlabel="kT/Δ₀",ylabel="Δ/Δ₀",xlims=(0.0,1.0),xticks=(0.0:0.1:1.0),
        ylims=(0.0,1.0))

p2=plot(xlabel="kT/Δ₀",ylabel="<H'> (MeV)",xlims=(0.0,1.0),xticks=(0.0:0.1:1.0))

p3=plot(xlabel="kT/Δ₀",ylabel="C",xlims=(0.0,1.0),xticks=(0.0:0.1:1.0))

#p1はD/D0 vs kT/D0,p2はE vs kT/D0,p3はC vs kT/D0
for A in 122:2:130
    kT,d,HT,T,C,d0=main(A)
    plot!(p1,kT/d0,d/d0,label="$(A)Sn")
    plot!(p2,kT/d0,HT/d0,label="$(A)Sn")
    plot!(p3,T/d0,C,label="$(A)Sn")
end

savefig(p1,"122d.pdf")
savefig(p2,"122H.pdf")
savefig(p3,"122C.pdf")

=#
#p=plot(xlabel="kT",ylabel="Δ",ylims=(0.0,1.5),xlims=(0.0,1.0),xticks=(0.0:0.1:1.0),title="Δ vs kT")
# nokTΔ₀_vals= Float64[]
# HT_vals  = Float64[]
# ΔΔ₀_vals   = Float64[]


# DD=Float64[]

# for A in 102:2:130
#     kT,Δ,HT,T,C,d0=main(A)
#     #plot!(kT/D1,H,label="$(A)Sn")
#     # plot!(kT,Δ,label="$(A)Sn")
#     push!(DD,d0)
# end
# A=52:2:80
# p=plot(ylims=(0.0,2.0))
# scatter!(A,DD)
# f(x)=12/sqrt(x)
# AA=102:2:130
# F=f.(AA)
# plot!(A,F)
# p=plot(xlims=(102,130),ylims=(0.0,1.0),xlabel="nucleon number A",
#     ylabel="kTc/Δ₀",title="Phase Transition Temperature vs. Nucleon Number")
# A=102:2:130
# scatter!(A,nokTΔ₀_vals,label="Transition Temp.")
# k=mean(nokTΔ₀_vals)
# hline!([k],label="Avg. Transition Temp.")
# println(k)
# display(p)
# savefig("A_kT.pdf")
# p=plot(xlabel="kT",ylabel="C",xlims=(0.0,1.0),xticks=(0.0:0.1:1.0),title="C vs kT")
# for A in 102:2:130
#     p=plot(xlabel="kT",ylabel="C",xlims=(0.0,1.0),xticks=(0.0:0.1:1.0),title="C vs kT")
#     kT,Δ,HT,T,C,nokT=main(A)
#     plot!(T,C,label="$(A)Sn")
#     display(p)
# end


# kT,Δ,HT,T,C,nokT=main(102)
# plot!(T,C,label="102Sn")
# display(p)
#savefig("Comp_FT_C.pdf")
# p=plot(xlabel="kT (MeV)",ylabel="Entropy S",title="Entropy vs kT",xlims=(0.0,1.0))

# kT,Δ,HT,T,C,a,F,S=main(116)
# plot!(kT,S,label="116Sn")
# display(p)

# savefig("S_vs_kT.pdf")
