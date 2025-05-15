using LinearAlgebra
using Plots

#jシェルの情報。構造体で定義する。
struct Shell
    j::Float64  #角運動量
    Omega::Int  #最大ペア数　j+1/2
    epsilon::Float64    #1粒子エネルギー
    N::Int  #充填数 2j+1個
end

#単粒子エネルギーのDict
epsilon_dict = Dict(
    100 => [-10.737, -9.231, -8.633, -7.941, -6.949],
    102 => [-10.605, -9.239, -8.522, -7.889, -6.888],
    104 => [-10.478, -9.246, -8.415, -7.838, -6.83],
    106 => [-10.365, -9.253, -8.313, -7.789, -6.776],
    108 => [-10.238, -9.26, -8.214, -7.742, -6.725],
    110 => [-10.125, -9.266, -8.118, -7.696, -6.677],
    112 => [-10.016, -9.272, -8.026, -7.652, -6.632],
    114 => [-9.911, -9.278, -7.938, -7.609, -6.59],
    116 => [-9.81, -9.284, -7.852, -7.568,-6.55],
    118 => [-9.713, -9.289, -7.77, -7.528, -6.512],
    120 => [-9.619, -9.294, -7.691, -7.49, -6.477],
    122 => [-9.528, -9.298, -7.614, -7.452, -6.444],
    124 => [-9.44, -9.303, -7.54, -7.416, -6.412],
    126 => [-9.356, -9.307, -7.469, -7.382, -6.383],
    128 => [-9.311, -9.274, -7.4, -7.348, -6.355],
    130 => [-9.315, -9.195, -7.333, -7.315, -6.329],
    132 => [-9.318, -9.119, -7.284, -7.269, -6.304]
)

#化学ポテンシャルのDict
chem_pt = Dict(
    102 => -12.898430409469947,
    104 => -11.678420563207569,
    106 => -11.058357158681536,
    108 => -10.614875054057995,
    110 => -10.269808853039285,
    112 => -9.981950183364672,
    114 => -9.733204164094236,
    116 => -9.510885290039697,
    118 => -9.303827472138105,
    120 => -9.10260473084743,
    122 => -8.897400025937786,
    124 => -8.680843801569413,
    126 => -8.446484703344876,
    128 => -8.191630416674602,
    130 => -7.9290544835707335,
    132 => -7.675286007457947
)

#seniority の取りうる値
function allowed_seniority(n::Int)
    start = iseven(n) ? 0 : 1
    return collect(start:2:n)
end

#核子から各シェルの中性子数N_jを決定する関数
function generate_shells(A::Int)
    Z = 50  #陽子数
    N = A - Z   #中性子数
    n = N - 50 #魔法数50から溢れた中性子数
    epsilons = get(epsilon_dict,A,nothing)
    if epsilons === nothing
        error("入力された核子数 A = $A に対応するパラメータが定義されていません")
    end
    #各shellのj値
    j_vals = [5/2, 7/2, 1/2, 3/2 ,11/2]
    shells=Shell[]
    for (i,eps) in enumerate(epsilons)
        j=j_vals[i]
        Omega = Int((2*j+1)/2)
        N_shell=min(2*Omega,n)
        push!(shells,Shell(j,Omega,eps,N_shell))
        n -=N_shell
        if n<=0
            break
        end
    end
    return shells
end
"
各shellについて、seniorityの値のリストを求め、
それの直積によって全シェルの組み合わせを生成。
"
function generate_configurations(shells::Vector{Shell})
    s_options = [allowed_seniority(shell.N) for shell in shells]
    configs = collect(Iterators.product(s_options...))
    #各configはタプルとして得られる。
    return configs
end

"
与えられた全核状態に対してエネルギーを計算
"
function configuration_energy(A,shells::Vector{Shell},s_config::NTuple{N,Int},g::Float64) where N
    E = 0.0
    for i in 1:length(shells)
        shell=shells[i]
        s=s_config[i]
        λ = get(chem_pt,A,nothing)
        # 各シェルのエネルギー寄与：1粒子項 + ペアリング項
        E += (shell.epsilon - λ)* shell.N -
            (g / 4) * (s^2 - 2*s*(shell.Omega + 1) + 2*shell.N*(shell.Omega + 1) - shell.N^2)
    end
    return E
end

function main(A::Int)
    g =0.11615 #ペアリング強度
    shells = generate_shells(A)
    configs = generate_configurations(shells)
    

    # 各全核状態（すなわち各シェルの seniority 値の組み合わせ）ごとにエネルギーを計算
    energies = [configuration_energy(A,shells, config, g) for config in configs]

    # ここからは、分配関数・平均エネルギーなどの熱力学量の計算例
    T_range = 1.e-6:0.0001:1.0
    Cv_vals = Float64[]
    E_vals = Float64[]
    F_vals = Float64[]
    S_vals = Float64[]
    for T in T_range
        beta = 1 / T
        boltz_factors = exp.(-beta .* energies)
        Z = sum(boltz_factors)
        p = boltz_factors ./ Z
        E_avg = sum(energies .* p)
        Cv=#数値微分を行って代入したい
        S = -sum(p .* log.(p))
        F = E_avg - T * S
        
        push!(E_vals, E_avg)
        push!(F_vals, F)
        push!(S_vals,S)
    end
    
    # 数値微分のための温度刻み幅（T_rangeは等間隔とする）
    deltaT = Float64(T_range.step)  # 例では 0.0001

    # 数値微分（中央差分）を用いて dE/dT を計算し、それを Cv とする
    Cv_vals = Float64[]
    nT = length(E_vals)
    for i in 1:nT
        if i == 1
            # 最初の点：前進差分
            dE_dT = (E_vals[i+1] - E_vals[i]) / deltaT
        elseif i == nT
            # 最後の点：後退差分
            dE_dT = (E_vals[i] - E_vals[i-1]) / deltaT
        else
            # 中央差分
            dE_dT = (E_vals[i+1] - E_vals[i-1]) / (2*deltaT)
        end
        push!(Cv_vals, dE_dT)
    end

    # プロット例
    p1 = plot(T_range, F_vals, xlabel="T", ylabel="Free Energy", title="Free Energy F")
    p2 = plot(T_range, E_vals, xlabel="T", ylabel="Average Energy", title="<E>")
    p3 = plot(T_range, Cv_vals, xlabel="T", ylabel="Specific Heat", title="Cv")
    p4 = plot(T_range, S_vals, xlabel="T",ylabel="Entropy",title="S")
    display(plot(p4))
    
    println("各シェルの情報:")
    for shell in shells
        println(shell)
    end
    println("全状態数 = ", configs)
end

main(116)