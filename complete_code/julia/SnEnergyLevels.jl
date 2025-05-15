using LinearAlgebra
using SpecialFunctions
using QuadGK
using Plots
using Distributed

const ALPHA = 137.035999

# 物理定数の定義
const Z = 50
const Mcc = 939  # MeV
const hbarc = 197.33  # MeV fm

function E_HO(n::Int, l::Int, hbaro::Float64)
    return hbaro * (l + 2n + 1.5)
end

function hyper(n::Int,γ::Float64,x::Float64)
    nn=-n
    if n==0
        return 1.0
    elseif nn==1
        return 1 - x/γ
    else
        Mnm2 = 1.0
        Mnm1 = 1 - x/γ
        Mn=0.0
        for k in 1:(nn-1)
            Mn = ((2k - 2 + γ - x)*Mnm1-(k - 1)*Mnm2)/(k - 1 + γ)
            Mnm2, Mnm1=Mnm1, Mn
        end
        return Mn
    end
end

function hyper1(n::Int, γ::Float64, x::Float64)
    M0 = 1.0
    n == 0 && return M0
    
    M1 = 1.0 - x / γ
    n == 1 && return M1
    
    M_prev2, M_prev1 = M0, M1
    for k in 2:n
        Mn = ((2k + γ - x - 2) * M_prev1 - (k - 1) * M_prev2) / (k + γ - 1)
        M_prev2, M_prev1 = M_prev1, Mn
    end
    
    return M_prev1
end

#基底
function radial(n::Int, l::Int, r::Float64, α::Float64)
    N = sqrt((2 * gamma(l + 1.5 + n)) / (factorial(n) * (gamma(l + 1.5))^2))
    q = α * r
    return sqrt(α) * N * q^(l + 1) * exp(-0.5 * q^2) * hyper1(n, l + 1.5, q^2)
end

@inline function f(r::Float64, R::Float64, a::Float64)
    return 1.0 / (1.0 + exp((r - R) / a))
end

@inline function df(r::Float64, R::Float64, a::Float64)
    X = exp((r - R) / a)
    return -X / (a * (1.0 + X)^2)
end

@inline function U_coul(Z::Int, r::Float64, R::Float64)
    if r < R
        return 0.5 * ((Z - 1) * ((3R^2 - r^2) / R^3)) / ALPHA
    else
        return (Z - 1) / (ALPHA * r)
    end
end

function total_pot(r::Float64, l::Int, j::Float64, u0::Float64, uls::Float64, r0::Float64, R::Float64, a::Float64, Z::Int, tau::Int, moo::Float64)
    HOval = 0.5 * moo * r^2
    fval = f(r, R, a)
    dfval = df(r, R, a)
    ls = 0.5 * (j * (j + 1) - l * (l + 1) - 0.75)
    soval = uls * r0^2 * (1/r) * dfval * ls
    coulval = tau == -1 ? U_coul(Z, r, R) : 0.0
    return u0 * fval + soval + coulval - HOval
end

function matrix_element(n1::Int, l1::Int, j1::Float64, n2::Int, l2::Int, j2::Float64;
    α::Float64, u0::Float64, uls::Float64, r0::Float64, R::Float64, a::Float64, Z::Int, tau::Int, hbaro::Float64, moo::Float64)
    
    (l1 != l2 || j1 != j2) && return 0.0
    
    T_val = (n1 == n2 && l1 == l2) ? E_HO(n1, l1, hbaro) : 0.0
    
    integrand_U(r) = radial(n1, l1, r, α) * total_pot(r, l1, j1, u0, uls, r0, R, a, Z, tau, moo) * radial(n2, l2, r, α)
    rmax = 10.0/α
    U_val, err = quadgk(integrand_U, 1e-6, rmax, rtol=1e-6)
    
    return T_val + U_val
end

function build_H(basis_list::Vector{Tuple{Int,Int,Float64}};
    α::Float64, u0::Float64, uls::Float64, r0::Float64, R::Float64, a::Float64, Z::Int, tau::Int, hbaro::Float64, moo::Float64)
    
    dim = length(basis_list)
    H = zeros(Float64, dim, dim)
    
    Threads.@threads for i in 1:dim
        n1, l1, j1 = basis_list[i]
        for j in i:dim  # 対称行列なので上三角部分のみ計算
            n2, l2, j2 = basis_list[j]
            val = matrix_element(n1, l1, j1, n2, l2, j2;
                α=α, u0=u0, uls=uls, r0=r0, R=R, a=a, Z=Z, tau=tau, hbaro=hbaro, moo=moo)
            H[i,j] = val
            if i != j
                H[j,i] = val  # 対称性を利用
            end
        end
    end
    return H
end

#陽子基底の生成
function C_proton_basis()
    return [
        (0, 0, 0.5),  # 0s_{1/2}
        (0, 1, 1.5),  # 0p_{3/2}
        (0, 1, 0.5),  # 0p_{1/2}
        (0, 2, 2.5),  # 0d_{5/2}
        (1, 0, 0.5),  # 1s_{1/2}
        (0, 2, 1.5),  # 0d_{3/2}
        (0, 3, 3.5),  # 0f_{7/2}
        (0, 3, 2.5),  # 0f_{5/2}
        (1, 1, 1.5),  # 1p_{3/2}
        (1, 1, 0.5),  # 1p_{1/2}
        (0, 4, 4.5)   # 0g_{9/2}
    ]
end

#中性子基底の生成
function C_neutron_basis(N::Int)
    all_neutron_orbitals = [
        (0, 0, 0.5),  # 0s_{1/2}
        (0, 1, 1.5),  # 0p_{3/2}
        (0, 1, 0.5),  # 0p_{1/2}
        (0, 2, 2.5),  # 0d_{5/2}
        (1, 0, 0.5),  # 1s_{1/2}
        (0, 2, 1.5),  # 0d_{3/2}
        (0, 3, 3.5),  # 0f_{7/2}
        (0, 3, 2.5),  # 0f_{5/2}
        (1, 1, 1.5),  # 1p_{3/2}
        (1, 1, 0.5),  # 1p_{1/2}
        (0, 4, 4.5),  # 0g_{9/2}
        (1, 2, 2.5),  # 1d_{5/2}
        (0, 4, 3.5),  # 0g_{7/2}
        (2, 0, 0.5),  # 2s_{1/2}
        (1, 2, 1.5),  # 1d_{3/2}
        (0, 5, 5.5)   # 0h_{11/2}
    ]

    capacities = [
        2,  # 0s_{1/2}
        4,  # 0p_{3/2}
        2,  # 0p_{1/2}
        6,  # 0d_{5/2}
        2,  # 1s_{1/2}
        4,  # 0d_{3/2}
        8,  # 0f_{7/2}
        6,  # 0f_{5/2}
        4,  # 1p_{3/2}
        2,  # 1p_{1/2}
        10, # 0g_{9/2}
        6,  # 1d_{5/2}
        8,  # 0g_{7/2}
        2,  # 2s_{1/2}
        4,  # 1d_{3/2}
        12  # 0h_{11/2}
    ]

    included_orbitals = Tuple{Int, Int, Float64}[]
    included_capacities = []
    remaining_neutrons = N

    for (orbital, capacity) in zip(all_neutron_orbitals, capacities)
        if remaining_neutrons > 0
            push!(included_orbitals, orbital)
            push!(included_capacities, capacity)
            remaining_neutrons -= capacity
        else
            break
        end
    end

    occupation = Dict{Tuple{Int, Int, Float64}, Int}()
    remaining_neutrons = N

    for (orbital, capacity) in zip(included_orbitals, included_capacities)
        if remaining_neutrons > 0
            occupy = min(remaining_neutrons, capacity)
            occupation[orbital] = occupy
            remaining_neutrons -= occupy
        else
            occupation[orbital] = 0
        end
    end

    return all_neutron_orbitals
end

function main()
    # 結果を保存するための配列
    results = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}()
    
    # 並列処理で各質量数の計算を実行
    Threads.@threads for A in 100:2:132
        N = A - Z
        
        # 物理パラメータの計算
        hbaro = 41 / (A^(1/3))  # MeV
        moo = (Mcc * hbaro^2) / (hbarc^2)
        α = sqrt(Mcc * hbaro) / hbarc
        
        # ポテンシャルパラメータ
        u0 = -51 + 33 * (N - Z) / A
        uls = 22 - 14 * (N - Z) / A
        r0 = 1.27
        R = r0 * A^(1/3)
        a = 0.67
        
        # 基底状態の計算
        basis_list_n = C_neutron_basis(N)
        basis_list_p = C_proton_basis()
        
        # ハミルトニアンの構築と対角化
        H_neutron = build_H(basis_list_n; α=α, u0=u0, uls=uls, r0=r0, R=R, a=a, Z=Z, tau=1, hbaro=hbaro, moo=moo)
        H_proton = build_H(basis_list_p; α=α, u0=u0, uls=uls, r0=r0, R=R, a=a, Z=Z, tau=-1, hbaro=hbaro, moo=moo)
        
        vals_n, vecs_n = eigen(Symmetric(H_neutron))
        vals_p, vecs_p = eigen(Symmetric(H_proton))
        
        # 結果の保存
        results[A] = (vals_n, vals_p)
        
        # 結果をファイルに出力
        open("$(A)_Sn_Energy_Levels.txt", "w") do file
            println(file, "Energy Levels for A=$(A) Sn:")
            
            println(file, "\nNeutron Levels (MeV):")
            for (energy, vec) in zip(vals_n, eachcol(vecs_n))
                max_idx = argmax(abs.(vec))
                n, l, j = basis_list_n[max_idx]
                println(file, "Energy: $(round(energy, digits=5)) MeV, Basis: (n=$(n), l=$(l), j=$(j))")
            end
            
            println(file, "\nProton Levels (MeV):")
            for (energy, vec) in zip(vals_p, eachcol(vecs_p))
                max_idx = argmax(abs.(vec))
                n, l, j = basis_list_p[max_idx]
                println(file, "Energy: $(round(energy, digits=5)) MeV, Basis: (n=$(n), l=$(l), j=$(j))")
            end
        end
    end
    
    println("計算が完了しました。")
    return results
end

# プログラムの実行

    main()