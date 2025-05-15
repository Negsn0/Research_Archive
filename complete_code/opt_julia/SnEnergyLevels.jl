# SnEnergyLevels.jl（最適化版）

using LinearAlgebra
using SpecialFunctions
using QuadGK
using Plots
using Base.Threads

const ALPHA = 137.035999

# —————————————————————————————
# ポテンシャル／軌道パラメータをまとめる構造体
struct PotParams
    α::Float64
    u0::Float64
    uls::Float64
    r0::Float64
    R::Float64
    a::Float64
    Z::Int
    tau::Int
    moo::Float64
    hbaro::Float64
end

# —————————————————————————————
# 1F1 型漸化式：hyper1（再帰版）
function hyper1(n::Int, γ::Float64, x::Float64)
    M0 = 1.0
    n == 0 && return M0
    M1 = 1.0 - x/γ
    n == 1 && return M1
    M_prev2, M_prev1 = M0, M1
    Mn = 0.0
    for k in 2:n
        Mn = ((2*k + γ - x - 2)*M_prev1 - (k - 1)*M_prev2) / (k + γ - 1)
        M_prev2, M_prev1 = M_prev1, Mn
    end
    return Mn
end

# —————————————————————————————
# norm_coeff を事前計算
function norm_coeffs(basis)
    [ sqrt((2*gamma(l+1.5+n)) / (factorial(n)*(gamma(l+1.5))^2))
      for (n,l,j) in basis ]
end

# —————————————————————————————
# radial 関数（inline 化＋norm_coeff キャッシュ利用）
@inline function radial(n::Int, l::Int, r::Float64, α::Float64, norm::Float64)
    q = α * r
    return sqrt(α) * norm * q^(l+1) * exp(-0.5*q^2) * hyper1(n, l+1.5, q^2)
end

# —————————————————————————————
# 全ポテンシャル（Woods–Saxon + spin–orbit + Coulomb + HO）
@inline function total_pot(r::Float64, l::Int, j::Float64, pp::PotParams)
    # Harmonic term
    HO = 0.5 * pp.moo * r^2
    # Woods–Saxon
    f    = 1.0 / (1 + exp((r - pp.R) / pp.a))
    df   = -exp((r - pp.R) / pp.a) / (pp.a * (1 + exp((r - pp.R) / pp.a))^2)
    # spin-orbit
    ls = 0.5*(j*(j+1) - l*(l+1) - 0.75)
    so = pp.uls * pp.r0^2 * (1/r) * df * ls
    # Coulomb (protonだけ)
    coul = pp.tau == -1 ? 0.5*((pp.Z-1)*((3*pp.R^2 - r^2)/(pp.R^3))) / ALPHA : 0.0
    return pp.u0 * f + so + coul - HO
end

# —————————————————————————————
# ハーモニックオシレータ固有エネルギー
@inline function E_HO(n::Int, l::Int, hbaro::Float64)
    return hbaro * (l + 2n + 1.5)
end

# —————————————————————————————
# 行列要素計算（上三角のみ・symmetry 利用・inline／キャッシュ）
function matrix_element(i::Int, j::Int,
                        basis::Vector{<:Tuple{Int,Int,Float64}},
                        norms::Vector{Float64},
                        pp::PotParams)
    n1,l1,j1 = basis[i]
    n2,l2,j2 = basis[j]
    # ℓ, j が異なれば 0
    if l1 != l2 || j1 != j2
        return 0.0
    end
    # 運動エネルギー（対角要素のみ）
    T = (i == j) ? E_HO(n1, l1, pp.hbaro) : 0.0
    # 積分項
    integrand(r) = radial(n1, l1, r, pp.α, norms[i]) *
                   total_pot(r, l1, j1, pp) *
                   radial(n2, l2, r, pp.α, norms[j])
    U, _ = quadgk(integrand, 1e-6, 10.0/pp.α; rtol=1e-6)
    return T + U
end

# —————————————————————————————
# 行列全体をスレッド並列で構築
function build_H(basis::Vector{<:Tuple{Int,Int,Float64}}, pp::PotParams)
    dim = length(basis)
    H = Matrix{Float64}(undef, dim, dim)
    norms = norm_coeffs(basis)
    @threads for i in 1:dim
        for j in i:dim
            val = matrix_element(i, j, basis, norms, pp)
            H[i,j] = val
            H[j,i] = val
        end
    end
    return H
end

# —————————————————————————————
# 基底生成
function C_proton_basis()
    return [
        (0,0,0.5),(0,1,1.5),(0,1,0.5),(0,2,2.5),
        (1,0,0.5),(0,2,1.5),(0,3,3.5),(0,3,2.5),
        (1,1,1.5),(1,1,0.5),(0,4,4.5)
    ]
end

function C_neutron_basis(N::Int)
    all = [
        (0,0,0.5),(0,1,1.5),(0,1,0.5),(0,2,2.5),
        (1,0,0.5),(0,2,1.5),(0,3,3.5),(0,3,2.5),
        (1,1,1.5),(1,1,0.5),(0,4,4.5),(1,2,2.5),
        (0,4,3.5),(2,0,0.5),(1,2,1.5),(0,5,5.5)
    ]
    # 容量リスト
    caps = [2,4,2,6,2,4,8,6,4,2,10,6,8,2,4,12]
    sel = Tuple{Int,Int,Float64}[]
    rem = N
    for (orb, cap) in zip(all, caps)
        if rem > 0
            push!(sel, orb)
            rem -= min(rem, cap)
        end
    end
    return sel
end

# —————————————————————————————
# メインルーチン
function main()
    Z = 50
    for A in 100:2:132
        N = A - Z
        # 物理・ポテンシャルパラメータ
        hbaro = 41 / A^(1/3)
        Mcc   = 939
        hbarc = 197.33
        u0    = -51 + 33*(N-Z)/A
        uls   = 22  - 14*(N-Z)/A
        r0    = 1.27
        R     = r0 * A^(1/3)
        a     = 0.67
        tau   = 1
        moo   = (Mcc * hbaro^2) / hbarc^2
        α     = sqrt(Mcc*hbaro) / hbarc

        # Params 構築
        pp = PotParams(α, u0, uls, r0, R, a, Z, tau, moo, hbaro)

        # Neutron
        basis_n = C_neutron_basis(N)
        Hn = build_H(basis_n, pp)
        vals_n, vecs_n = eigen(Hn)

        # Proton (tau = -1)
        pp = PotParams(α, u0, uls, r0, R, a, Z, -1, moo, hbaro)
        basis_p = C_proton_basis()
        Hp = build_H(basis_p, pp)
        vals_p, vecs_p = eigen(Hp)

        # ファイル出力
        open("0$(A)_Sn_Energy_Levels.txt","w") do io
            println(io, "A = $A")
            println(io, "--- Neutron Levels (MeV) ---")
            for (e, v) in zip(vals_n, eachcol(vecs_n))
                idx = argmax(abs.(v))
                println(io, "E = $(round(e, digits=5)), basis = $(basis_n[idx])")
            end
            println(io, "\n--- Proton Levels (MeV) ---")
            for (e, v) in zip(vals_p, eachcol(vecs_p))
                idx = argmax(abs.(v))
                println(io, "E = $(round(e, digits=5)), basis = $(basis_p[idx])")
            end
        end
    end
    println("Calculation complete.")
    return nothing
end

main()
