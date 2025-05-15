using LinearAlgebra
using SpecialFunctions
using QuadGK
using Plots
const ALPHA=137.035999

function E_HO(n::Int,l::Int,hbaro::Float64)
    return hbaro*(l+2*n+1.5)
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
    # 初期条件
    M0 = 1.0
    if n == 0
        return M0
    end
    
    M1 = 1.0 - x / γ
    if n == 1
        return M1
    end
    
    # 漸化式で計算
    M_prev2 = M0
    M_prev1 = M1
    Mn = 0.0
    for k in 2:n
        Mn = ((2*k + γ - x - 2) * M_prev1 - (k - 1) * M_prev2) / (k + γ - 1)
        M_prev2 = M_prev1
        M_prev1 = Mn
    end
    
    return Mn
end


#基底
function radial(n::Int,l::Int,r::Float64,α::Float64)
    N=sqrt((2*gamma(l+1.5+n))/(factorial(n)*(gamma(l+1.5))^2))
    q = α * r
    return sqrt(α)*N*q^(l+1)*exp(-0.5*q^2)*hyper1(n,l+1.5,q^2)
end

function f(r::Float64,R::Float64,a::Float64)
    return 1.0 / (1.0 + exp((r - R) / a))
end

function df(r::Float64,R::Float64,a::Float64)
    X = exp((r - R) / a)
    return -1.0 * (1.0 / (1.0 + X)^2) * (X / a)
end

function U_coul(Z::Int,r::Float64,R::Float64)
    if(r<R)
        return 0.5*((Z-1)*((3R^2-r^2)/(R^3)))/ALPHA
    else
        return ((Z-1)*r)/(ALPHA)
    end
end

function total_pot(r::Float64,l::Int,j::Float64,u0::Float64,uls::Float64,r0::Float64,R::Float64,a::Float64,Z::Int,tau::Int,moo::Float64)
    HOval=0.5*moo*r^2
    fval=f(r,R,a)
    dfval=df(r,R,a)
    ls=0.5*(j*(j+1)-l*(l+1)-0.75)
    soval=uls*r0^2*(1/r)*dfval*ls
    coulval=0.0
    if tau==-1
        coulval=U_coul(Z,r,R)
    elseif tau==+1
        coulval=0.0
    end
    return u0*fval + soval + coulval - HOval
end

function matrix_element(n1::Int,l1::Int,j1::Float64,n2::Int,l2::Int,j2::Float64;
    α::Float64,u0::Float64,uls::Float64,r0::Float64,R::Float64,a::Float64,Z::Int,tau::Int,hbaro::Float64,moo::Float64)
    if (l1!=l2)||(j1!=j2)
        return 0.0
    end
    integrand_U(r)=radial(n1,l1,r,α)*total_pot(r,l1,j1,u0,uls,r0,R,a,Z,tau,moo)*radial(n2,l2,r,α)
    T_val = 0.0
    U_val = 0.0
    if n1==n2 && l1==l2
        T_val = E_HO(n1,l1,hbaro)
    end
    rmax=10.0/α
    U_val ,err =quadgk(r->integrand_U(r),1e-6,rmax;rtol=1e-6)
    return T_val + U_val
end

function build_H(basis_list::Array{Tuple{Int,Int,Float64},1};
    α::Float64,u0::Float64,uls::Float64,r0::Float64,R::Float64,a::Float64,Z::Int,tau::Int,hbaro::Float64,moo::Float64)

    dim = length(basis_list)
    H = Matrix{Float64}(undef, dim, dim)

    for i in 1:dim
        n1,l1,j1=basis_list[i]
        for j in 1:dim
            n2,l2,j2=basis_list[j]
            H[i,j]=matrix_element(n1,l1,j1,n2,l2,j2;α=α,u0=u0,uls=uls,r0=r0,R=R,a=a,Z=Z,tau=tau,hbaro=hbaro,moo=moo)
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

    return included_orbitals, occupation
end


function main()
    #核子数
    Z=50
    for A in (100:2:132)
        N=A-Z

        #物理量
        hbaro=41/(A^(1/3)) #Mev
        Mcc=939 #MeV
        hbarc=197.33 #MeV fm

        #ポテンシャルパラメーター
        u0=-51+33*(N-Z)/A
        uls=22-14*(N-Z)/A
        r0=1.27
        R=r0*A^(1/3)
        a=0.67
        tau=1 #中性子1,陽子-1
        moo=(Mcc*hbaro^2)/(hbarc^2)
        α=sqrt(Mcc*hbaro)/hbarc


        basis_list_n,occupation_n = C_neutron_basis(N)
        basis_list_p = C_proton_basis()
        tau=1
        H_neutron = build_H(basis_list_n;α=α,u0=u0,uls=uls,r0=r0,R=R,a=a,Z=Z,tau=tau,hbaro=hbaro,moo=moo)
        vals_n,vecs_n = eigen(H_neutron)
        tau=-1
        H_proton = build_H(basis_list_p;α=α,u0=u0,uls=uls,r0=r0,R=R,a=a,Z=Z,tau=tau,hbaro=hbaro,moo=moo)
        vals_p,vecs_p = eigen(H_proton)
        
        #エネルギーの確認
        println("$(A)Eigenvalues (Single-particle energies):")
        println(vals_n)
        


        # --- 以下、エネルギー準位図の描画 ---
        # vals は昇順になっているはずなので確認しつつ使う
        # もし降順であれば sort(vals) などで並び替えてください
        sorted_vals_n = sort(vals_n)
        sorted_vals_p = sort(vals_p)

        # plotの初期化
        # seriestype=:scatter で点を打つ方法もありますが、
        # バーのように描きたい場合は seriestype=:stem や :sticks などでもOK
        # 今回は細い横線だけが見えるようにする例を示します。
        fig = plot(title="$(A)Sn Energy Levels Comparison: Neutrons vs Protons",
                    xlabel="", ylabel="Energy (MeV)",
                    legend=false, size=(800, 600),xlims=(0.0,0.5))
        

        #中性子のエネルギー準位プロット
        for (idx, energy) in enumerate(sorted_vals_n)
            plot!(fig, [0.0, 0.2], [energy, energy], color=:blue, linewidth=2)
        end

        #陽子のエネルギー準位プロット
        for (idx,energy) in enumerate(sorted_vals_p)
            plot!(fig, [0.3, 0.5], [energy, energy], color=:red, linewidth=2)
        end

        # 軸を消したい場合は xaxis=:none の指定なども可能
        plot!(xaxis=false)
        
        savefig("aFixed$(A)_pn1.png")
        end
    return nothing
end

main()