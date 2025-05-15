using PyPlot
gr()
#核子数
Z=50
A=100
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

#woods
function f(r)
    return u0 / (1.0 + exp((r - R) / a))
end

#調和振動子
function g(r)
    return 0.5*moo*(r^2-R^2) + f(R)
end

#井戸
function h(r)
    if r<R
        return f(0)
    else
        return 0
    end
end

plot(xlims=(0.0,10.0),title="Comparison of H.O., F.S.Well, and Woods-Saxon",
    xlabel="r (fm)",ylabel="V (MeV)")
x=0.0:0.001:10

plot!(x,f,label="Woods-Saxon")
plot!(x,g,label="H.O.")
plot!(x,h,label="F.S.well")

PyPlot.savefig("Comp_pt.eps",format="eps")