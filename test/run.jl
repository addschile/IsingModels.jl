using Plots
using IsingModels

function mag_callback!(model,m)
  m[] += calc_avg_magnetization(model)
end

h = 0.0
J = 1.0
Tc = 2.269185*J
Tog = 0.1*Tc
n = 100

model = IsingModel2D(h,J,Tog,n,seed=convert(UInt64,42))
randomdata!(model)
anim = @animate for i in 1:100
  mcsweep(model)
  heatmap!(model.data, c = :cool)
end
gif(anim, "low_T_ising.gif")

allup!(model)
mcsweeps(1000,model)

Tfacs = range(0.1,1.;step=0.1)
ms = zeros(length(Tfacs))
for i in 1:length(Tfacs)
  T = Tfacs[i]*Tc
  model.T = T
  mcsweeps(100,model)
  m = Ref{Float64}(ms[i])
  mcsweeps(1000,model;cb=mag_callback!,cbargs=m)
  ms[i] = abs(m[] / 1000.0)
  println(ms[i])
end

plt = plot(Tfacs,ms)
savefig("avg_mag.png")

anim = @animate for i in 1:100
  mcsweep(model)
  heatmap!(model.data, c = :cool)
end
gif(anim, "critical_T_ising.gif")
