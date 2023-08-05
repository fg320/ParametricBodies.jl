using WaterLily
function circle(n,m;Re=250,U=1)
    radius, center = m/8, m/2
    body = AutoBody((x,t)->√sum(abs2, x .- center) - radius)
    Simulation((n,m), (U,0), radius; ν=U*radius/Re, body)
end

# include("TwoD_plots.jl")
# sim_gif!(circle(3*2^6,2^7),duration=10,clims=(-5,5),plotbody=true)

sim = circle(3*2^6,2^7);
sim_step!(sim,π)
include("TwoD_plots.jl")
a = sim.flow.σ;
@inside a[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
flood(a[inside(a)],clims=(-5,5))
body_plot!(sim)