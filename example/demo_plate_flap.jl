using WaterLily,StaticArrays
using ParametricBodies # New package
function make_sim(;L=32,Re=1e3,St=0.3,αₘ=-π/18,U=1,n=8,m=4,T=Float32,mem=Array)
    # Map from simulation coordinate x to surface coordinate ξ
    nose,pivot = SA[L,0.5f0m*L],SA[0.25f0L,0]
    θ₀ = T(αₘ+atan(π*St)); h₀=T(L); ω=T(π*St*U/h₀)
    function map(x,t)
        θ = θ₀*cos(ω*t); R = SA[cos(θ) -sin(θ); sin(θ) cos(θ)]
        h = SA[0,h₀*sin(ω*t)]
        ξ = R*(x-nose-h-pivot)+pivot # move to origin and align with x-axis
        return ξ   # reflect to positive y
    end

    function sdf(ξ,t) # Line segment SDF
        p = ξ-SA[clamp(ξ[1],0,L),0] # vector from closest point on [0,L] segment to ξ 
        p'*p-2                      # distance (with thickness offset)
    end
    body = AutoBody(sdf,map)

    Simulation((n*L,m*L),(U,0),L;ν=U*L/Re,body,T,mem)
end
using CUDA; @assert CUDA.functional()
include("viz.jl");Makie.inline!(false);

# sim = make_sim(mem=CuArray);
# fig,viz = body_omega_fig(sim);
# display(fig)
# for _ in 1:200
#     sim_step!(sim,sim_time(sim)+0.1)
#     # update!(viz,sim)
#     fig,viz = body_omega_fig(sim);
#     display(fig)
# end

function make_video(name="out_plate_flap.mp4")
    cycle = range(0,2,100)
    sim = make_sim(mem=CuArray)
    fig,viz = body_omega_fig(sim)
    Makie.record(fig,name,2cycle) do t
        sim_step!(sim,t)
        # update!(viz,sim)
        fig,viz = body_omega_fig(sim);
        display(fig)
    end
end