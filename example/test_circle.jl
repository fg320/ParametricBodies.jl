using WaterLily
using StaticArrays
using ParametricBodies # New package
using Plots; gr()
using GLMakie

# Note: measurements (m) are actually in computational cells

# Domain size
n = 3*2^6 # m
m = 2^7 # m

# Physics conditions
Re = 250;
U = 3; # m/s 
alpha = 0; # inflow angle, degrees
U_comp = (U*cos(alpha*π/180), U*sin(alpha*π/180));

# Circle geometry
radius, center = m/8, m/2;
sdf = (x,t)->√sum(abs2, x .- center) - radius
body = AutoBody(sdf);

# Memory
using CUDA; 
@assert CUDA.functional()
memory = Array; # Array (CPU) or CuArray (GPU)

# Simulation
dims_domain = (n,m);
u_BC = U_comp;
L = radius; # m
# U=√sum(abs2,u_BC);

sim = Simulation(dims_domain,
                u_BC, 
                L;
                Δt=0.25, # Real-time
                ν=U*L/Re,
                U=U,
                ϵ=1,
                uλ=(i,x)->u_BC[i],
                body=body,
                T=Float32,
                mem=memory);

# Defaults  # (dims,   
            # u_BC, 
            # L;
            # Δt=0.25,
            # ν=0.,
            # U=√sum(abs2,u_BC),
            # ϵ=1,
            # uλ=(i,x)->u_BC[i],
            # body=NoBody(),
            # T=Float32,
            # mem=Array)

# Flow 
# Fluid fields
sim.flow.u; # velocity vector field
sim.flow.u⁰; # previous velocity
sim.flow.f; # force vector
sim.flow.p; # pressure scalar field
sim.flow.σ; # divergence scalar
# BDIM fields
sim.flow.V; # body velocity vector
sim.flow.σᵥ; # body velocity divergence
sim.flow.μ₀; # zeroth-moment vector
sim.flow.μ₁; # first-moment tensor field
# Non-fields
sim.flow.U; # domain boundary values
sim.flow.Δt; # time step (stored in CPU memory)
sim.flow.ν; # kinematic viscosity

# Matrices
WaterLily.∮nds(sim.flow.p,sim.flow.V,sim.body,0) # Thrust and Side forces (last argument is time)
# others

# Running in time
time_tot = 10 # seconds
time_tot_nd = time_tot*U/L
sim.U;
sim.L;
sim_step!(sim,time_tot_nd; remeasure=true,verbose=true) # Second arg: total non-dimentional time
sim_time(sim) # non-dimentional time
print("Total time (seconds): ")
# print(sim_time(sim)*L/U)
# print(" or ")
print(WaterLily.time(sim))

# Plotting - body at the current time
include("methods.jl")
plot_body!(sim);
# plot_vorticity!(sim);
# plot_pressure!(sim);
# plot_u_x!(sim);
# plot_u_y!(sim);
# plot_u_mag!(sim);
# plot_u_mag_norm!(sim); # Specified clims in methods.jl

# Thrust force and Side force
# begin
# 	function get_force(sim,t)
# 		sim_step!(sim,t*U/L,remeasure=true) # Non-dimensional time
# 		return WaterLily.∮nds(sim.flow.p,sim.flow.V,sim.body,t) # Real time
# 	end
#     time_period = 200 # seconds
#     time_init = WaterLily.time(sim)
# 	forces = [get_force(sim,t) for t ∈ 
#               range(time_init,time_init+time_period,length=time_period)]; # Real time
#     "Got Forces"

#     # Split forces into thrust and side force arrays
#     thrust_forces = [force[1] for force in forces]
#     side_forces = [force[2] for force in forces]
# end

# # Thrust and Side Force Plot
# fig = GLMakie.Figure()
# ax1 = fig[1,1] = GLMakie.Axis(fig, xlabel = "Real Time [sec]", ylabel = "Thrust [N]")
# ax2 = fig[1,1] = GLMakie.Axis(fig, xlabel = "Real Time [sec]", ylabel = "Side Force [N]")
# # ax = GLMakie.Axis(fig[1, 1])
# GLMakie.scatter!(ax1, range(time_init,time_init+time_period,length=time_period), 
#                 thrust_forces, color = :red)
# GLMakie.scatter!(ax2, range(time_init,time_init+time_period,length=time_period), 
#                 side_forces, color = :blue)
# ax2.yaxisposition = :right
# ax2.yticklabelalign = (:left, :center)
# ax2.xticklabelsvisible = false
# ax2.xticklabelsvisible = false
# ax2.xlabelvisible = false
# GLMakie.linkxaxes!(ax1,ax2)
# display(fig)
