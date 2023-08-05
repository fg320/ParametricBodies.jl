using WaterLily
using StaticArrays
using ParametricBodies # New package
using Plots; gr()
using GLMakie

# Note: measurements (m) are actually in computational cells

# Memory and settings
using CUDA; 
@assert CUDA.functional()
mem = Array; # Array (CPU) or CuArray (GPU)
T = Float32

# Domain size (normalised by L)
n = 8
m = 4

# Physics conditions
Re = 10^4;
U = 2; # m/s 
alpha = 0; # inflow angle, degrees
U_comp = (U*cos(alpha*π/180), U*sin(alpha*π/180));

# # Circle geometry
# radius, center = m/8, m/2;
# sdf = (x,t)->√sum(abs2, x .- center) - radius
# body = AutoBody(sdf);

# Geometry
L = 32 # cell units
# Map from simulation coordinate x to surface coordinate ξ
nose = SA[L,0.5f0m*L]
θ = T(π/18) # radians - geometry angle
# θ = T(0) # radians - geometry angle
function map(x,t)
    R = SA[cos(θ) -sin(θ); sin(θ) cos(θ)]
    ξ = R*(x-nose) # move to origin and align with x-axis
    return SA[ξ[1],abs(ξ[2])]    # reflect to positive y
end

# Define foil using NACA0012 profile equation: https://tinyurl.com/NACA00xx
NACA(s) = 0.6f0*(0.2969f0s-0.126f0s^2-0.3516f0s^4+0.2843f0s^6-0.1036f0s^8)
foil(s,t) = L*SA[(1-s)^2,NACA(1-s)]
body = ParametricBody(foil,(0,1);map,T,mem)

# Simulation
dims_domain = (n*L,m*L);
u_BC = U_comp;
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
                mem=mem);

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
# WaterLily.∮nds(sim.flow.p,sim.flow.V,sim.body,0) # Thrust and Side forces (last argument is time)
WaterLily.∮nds_param(sim.flow.p,sim.flow.V,sim.body,ParametricBodies.sdf,0) # For parametric bodies
# others

# Running in time
time_tot = 2 # seconds 
time_tot_nd = time_tot*U/L # Flow passes X lengths of the body L
sim.U;
sim.L;
sim_step!(sim,time_tot_nd; remeasure=true,verbose=true) # Second arg: total non-dimentional time
print("Non-dimentional time")
print(sim_time(sim)) # non-dimentional time
print("Real time")
print(sim_time(sim)*L/U)
# print(" or ")
# print(WaterLily.time(sim))

# Plotting - body at the current time
include("methods.jl")
# plot_body!(sim);
# plot_vorticity!(sim);
# plot_pressure!(sim);
# plot_u_x!(sim);
# plot_u_y!(sim);
# plot_u_mag!(sim);
# plot_u_mag_norm!(sim); # Specified clims in methods.jl

# Thrust force and Side force
begin
	function get_force(sim,t)
		sim_step!(sim,t*U/L,remeasure=true,verbose=true) # Non-dimensional time
		return WaterLily.∮nds_param(sim.flow.p,sim.flow.V,sim.body,ParametricBodies.sdf,t) # Real time
	end
    time_period = 500 # seconds
    time_period_nd = time_period*U/L # Flow passes X lengths of the body L
    time_init = WaterLily.time(sim)
	forces = [get_force(sim,t) for t ∈ 
              range(time_init,time_init+time_period,length=time_period)]; # Real time
    "Got Forces"

    # Split forces into thrust and side force arrays
    thrust_forces = [force[1] for force in forces]
    side_forces = [force[2] for force in forces]
end

# Thrust and Side Force Plot
fig = GLMakie.Figure()
ax1 = fig[1,1] = GLMakie.Axis(fig, xlabel = "Real Time [sec]", ylabel = "Thrust [N]")
ax2 = fig[1,1] = GLMakie.Axis(fig, xlabel = "Real Time [sec]", ylabel = "Side Force [N]")
# ax = GLMakie.Axis(fig[1, 1])
GLMakie.scatter!(ax1, range(time_init,time_init+time_period,length=time_period), 
                thrust_forces, color = :red)
GLMakie.scatter!(ax2, range(time_init,time_init+time_period,length=time_period), 
                side_forces, color = :blue)
ax2.yaxisposition = :right
ax2.yticklabelalign = (:left, :center)
ax2.xticklabelsvisible = false
ax2.xticklabelsvisible = false
ax2.xlabelvisible = false
GLMakie.linkxaxes!(ax1,ax2)
display(fig)

# Assuming aligned inflow
drag_coeff = WaterLily.∮nds_param(sim.flow.p,sim.flow.V,sim.body,ParametricBodies.sdf,WaterLily.time(sim))[1]/(0.5*sim.L*sim.U^2)
lift_coeff = WaterLily.∮nds_param(sim.flow.p,sim.flow.V,sim.body,ParametricBodies.sdf,WaterLily.time(sim))[2]/(0.5*sim.L*sim.U^2)
