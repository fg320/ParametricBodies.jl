import GLMakie

# Body plot function
function plot_body!(sim,t=WaterLily.time(sim),resolution=(1400,700))
    #Set up figure
    fig = GLMakie.Figure(;resolution)
    ax = GLMakie.Axis(fig[1, 1]; autolimitaspect=1)
    # GLMakie.hidedecorations!(ax); GLMakie.hidespines!(ax)

    # Plot geometry
    a = sim.flow.σ;
    WaterLily.measure_sdf!(a,sim.body,t)
    f = a[inside(a)]
    colormap = GLMakie.to_colormap([:grey30,(:grey,0.5)])
    display(GLMakie.contourf!(ax,f,levels=[-100,0,1];colormap))
            
    # Get minimum and maximum values of x and y
    x_min, x_max = minimum(axes(f, 1)), maximum(axes(f, 1))
    y_min, y_max = minimum(axes(f, 2)), maximum(axes(f, 2))
    
    GLMakie.xlims!(x_min, x_max) 
    GLMakie.ylims!(y_min, y_max)

    display(fig)
end

# Pressure plot function
function plot_pressure!(sim,resolution=(1400,700))
    #Set up figure
    fig = GLMakie.Figure(;resolution)
    ax = GLMakie.Axis(fig[1, 1]; autolimitaspect=1)
    GLMakie.hidedecorations!(ax); GLMakie.hidespines!(ax)

    # Plot pressure
    a = sim.flow.p;
    f = a[inside(a)]
    
    clims=()
    if length(clims)==2
        @assert clims[1]<clims[2]
        @. f=min(clims[2],max(clims[1],f))
    else
        clims = (minimum(f),maximum(f))
    end
    pltobj = GLMakie.heatmap!(ax,f,colorrange=clims,
                    colormap=:bluesreds,interpolate=true)
            
    # Get minimum and maximum values of x and y
    x_min, x_max = minimum(axes(f, 1)), maximum(axes(f, 1))
    y_min, y_max = minimum(axes(f, 2)), maximum(axes(f, 2))
    
    GLMakie.xlims!(x_min, x_max) 
    GLMakie.ylims!(y_min, y_max)

    # Plot geometry
    a = sim.flow.σ;
    WaterLily.measure_sdf!(a,sim.body,WaterLily.time(sim))
    f = a[inside(a)]
    colormap = GLMakie.to_colormap([:grey30,(:grey,0.5)])
    display(GLMakie.contourf!(ax,f,levels=[-100,0,1];colormap))

    GLMakie.Colorbar(fig[1, 2], pltobj, height=GLMakie.Relative(0.5))

    display(fig)
end

# Velocity x-dir plot function
function plot_u_x!(sim,resolution=(1400,700))
    #Set up figure
    fig = GLMakie.Figure(;resolution)
    ax = GLMakie.Axis(fig[1, 1]; autolimitaspect=1)
    GLMakie.hidedecorations!(ax); GLMakie.hidespines!(ax)

    # Plot velocity x-dir
    a = sim.flow.u[:,:,1];
    f = a[inside(a)]
    
    clims=()
    if length(clims)==2
        @assert clims[1]<clims[2]
        @. f=min(clims[2],max(clims[1],f))
    else
        clims = (minimum(f),maximum(f))
    end
    pltobj = GLMakie.heatmap!(ax,f,colorrange=clims,
                    colormap=:bluesreds,interpolate=true)
            
    # Get minimum and maximum values of x and y
    x_min, x_max = minimum(axes(f, 1)), maximum(axes(f, 1))
    y_min, y_max = minimum(axes(f, 2)), maximum(axes(f, 2))
    
    GLMakie.xlims!(x_min, x_max) 
    GLMakie.ylims!(y_min, y_max)

    # Plot geometry
    a = sim.flow.σ;
    WaterLily.measure_sdf!(a,sim.body,WaterLily.time(sim))
    f = a[inside(a)]
    colormap = GLMakie.to_colormap([:grey30,(:grey,0.5)])
    display(GLMakie.contourf!(ax,f,levels=[-100,0,1];colormap))

    GLMakie.Colorbar(fig[1, 2], pltobj, height=GLMakie.Relative(0.5))

    display(fig)
end

# Velocity y-dir plot function
function plot_u_y!(sim,resolution=(1400,700))
    #Set up figure
    fig = GLMakie.Figure(;resolution)
    ax = GLMakie.Axis(fig[1, 1]; autolimitaspect=1)
    GLMakie.hidedecorations!(ax); GLMakie.hidespines!(ax)

    # Plot velocity y-dir
    a = sim.flow.u[:,:,2];
    f = a[inside(a)]
    
    clims=()
    if length(clims)==2
        @assert clims[1]<clims[2]
        @. f=min(clims[2],max(clims[1],f))
    else
        clims = (minimum(f),maximum(f))
    end
    pltobj = GLMakie.heatmap!(ax,f,colorrange=clims,
                    colormap=:bluesreds,interpolate=true)
            
    # Get minimum and maximum values of x and y
    x_min, x_max = minimum(axes(f, 1)), maximum(axes(f, 1))
    y_min, y_max = minimum(axes(f, 2)), maximum(axes(f, 2))
    
    GLMakie.xlims!(x_min, x_max) 
    GLMakie.ylims!(y_min, y_max)

    # Plot geometry
    a = sim.flow.σ;
    WaterLily.measure_sdf!(a,sim.body,WaterLily.time(sim))
    f = a[inside(a)]
    colormap = GLMakie.to_colormap([:grey30,(:grey,0.5)])
    display(GLMakie.contourf!(ax,f,levels=[-100,0,1];colormap))

    GLMakie.Colorbar(fig[1, 2], pltobj, height=GLMakie.Relative(0.5))

    display(fig)
end

# Velocity magnitude plot function
function plot_u_mag!(sim,resolution=(1400,700))
    #Set up figure
    fig = GLMakie.Figure(;resolution)
    ax = GLMakie.Axis(fig[1, 1]; autolimitaspect=1)
    GLMakie.hidedecorations!(ax); GLMakie.hidespines!(ax)

    # Plot velocity magnitude
    a = .√((sim.flow.u[:,:,1].^2).+(sim.flow.u[:,:,2].^2));
    f = a[inside(a)]
    
    clims=()
    if length(clims)==2
        @assert clims[1]<clims[2]
        @. f=min(clims[2],max(clims[1],f))
    else
        clims = (minimum(f),maximum(f))
    end
    pltobj = GLMakie.heatmap!(ax,f,colorrange=clims,
                    colormap=:bluesreds,interpolate=true)
            
    # Get minimum and maximum values of x and y
    x_min, x_max = minimum(axes(f, 1)), maximum(axes(f, 1))
    y_min, y_max = minimum(axes(f, 2)), maximum(axes(f, 2))
    
    GLMakie.xlims!(x_min, x_max) 
    GLMakie.ylims!(y_min, y_max)

    # Plot geometry
    a = sim.flow.σ;
    WaterLily.measure_sdf!(a,sim.body,WaterLily.time(sim))
    f = a[inside(a)]
    colormap = GLMakie.to_colormap([:grey30,(:grey,0.5)])
    display(GLMakie.contourf!(ax,f,levels=[-100,0,1];colormap))

    GLMakie.Colorbar(fig[1, 2], pltobj, height=GLMakie.Relative(0.5))

    display(fig)
end

# Velocity magnitude normalised plot function (norm by U, inflow magnitude)
function plot_u_mag_norm!(sim,resolution=(1400,700))
    #Set up figure
    fig = GLMakie.Figure(;resolution)
    ax = GLMakie.Axis(fig[1, 1]; autolimitaspect=1)
    GLMakie.hidedecorations!(ax); GLMakie.hidespines!(ax)

    # Plot velocity magnitude normalised
    a = (.√((sim.flow.u[:,:,1].^2).+(sim.flow.u[:,:,2].^2)))./sim.U;
    f = a[inside(a)]
    
    clims=(-1,3)
    if length(clims)==2
        @assert clims[1]<clims[2]
        @. f=min(clims[2],max(clims[1],f))
    else
        clims = (minimum(f),maximum(f))
    end
    pltobj = GLMakie.heatmap!(ax,f,colorrange=clims,
                    colormap=:bluesreds,interpolate=true)
            
    # Get minimum and maximum values of x and y
    x_min, x_max = minimum(axes(f, 1)), maximum(axes(f, 1))
    y_min, y_max = minimum(axes(f, 2)), maximum(axes(f, 2))
    
    GLMakie.xlims!(x_min, x_max) 
    GLMakie.ylims!(y_min, y_max)

    # Plot geometry
    a = sim.flow.σ;
    WaterLily.measure_sdf!(a,sim.body,WaterLily.time(sim))
    f = a[inside(a)]
    colormap = GLMakie.to_colormap([:grey30,(:grey,0.5)])
    display(GLMakie.contourf!(ax,f,levels=[-100,0,1];colormap))

    GLMakie.Colorbar(fig[1, 2], pltobj, height=GLMakie.Relative(0.5))

    display(fig)
end

# Vorticity plot function
function plot_vorticity!(sim,resolution=(1400,700))
    #Set up figure
    fig = GLMakie.Figure(;resolution)
    ax = GLMakie.Axis(fig[1, 1]; autolimitaspect=1)
    GLMakie.hidedecorations!(ax); GLMakie.hidespines!(ax)

    # Plot vorticity
    a = sim.flow.σ;
    @inside a[I] = WaterLily.curl(3,I,sim.flow.u)
    f = a[inside(a)]
    
    clims=()
    if length(clims)==2
        @assert clims[1]<clims[2]
        @. f=min(clims[2],max(clims[1],f))
    else
        clims = (minimum(f),maximum(f))
    end
    pltobj = GLMakie.heatmap!(ax,a[inside(a)],colorrange=clims,
                     colormap=:curl,interpolate=true)
            
    # Get minimum and maximum values of x and y
    x_min, x_max = minimum(axes(f, 1)), maximum(axes(f, 1))
    y_min, y_max = minimum(axes(f, 2)), maximum(axes(f, 2))
    
    GLMakie.xlims!(x_min, x_max) 
    GLMakie.ylims!(y_min, y_max)

    # Plot geometry
    a = sim.flow.σ;
    WaterLily.measure_sdf!(a,sim.body,WaterLily.time(sim))
    f = a[inside(a)]
    colormap = GLMakie.to_colormap([:grey30,(:grey,0.5)])
    display(GLMakie.contourf!(ax,f,levels=[-100,0,1];colormap))

    GLMakie.Colorbar(fig[1, 2], pltobj, height=GLMakie.Relative(0.5))

    display(fig)
end