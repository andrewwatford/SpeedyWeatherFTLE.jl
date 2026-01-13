using LinearAlgebra

Re = 6.371e6 # Average Earth radius in meters

function displacement_gradient_matrix_central(plonds, platds, dist_km)
    """
    Compute the displacement gradient matrix given particle positions at a fixed time
    Uses a central difference scheme to do this, with four points per grid cell

    Inputs:
        plonds: longitudes of particles
        platds: latitudes of particles
        dist_km: perturbation in position applied before starting simulation

    Outputs:
        B: (2,2,N) array where N is the number of grid points in the simulation.
        B[:,:,k] is the displacement gradient matrix for the k'th grid point. 
    """

    Ngpoints = length(plonds) ÷ 4 # Number of grid points

    dfac = Re / dist_km / 2

    cos_factor = cos.(deg2rad.(platds)) # Cos(latitude)

    # Derivative of x w.r.t. X
    xX = [deg2rad((plonds[4i-3] - plonds[4i-2]))*cos_factor[i]*dfac for i in 1:Ngpoints];
    # Derivative of y w.r.t. X 
    yX = [deg2rad((platds[4i-3] - platds[4i-2]))*dfac for i in 1:Ngpoints];
    # Derivative of x w.r.t. Y
    xY = [deg2rad((plonds[4i-1] - plonds[4i]))*cos_factor[i]*dfac for i in 1:Ngpoints];
    # Derivative of y w.r.t. Y 
    yY = [deg2rad((platds[4i-1] - platds[4i]))*dfac for i in 1:Ngpoints];

    # Displacement gradient matrix
    B = Array{eltype(xX)}(undef, 2, 2, Ngpoints)
    B[1,1,:] .= xX
    B[1,2,:] .= xY
    B[2,1,:] .= yX
    B[2,2,:] .= yY 
    # ^ Not the most elegant, but explicit 
    
    return B
end

function FTLE_over_grid(B, T)
    """
    Compute FTLE over a grid

    Inputs:
        B: displacement gradient matrix
        T: time after particle release which B corresponds to

    Outputs: 
        FTLE_grid: FTLE at each grid point. Units are 
    """

    Ngpoints = size(B, 3) # Number of grid points

    FTLE_grid = Vector{Float64}(undef, Ngpoints) 

    for k in 1:Ngpoints
        # Gradient tensor
        Bk = B[:,:,k];
        # Right Cauchy-Green tensor
        CG = Bk' * Bk;
        # Largest eigenvalue
        lmax = maximum(eigvals(CG))
        # FTLE - in units of days^-1
        FTLE_grid[k] = log(lmax) / (2 * T)
    end

    return FTLE_grid
end

export Re
