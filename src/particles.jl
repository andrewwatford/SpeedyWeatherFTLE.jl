function _check_stencil_centers(londs, latds, dist_km)
    length(londs) == length(latds) ||
        throw(DimensionMismatch("londs and latds must have the same length"))

    dist_km = _check_dist_km(dist_km)
    delta_lat = rad2deg(dist_km * 1000 / Re)
    delta_lat <= _MAX_STENCIL_DEGREES ||
        throw(ArgumentError("dist_km=$dist_km is too large for the local spherical stencil"))

    for i in eachindex(londs, latds)
        lon = londs[i]
        lat = latds[i]
        lon isa Real && lat isa Real && isfinite(lon) && isfinite(lat) ||
            throw(ArgumentError("stencil centers must be finite real values"))
        abs(lat) + delta_lat < 90 ||
            throw(ArgumentError("dist_km=$dist_km makes the stencil cross a pole at latitude $lat"))
        delta_lon = delta_lat / abs(cosd(lat))
        isfinite(delta_lon) && delta_lon <= _MAX_STENCIL_DEGREES ||
            throw(ArgumentError("dist_km=$dist_km creates a nonlocal east/west stencil at latitude $lat"))
    end

    return length(londs), delta_lat
end

function _initial_particle_positions!(plonds, platds, londs, latds, dist_km)
    npoints, delta_lat = _check_stencil_centers(londs, latds, dist_km)
    expected_length = 4 * npoints
    length(plonds) == expected_length || throw(DimensionMismatch("plonds must have length $expected_length"))
    length(platds) == expected_length || throw(DimensionMismatch("platds must have length $expected_length"))

    @inbounds for i in 1:npoints
        p = 4 * i
        delta_lon = delta_lat / cosd(latds[i])
        plonds[p - 3] = londs[i] + delta_lon
        platds[p - 3] = latds[i]
        plonds[p - 2] = londs[i] - delta_lon
        platds[p - 2] = latds[i]
        plonds[p - 1] = londs[i]
        platds[p - 1] = latds[i] + delta_lat
        plonds[p] = londs[i]
        platds[p] = latds[i] - delta_lat
    end

    return plonds, platds
end

function _initial_particle_positions!(plonds, platds, grid_or_spectral_grid, dist_km)
    londs, latds = RingGrids.get_londlatds(_spatial_grid(grid_or_spectral_grid))
    return _initial_particle_positions!(plonds, platds, londs, latds, dist_km)
end

function _initial_particle_positions(londs, latds, dist_km)
    plonds = Vector{Float64}(undef, 4 * length(londs))
    platds = similar(plonds)
    return _initial_particle_positions!(plonds, platds, londs, latds, dist_km)
end

function _initial_particle_positions(grid_or_spectral_grid, dist_km)
    npoints = _grid_npoints(grid_or_spectral_grid)
    plonds = Vector{Float64}(undef, 4 * npoints)
    platds = similar(plonds)
    return _initial_particle_positions!(plonds, platds, grid_or_spectral_grid, dist_km)
end

function _perturb_positions!(particles, londs, latds, dist_km)
    npoints, delta_lat = _check_stencil_centers(londs, latds, dist_km)
    expected_length = 4 * npoints
    length(particles) == expected_length ||
        throw(DimensionMismatch("particles must have length $expected_length"))

    @inbounds for i in 1:npoints
        p = 4 * i
        delta_lon = delta_lat / cosd(latds[i])
        particles[p - 3] = Particle(londs[i] + delta_lon, latds[i])
        particles[p - 2] = Particle(londs[i] - delta_lon, latds[i])
        particles[p - 1] = Particle(londs[i], latds[i] + delta_lat)
        particles[p] = Particle(londs[i], latds[i] - delta_lat)
    end

    return particles
end
