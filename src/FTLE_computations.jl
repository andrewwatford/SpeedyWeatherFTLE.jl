function displacement_gradient_matrix_central!(B, plonds, platds, dist_km; radius=SpeedyWeather.DEFAULT_RADIUS)
    _check_dist_km(dist_km)
    radius = _check_radius(radius)
    length(plonds) == length(platds) ||
        throw(DimensionMismatch("plonds and platds must have the same length"))
    length(plonds) % 4 == 0 ||
        throw(ArgumentError("particle vectors must contain four particles per grid point"))

    npoints = length(plonds) ÷ 4
    size(B) == (2, 2, npoints) ||
        throw(DimensionMismatch("B must have size (2, 2, $npoints)"))

    dfac = radius / (1000 * dist_km) / 2
    @inbounds for i in 1:npoints
        p = 4 * i
        B[1, 1, i] = deg2rad(_wrapped_lon_diff(plonds[p - 3], plonds[p - 2])) *
                     cosd((platds[p - 3] + platds[p - 2]) / 2) * dfac
        B[2, 1, i] = deg2rad(platds[p - 3] - platds[p - 2]) * dfac
        B[1, 2, i] = deg2rad(_wrapped_lon_diff(plonds[p - 1], plonds[p])) *
                     cosd((platds[p - 1] + platds[p]) / 2) * dfac
        B[2, 2, i] = deg2rad(platds[p - 1] - platds[p]) * dfac
    end

    return B
end

function displacement_gradient_matrix_central(plonds, platds, dist_km; radius=SpeedyWeather.DEFAULT_RADIUS)
    B = Array{Float64}(undef, 2, 2, length(plonds) ÷ 4)
    return displacement_gradient_matrix_central!(B, plonds, platds, dist_km; radius)
end

@inline function _largest_cauchy_green_eigenvalue(a, b, c, d)
    c11 = a*a + c*c
    c12 = a*b + c*d
    c22 = b*b + d*d
    return (c11 + c22 + sqrt((c11 - c22)^2 + 4 * c12^2)) / 2
end

function FTLE_over_grid!(ftle, B, time)
    npoints = size(B, 3)
    length(ftle) == npoints || throw(DimensionMismatch("ftle must have length $npoints"))

    duration = _duration_magnitude(time, "time")
    if iszero(duration)
        fill!(ftle, NaN)
        return ftle
    end

    @inbounds for i in 1:npoints
        lmax = _largest_cauchy_green_eigenvalue(B[1, 1, i], B[1, 2, i], B[2, 1, i], B[2, 2, i])
        ftle[i] = log(lmax) / (2 * duration)
    end

    return ftle
end

function FTLE_over_grid(B, time)
    ftle = Vector{Float64}(undef, size(B, 3))
    return FTLE_over_grid!(ftle, B, time)
end

_particle_column(A::Union{NCDatasets.Variable, NCDatasets.CFVariable}, index) = A[:, index]
_particle_column(A::AbstractArray, index) = view(A, :, index)
_particle_column(A, index) = A[:, index]

function _check_particle_column(plonds, platds, time_index)
    length(plonds) == length(platds) ||
        throw(DimensionMismatch("lon and lat columns at time index $time_index have different lengths"))

    for i in eachindex(plonds, platds)
        lon = plonds[i]
        lat = platds[i]
        (ismissing(lon) || ismissing(lat)) &&
            throw(ArgumentError("particle positions at time index $time_index contain missing values"))
        lon isa Real && lat isa Real && isfinite(lon) && isfinite(lat) ||
            throw(ArgumentError("particle positions at time index $time_index must be finite real values"))
    end

    return nothing
end

function _FTLE_from_particles!(
    ftle,
    B,
    plonds_time,
    platds_time,
    time_hours::AbstractVector{<:Real},
    grid_or_npoints,
    dist_km;
    time_indices=Colon(),
    radius=SpeedyWeather.DEFAULT_RADIUS,
)
    _check_dist_km(dist_km)
    radius = _check_radius(radius)
    size(plonds_time) == size(platds_time) ||
        throw(DimensionMismatch("plonds_time and platds_time must have the same size"))

    npoints = _grid_npoints(grid_or_npoints)
    expected_particles = 4 * npoints
    size(plonds_time, 1) == expected_particles ||
        throw(DimensionMismatch("particle arrays must have $expected_particles rows"))
    size(plonds_time, 2) == length(time_hours) ||
        throw(DimensionMismatch("particle arrays and time_hours must have the same number of time samples"))

    _check_time_hours(time_hours)
    selected = _checked_time_indices(time_indices, time_hours)
    size(ftle) == (npoints, length(selected)) ||
        throw(DimensionMismatch("ftle must have size ($npoints, $(length(selected)))"))
    size(B) == (2, 2, npoints) ||
        throw(DimensionMismatch("B must have size (2, 2, $npoints)"))

    for (out_index, time_index) in enumerate(selected)
        plonds = _particle_column(plonds_time, time_index)
        platds = _particle_column(platds_time, time_index)
        _check_particle_column(plonds, platds, time_index)
        displacement_gradient_matrix_central!(B, plonds, platds, dist_km; radius)
        FTLE_over_grid!(view(ftle, :, out_index), B, time_hours[time_index])
    end

    return ftle
end

function _FTLE_from_particles(
    plonds_time,
    platds_time,
    time_hours::AbstractVector{<:Real},
    grid_or_npoints,
    dist_km;
    time_indices=Colon(),
    radius=SpeedyWeather.DEFAULT_RADIUS,
)
    npoints = _grid_npoints(grid_or_npoints)
    _check_time_hours(time_hours)
    selected = _checked_time_indices(time_indices, time_hours)
    ftle = Array{Float64}(undef, npoints, length(selected))
    B = Array{Float64}(undef, 2, 2, npoints)
    _FTLE_from_particles!(ftle, B, plonds_time, platds_time, time_hours, npoints, dist_km; time_indices=selected, radius)
    return ftle, _selected_time_hours(time_hours, selected)
end

function _particle_file_time_hours(ds)
    times = ds["time"][:]
    isempty(times) && throw(ArgumentError("particle file time variable is empty"))

    first_time = first(times)
    if first_time isa Real
        _check_time_hours(times)
        return Float64.(times .- first_time)
    end

    try
        return Float64.((times .- first_time) ./ Hour(1))
    catch
        throw(ArgumentError("particle file time variable must contain real elapsed times or date/time values"))
    end
end

function _validate_particle_file_metadata(ds, dist_km; radius=SpeedyWeather.DEFAULT_RADIUS)
    if haskey(ds.attrib, _DIST_KM_ATTRIBUTE)
        file_dist_km = ds.attrib[_DIST_KM_ATTRIBUTE]
        file_dist_km isa Real ||
            throw(ArgumentError("particle file dist_km metadata must be real-valued"))
        isapprox(Float64(file_dist_km), Float64(dist_km); rtol=0, atol=eps(Float64) * max(abs(dist_km), 1)) ||
            throw(ArgumentError("particle file was generated with dist_km=$file_dist_km, not $dist_km"))
    end
    if haskey(ds.attrib, _RADIUS_ATTRIBUTE)
        file_radius = ds.attrib[_RADIUS_ATTRIBUTE]
        file_radius isa Real ||
            throw(ArgumentError("particle file radius metadata must be real-valued"))
        isapprox(Float64(file_radius), Float64(radius); rtol=0, atol=eps(Float64) * max(abs(radius), 1)) ||
            throw(ArgumentError("particle file was generated with radius=$file_radius, not $radius"))
    end
    if haskey(ds.attrib, _PARTICLE_ORDER_ATTRIBUTE)
        String(ds.attrib[_PARTICLE_ORDER_ATTRIBUTE]) == _PARTICLE_ORDER ||
            throw(ArgumentError("particle file has unsupported FTLE particle order"))
    end
    return nothing
end

function _initial_position_tolerance(expected_plonds, expected_platds)
    max_offset = 0.0
    for i in firstindex(expected_plonds):4:lastindex(expected_plonds)
        max_offset = max(
            max_offset,
            abs(_wrapped_lon_diff(expected_plonds[i], expected_plonds[i + 1])) / 2,
            abs(expected_platds[i + 2] - expected_platds[i + 3]) / 2,
        )
    end
    return max(1.2e-2, 0.05 * max_offset)
end

function _validate_initial_positions(ds, grid_or_spectral_grid, dist_km; radius=SpeedyWeather.DEFAULT_RADIUS)
    grid_or_spectral_grid isa Integer &&
        throw(ArgumentError("initial position validation needs a grid; pass validate_initial_positions=false for npoint-only files"))

    expected_plonds, expected_platds = _initial_particle_positions(grid_or_spectral_grid, dist_km; radius)
    actual_plonds = ds["lon"][:, 1]
    actual_platds = ds["lat"][:, 1]
    tolerance = _initial_position_tolerance(expected_plonds, expected_platds)

    for i in eachindex(expected_plonds)
        lon_error = abs(_wrapped_lon_diff(Float64(actual_plonds[i]), expected_plonds[i]))
        lat_error = abs(Float64(actual_platds[i]) - expected_platds[i])
        max(lon_error, lat_error) <= tolerance ||
            throw(ArgumentError("particle file initial positions do not match the FTLE stencil"))
    end

    return nothing
end

function _validate_particle_file(
    ds,
    grid_or_spectral_grid,
    dist_km;
    validate_initial_positions,
    radius=SpeedyWeather.DEFAULT_RADIUS,
)
    haskey(ds.dim, "particle") || throw(ArgumentError("particle file must contain a particle dimension"))
    haskey(ds.dim, "time") || throw(ArgumentError("particle file must contain a time dimension"))
    for name in ("time", "lon", "lat")
        haskey(ds, name) || throw(ArgumentError("particle file must contain variable $name"))
    end

    npoints = _grid_npoints(grid_or_spectral_grid)
    nparticles = ds.dim["particle"]
    ntimes = ds.dim["time"]
    expected_particles = 4 * npoints
    nparticles == expected_particles ||
        throw(DimensionMismatch("particle file has $nparticles particles, expected $expected_particles"))
    ntimes > 0 || throw(ArgumentError("particle file time dimension is empty"))

    dimnames(ds["time"]) == ("time",) ||
        throw(DimensionMismatch("particle file time variable must have dimension (time)"))
    for name in ("lon", "lat")
        dimnames(ds[name]) == ("particle", "time") ||
            throw(DimensionMismatch("particle file $name variable must have dimensions (particle, time)"))
        Tuple(dimsize(ds[name])) == (nparticles, ntimes) ||
            throw(DimensionMismatch("particle file $name variable has the wrong size"))
    end

    _validate_particle_file_metadata(ds, dist_km; radius)
    validate_initial_positions && _validate_initial_positions(ds, grid_or_spectral_grid, dist_km; radius)
    return nothing
end

function _write_particle_file_metadata(path; dist_km, radius=SpeedyWeather.DEFAULT_RADIUS)
    ds = NCDataset(path, "a")
    try
        ds.attrib[_DIST_KM_ATTRIBUTE] = Float64(dist_km)
        ds.attrib[_RADIUS_ATTRIBUTE] = Float64(radius)
        ds.attrib[_PARTICLE_ORDER_ATTRIBUTE] = _PARTICLE_ORDER
    finally
        close(ds)
    end
    return path
end

function _FTLE_from_particle_file!(
    ftle,
    B,
    path::AbstractString,
    grid_or_spectral_grid,
    dist_km;
    time_indices=Colon(),
    validate_initial_positions=true,
    radius=SpeedyWeather.DEFAULT_RADIUS,
)
    radius = _check_radius(radius)
    ds = NCDataset(path, "r")
    try
        _validate_particle_file(ds, grid_or_spectral_grid, dist_km; validate_initial_positions, radius)
        time_hours = _particle_file_time_hours(ds)
        selected = _checked_time_indices(time_indices, time_hours)
        _FTLE_from_particles!(ftle, B, ds["lon"], ds["lat"], time_hours, grid_or_spectral_grid, dist_km; time_indices=selected, radius)
        return ftle, _selected_time_hours(time_hours, selected)
    finally
        close(ds)
    end
end

function _FTLE_from_particle_file(
    path::AbstractString,
    grid_or_spectral_grid,
    dist_km;
    time_indices=Colon(),
    validate_initial_positions=true,
    radius=SpeedyWeather.DEFAULT_RADIUS,
)
    radius = _check_radius(radius)
    ds = NCDataset(path, "r")
    try
        _validate_particle_file(ds, grid_or_spectral_grid, dist_km; validate_initial_positions, radius)
        time_hours = _particle_file_time_hours(ds)
        selected = _checked_time_indices(time_indices, time_hours)
        npoints = _grid_npoints(grid_or_spectral_grid)
        ftle = Array{Float64}(undef, npoints, length(selected))
        B = Array{Float64}(undef, 2, 2, npoints)
        _FTLE_from_particles!(ftle, B, ds["lon"], ds["lat"], time_hours, npoints, dist_km; time_indices=selected, radius)
        return ftle, _selected_time_hours(time_hours, selected)
    finally
        close(ds)
    end
end
