@testset "particle_advection.jl" begin
    using SpeedyWeather
    for backwards in (true, false)
        for dynamics in (true, false)
            # initialize spectral grid, particle advection scheme, model, and simulation
            spectral_grid = SpectralGrid(nlayers=1, nparticles=3)
            particle_advection = ParticleAdvection2D(spectral_grid, backwards=backwards)
            model = BarotropicModel(spectral_grid, dynamics=dynamics; particle_advection)
            simulation = initialize!(model)

            # run particle advection
            run!(simulation)

            @test simulation.prognostic_variables.clock.timestep_counter == simulation.prognostic_variables.clock.n_timesteps
        end
    end
end