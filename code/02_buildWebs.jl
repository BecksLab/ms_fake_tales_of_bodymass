using CSV
using DataFrames
using FoodWebTools
using ProgressMeter
using SpeciesInteractionNetworks
using Random
using Statistics

# --------------------------------------------------
# Settings
# --------------------------------------------------

n_runs = 100
Random.seed!(66)

Ropt_values = collect(50.0:20.0:200.0) 
gamma_values = collect(0.5:0.5:5.5)    
threshold_values = collect(range(0.001, 0.1; length=10))

# --------------------------------------------------
# Load data
# --------------------------------------------------

gb_meta = DataFrame(CSV.File("data/subset/metaweb_GB.csv"))
gb_body = DataFrame(CSV.File("data/subset/bodysizes_GB.csv"))

can_meta = DataFrame(CSV.File("data/subset/metaweb_CAN.csv"))
can_body = DataFrame(CSV.File("data/subset/bodysizes_CAN.csv"))

locations = [
    (name = "GB", meta = gb_meta, body = gb_body,
     probabilistic_metaweb = false),
    (name = "Canada", meta = can_meta, body = can_body,
     probabilistic_metaweb = true)
]


# --------------------------------------------------
# Build metaweb
# --------------------------------------------------

function build_metaweb(meta, probabilistic::Bool)

    species = unique(vcat(meta.resource, meta.consumer))
    nodes = Unipartite(Symbol.(species))

    edges = if probabilistic
        Probabilistic(zeros(Float64, length(species), length(species)))
    else
        Binary(zeros(Bool, length(species), length(species)))
    end

    network = SpeciesInteractionNetwork(nodes, edges)

    for j in 1:nrow(meta)
        interaction = (
            Symbol(meta.consumer[j]),
            Symbol(meta.resource[j])
        )

        if probabilistic
            network[interaction...] = meta.score[j]
        else
            network[interaction...] = true
        end
    end

    return network
end


# --------------------------------------------------
# Run simulations for one parameter combination
# --------------------------------------------------

function run_analysis(
    meta,
    body;
    location,
    n_runs = 100,
    probabilistic_metaweb = false,
    atn_mode = :probabilistic,
    Ropt = 100.0,
    gamma = 2.0,
    threshold = 0.01
)

    metaweb = build_metaweb(meta, probabilistic_metaweb)

    results = DataFrame(
        location = String[],
        atn_mode = String[],
        run = Int[],
        Ropt = Float64[],
        gamma = Float64[],
        threshold = Union{Missing, Float64}[],
        n_atn = Int[],
        n_forbidden = Int[],
        forbidden_pct = Float64[]
    )

    for r in 1:n_runs

        # Sample the metaweb if it is probabilistic
        sampled_metaweb = probabilistic_metaweb ?
            randomdraws(metaweb) : metaweb

        # Identify producers, aligned to body.species
        gen = generality(sampled_metaweb)

        producer_lookup = Dict(
            Symbol(s) => iszero(g)
            for (s, g) in zip(keys(gen), values(gen))
        )

        is_producer = [
            get(producer_lookup, Symbol(s), false)
            for s in body.species
        ]

        # Generate the ATN using the chosen parameters
        if atn_mode == :probabilistic

            A_atn = lmatrix(
                body.species,
                Float64.(body.bodymass),
                is_producer;
                Ropt = Ropt,
                γ = gamma,
                probabilistic = true
            )

            atn = SpeciesInteractionNetwork(
                Unipartite(Symbol.(body.species)),
                Probabilistic(A_atn)
            )

            atn = randomdraws(atn)

        elseif atn_mode == :binary

            A_atn = lmatrix(
                body.species,
                Float64.(body.bodymass),
                is_producer;
                Ropt = Ropt,
                γ = gamma,
                threshold = threshold,
                probabilistic = false
            )

            atn = SpeciesInteractionNetwork(
                Unipartite(Symbol.(body.species)),
                Binary(Bool.(A_atn .!= 0))
            )

        else
            error("Unknown ATN mode: $atn_mode")
        end

        # Calculate forbidden interactions
        atn_interactions = Set(interactions(atn))
        meta_interactions = Set(interactions(sampled_metaweb))

        forbidden = setdiff(atn_interactions, meta_interactions)

        n_atn = length(atn_interactions)
        n_forbidden = length(forbidden)

        forbidden_pct = n_atn > 0 ?
            100 * n_forbidden / n_atn : NaN

        push!(results, (
            location,
            string(atn_mode),
            r,
            Ropt,
            gamma,
            atn_mode == :binary ? threshold : missing,
            n_atn,
            n_forbidden,
            forbidden_pct
        ))
    end

    return results
end


# --------------------------------------------------
# Run parameter sensitivity analysis
# --------------------------------------------------

all_results = DataFrame[]

@showprogress "Running sensitivity analysis" for loc in locations
    # A. Probabilistic ATN: vary Ropt and gamma
    for Ropt in Ropt_values, gamma in gamma_values
        push!(all_results, run_analysis(
            loc.meta, loc.body;
            location = loc.name,
            n_runs = n_runs,
            probabilistic_metaweb = loc.probabilistic_metaweb,
            atn_mode = :probabilistic,
            Ropt = Ropt,
            gamma = gamma
        ))
    end

    # B. Binary ATN: additionally vary threshold
    for Ropt in Ropt_values,
        gamma in gamma_values,
        threshold in threshold_values
        push!(all_results, run_analysis(
            loc.meta, loc.body;
            location = loc.name,
            n_runs = n_runs,
            probabilistic_metaweb = loc.probabilistic_metaweb,
            atn_mode = :binary,
            Ropt = Ropt,
            gamma = gamma,
            threshold = threshold
        ))
    end
end
results = vcat(all_results...)

# Save all simulation-level results
CSV.write("data/forbidden_interactions.csv", results)


# --------------------------------------------------
# Summarise sensitivity results
# --------------------------------------------------

summary = combine(
    groupby(
        results,
        [:location, :atn_mode, :Ropt, :gamma, :threshold]
    ),
    :forbidden_pct => mean => :mean_forbidden_pct,
    :forbidden_pct => std => :sd_forbidden_pct,
    :forbidden_pct => median => :median_forbidden_pct,
    :forbidden_pct => minimum => :min_forbidden_pct,
    :forbidden_pct => maximum => :max_forbidden_pct
)

CSV.write("data/forbidden_summary.csv", summary)