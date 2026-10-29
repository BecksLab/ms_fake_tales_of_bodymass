using CSV
using DataFrames
using FoodWebTools
using SpeciesInteractionNetworks

using Random
using Statistics

# --------------------------------------------------
# Settings
# --------------------------------------------------

n_runs = 100
Random.seed!(1234)  # For reproducibility

# --------------------------------------------------
# Load data
# --------------------------------------------------

gb_meta = DataFrame(CSV.File("data/subset/metaweb_GB.csv"))
gb_body = DataFrame(CSV.File("data/subset/bodysizes_GB.csv"))

can_meta = DataFrame(CSV.File("data/subset/metaweb_CAN.csv"))
can_body = DataFrame(CSV.File("data/subset/bodysizes_CAN.csv"))


# --------------------------------------------------
# Build metaweb
# --------------------------------------------------

function build_metaweb(meta, probabilistic::Bool)

    species = unique(vcat(meta.resource, meta.consumer))
    nodes = Unipartite(Symbol.(species))

    if probabilistic
        edges = Probabilistic(
            zeros(Float64, length(species), length(species))
        )
    else
        edges = Binary(
            zeros(Bool, length(species), length(species))
        )
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
# Run ATN and calculate forbidden interactions
# --------------------------------------------------

function run_analysis(meta, body; n_runs=100,
                      probabilistic_metaweb=false)

    # Build the metaweb once
    metaweb = build_metaweb(meta, probabilistic_metaweb)

    results = DataFrame(
        run = Int[],
        n_atn = Int[],
        n_forbidden = Int[],
        forbidden_pct = Float64[]
    )

    for r in 1:n_runs

        # Sample the metaweb if it contains probabilistic links
        sampled_metaweb = probabilistic_metaweb ?
            randomdraws(metaweb) : metaweb

        # Identify basal species (producers)
        gen = values(generality(sampled_metaweb))
        is_producer = collect(gen .== 0)

        # Generate and sample the ATN
        A_atn = lmatrix(
            body.species,
            body.bodymass,
            is_producer;
            probabilistic = true
        )

        atn = SpeciesInteractionNetwork(
            Unipartite(Symbol.(body.species)),
            Probabilistic(A_atn)
        )

        atn = randomdraws(atn)

        # Compare interactions
        atn_interactions = Set(interactions(atn))
        meta_interactions = Set(interactions(sampled_metaweb))

        forbidden = setdiff(
            atn_interactions,
            meta_interactions
        )

        n_atn = length(atn_interactions)
        n_forbidden = length(forbidden)

        forbidden_pct = n_atn > 0 ?
            100 * n_forbidden / n_atn : NaN

        push!(results, (
            r,
            n_atn,
            n_forbidden,
            forbidden_pct
        ))
    end

    return results
end


# --------------------------------------------------
# Run 100 simulations for each location
# --------------------------------------------------

gb_results = run_analysis(
    gb_meta,
    gb_body;
    n_runs = n_runs,
    probabilistic_metaweb = false
)

can_results = run_analysis(
    can_meta,
    can_body;
    n_runs = n_runs,
    probabilistic_metaweb = true
)

# Combine results
gb_results.location .= "GB"
can_results.location .= "Canada"

results = vcat(gb_results, can_results)

# Put location first
select!(results, :location, :)

# Save results
CSV.write("data/forbidden_interactions.csv", results)

# Summary statistics for each location
summary = combine(
    groupby(results, :location),
    :forbidden_pct => mean => :mean_forbidden_pct,
    :forbidden_pct => std => :sd_forbidden_pct,
    :forbidden_pct => median => :median_forbidden_pct
)

println(summary)