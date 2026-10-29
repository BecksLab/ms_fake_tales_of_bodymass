using CSV
using DataFrames
using FoodWebTools
using JLD2
using ProgressMeter
using SpeciesInteractionNetworks
using Random
using Statistics

# --------------------------------------------------
# Configuration
# --------------------------------------------------

n_runs = 10
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
    (
        name = "GB",
        meta = gb_meta,
        body = gb_body,
        probabilistic_metaweb = false,
    ),
    (
        name = "CAN",
        meta = can_meta,
        body = can_body,
        probabilistic_metaweb = true,
    ),
]

# --------------------------------------------------
# Build metawebs
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
            Symbol(meta.resource[j]),
        )

        if probabilistic
            network[interaction...] = meta.score[j]
        else
            network[interaction...] = true
        end
    end

    return network
end

# Build each original metaweb once.
metawebs = Dict(
    loc.name => build_metaweb(loc.meta, loc.probabilistic_metaweb)
    for loc in locations
)

# --------------------------------------------------
# Output DataFrames
# --------------------------------------------------

# Network objects and their metadata are kept together.
function empty_network_df()
    return DataFrame(
        location = String[],
        network_type = String[],
        run = Int[],
        Ropt = Union{Missing, Float64}[],
        gamma = Union{Missing, Float64}[],
        threshold = Union{Missing, Float64}[],
        probabilistic_metaweb = Bool[],
        is_original = Bool[],
        network = Any[],
    )
end

# The six groups saved to JLD2.
CAN_metaweb = empty_network_df()
GB_metaweb = empty_network_df()

CAN_binary_atn = empty_network_df()
GB_binary_atn = empty_network_df()

CAN_probabilistic_atn = empty_network_df()
GB_probabilistic_atn = empty_network_df()

# Store each original metaweb once.
push!(CAN_metaweb, (
    "CAN", "metaweb", 0, missing, missing, missing,
    true, true, metawebs["CAN"],
))

push!(GB_metaweb, (
    "GB", "metaweb", 0, missing, missing, missing,
    false, true, metawebs["GB"],
))

# --------------------------------------------------
# Summary results
# --------------------------------------------------

all_results = DataFrame[]

# --------------------------------------------------
# Run simulations
# --------------------------------------------------

function run_analysis(
    metaweb,
    body;
    location,
    n_runs = 100,
    probabilistic_metaweb = false,
    atn_mode = :probabilistic,
    Ropt = 100.0,
    gamma = 2.0,
    threshold = 0.01,
    network_results,
    metaweb_results,
)
    results = DataFrame(
        location = String[],
        atn_mode = String[],
        run = Int[],
        Ropt = Float64[],
        gamma = Float64[],
        threshold = Union{Missing, Float64}[],
        n_atn = Int[],
        n_forbidden = Int[],
        forbidden_pct = Float64[],
    )

    for r in 1:n_runs

        # Draw a metaweb for this run when the metaweb is probabilistic.
        sampled_metaweb = probabilistic_metaweb ?
            randomdraws(metaweb) : metaweb

        # Save each sampled metaweb alongside its run metadata.
        # The fixed GB metaweb is already saved once above.
        if probabilistic_metaweb
            push!(metaweb_results, (
                location,
                "metaweb",
                r,
                Ropt,
                gamma,
                atn_mode == :binary ? threshold : missing,
                probabilistic_metaweb,
                false,
                sampled_metaweb,
            ))
        end

        # Identify producers from the sampled metaweb.
        gen = generality(sampled_metaweb)

        producer_lookup = Dict(
            Symbol(s) => iszero(g)
            for (s, g) in zip(keys(gen), values(gen))
        )

        is_producer = [
            get(producer_lookup, Symbol(s), false)
            for s in body.species
        ]

        # Generate the ATN.
        if atn_mode == :probabilistic
            A_atn = lmatrix(
                body.species,
                Float64.(body.bodymass),
                is_producer;
                Ropt = Ropt,
                γ = gamma,
                probabilistic = true,
            )

            atn = SpeciesInteractionNetwork(
                Unipartite(Symbol.(body.species)),
                Probabilistic(A_atn),
            )

            # Save the realised ATN, as in the original analysis.
            atn = randomdraws(atn)

        elseif atn_mode == :binary
            A_atn = lmatrix(
                body.species,
                Float64.(body.bodymass),
                is_producer;
                Ropt = Ropt,
                γ = gamma,
                threshold = threshold,
                probabilistic = false,
            )

            atn = SpeciesInteractionNetwork(
                Unipartite(Symbol.(body.species)),
                Binary(Bool.(A_atn .!= 0)),
            )

        else
            error("Unknown ATN mode: $atn_mode")
        end

        # Save the ATN and its metadata in the appropriate DataFrame.
        push!(network_results, (
            location,
            string(atn_mode),
            r,
            Ropt,
            gamma,
            atn_mode == :binary ? threshold : missing,
            probabilistic_metaweb,
            false,
            atn,
        ))

        # Compare ATN interactions against the metaweb used in this run.
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
            forbidden_pct,
        ))
    end

    return results
end

# --------------------------------------------------
# Run the full parameter grid
# --------------------------------------------------

total_combinations = length(locations) * (
    length(Ropt_values) * length(gamma_values) +
    length(Ropt_values) * length(gamma_values) *
    length(threshold_values)
)

@showprogress "Running sensitivity analysis" for loc in locations

    # Select the output DataFrames for this location.
    if loc.name == "CAN"
        metaweb_results = CAN_metaweb
        binary_results = CAN_binary_atn
        probabilistic_results = CAN_probabilistic_atn
    else
        metaweb_results = GB_metaweb
        binary_results = GB_binary_atn
        probabilistic_results = GB_probabilistic_atn
    end

    # Probabilistic ATNs.
    for Ropt in Ropt_values, gamma in gamma_values
        push!(all_results, run_analysis(
            metawebs[loc.name],
            loc.body;
            location = loc.name,
            n_runs = n_runs,
            probabilistic_metaweb = loc.probabilistic_metaweb,
            atn_mode = :probabilistic,
            Ropt = Ropt,
            gamma = gamma,
            network_results = probabilistic_results,
            metaweb_results = metaweb_results,
        ))
    end

    # Binary ATNs.
    for Ropt in Ropt_values,
        gamma in gamma_values,
        threshold in threshold_values

        push!(all_results, run_analysis(
            metawebs[loc.name],
            loc.body;
            location = loc.name,
            n_runs = n_runs,
            probabilistic_metaweb = loc.probabilistic_metaweb,
            atn_mode = :binary,
            Ropt = Ropt,
            gamma = gamma,
            threshold = threshold,
            network_results = binary_results,
            metaweb_results = metaweb_results,
        ))
    end
end

# --------------------------------------------------
# Export summary CSVs
# --------------------------------------------------

results = vcat(all_results...)

CSV.write("data/forbidden_interactions.csv", results)

summary = combine(
    groupby(
        results,
        [:location, :atn_mode, :Ropt, :gamma, :threshold],
    ),
    :forbidden_pct => mean => :mean_forbidden_pct,
    :forbidden_pct => std => :sd_forbidden_pct,
    :forbidden_pct => median => :median_forbidden_pct,
    :forbidden_pct => minimum => :min_forbidden_pct,
    :forbidden_pct => maximum => :max_forbidden_pct,
)

CSV.write("data/forbidden_summary.csv", summary)

# --------------------------------------------------
# Save all six network DataFrames to one JLD2 file
# --------------------------------------------------

JLD2.jldsave(
    "data/networks.jld2";
    CAN_metaweb,
    GB_metaweb,
    CAN_binary_atn,
    GB_binary_atn,
    CAN_probabilistic_atn,
    GB_probabilistic_atn,
)