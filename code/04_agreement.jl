
using CSV
using DataFrames
using JLD2
using SpeciesInteractionNetworks
using Statistics


# --------------------------------------------------
# Load saved network objects
# --------------------------------------------------

saved = JLD2.load("data/networks.jld2")

required_keys = [
    "CAN_metaweb",
    "GB_metaweb",
    "CAN_binary_atn",
    "GB_binary_atn",
    "CAN_probabilistic_atn",
    "GB_probabilistic_atn",
]

for key in required_keys
    haskey(saved, key) ||
        error("Missing key '$key' in $input_file")
end

# --------------------------------------------------
# Helpers
# --------------------------------------------------

# Convert interactions to a set of directed consumer-resource
# pairs. Species identifiers are retained as stored in the network.
interaction_set(network) = Set(interactions(network))

function agreement_metrics(atn, metaweb)
    atn_links = interaction_set(atn)
    meta_links = interaction_set(metaweb)

    shared = intersect(atn_links, meta_links)
    forbidden = setdiff(atn_links, meta_links)
    missed = setdiff(meta_links, atn_links)
    union_links = union(atn_links, meta_links)

    n_atn = length(atn_links)
    n_metaweb = length(meta_links)
    n_shared = length(shared)
    n_forbidden = length(forbidden)
    n_missed = length(missed)
    n_union = length(union_links)

    return (
        n_atn = n_atn,
        n_metaweb = n_metaweb,
        n_shared = n_shared,
        n_forbidden = n_forbidden,
        n_missed = n_missed,

        forbidden_pct = n_atn > 0 ?
            100 * n_forbidden / n_atn : NaN,

        precision = n_atn > 0 ?
            n_shared / n_atn : NaN,

        recall = n_metaweb > 0 ?
            n_shared / n_metaweb : NaN,

        jaccard = n_union > 0 ?
            n_shared / n_union : NaN,
    )
end

# Match a saved sampled Canadian metaweb to the ATN that used it.
# The original workflow saves the run and parameter values alongside
# each sampled metaweb.
function matching_metaweb(atn_row, metaweb_df)
    candidates = filter(metaweb_df) do row
        row.location == atn_row.location &&
        row.run == atn_row.run &&
        row.is_original == false &&
        row.Ropt == atn_row.Ropt &&
        row.gamma == atn_row.gamma &&
        isequal(row.threshold, atn_row.threshold)
    end

    nrow(candidates) == 1 ||
        error(
            "Expected exactly one sampled metaweb for " *
            "$(atn_row.location), run=$(atn_row.run), " *
            "Ropt=$(atn_row.Ropt), gamma=$(atn_row.gamma), " *
            "threshold=$(atn_row.threshold); found $(nrow(candidates))."
        )

    return candidates.network[1]
end

# --------------------------------------------------
# Collect saved network groups
# --------------------------------------------------

metawebs = vcat(
    saved["GB_metaweb"],
    saved["CAN_metaweb"],
    cols = :union,
)

atns = vcat(
    saved["GB_binary_atn"],
    saved["CAN_binary_atn"],
    saved["GB_probabilistic_atn"],
    saved["CAN_probabilistic_atn"],
    cols = :union,
)

# --------------------------------------------------
# Calculate agreement for every saved ATN
# --------------------------------------------------

results = DataFrame(
    location = String[],
    atn_mode = String[],
    run = Int[],
    Ropt = Float64[],
    gamma = Float64[],
    threshold = Union{Missing, Float64}[],
    n_atn = Int[],
    n_metaweb = Int[],
    n_shared = Int[],
    n_forbidden = Int[],
    n_missed = Int[],
    forbidden_pct = Float64[],
    precision = Float64[],
    recall = Float64[],
    jaccard = Float64[],
)

for row in eachrow(atns)
    # GB uses its fixed original binary metaweb.
    # Canada uses the sampled metaweb saved for this run.
    reference = if row.location == "GB"
        original = filter(
            r -> r.location == "GB" && r.is_original,
            metawebs,
        )
        nrow(original) == 1 ||
            error("Expected one original GB metaweb.")
        original.network[1]

    elseif row.location == "CAN"
        matching_metaweb(row, metawebs)

    else
        error("Unknown location: $(row.location)")
    end

    metrics = agreement_metrics(row.network, reference)

    push!(results, (
        row.location,
        row.network_type,
        row.run,
        row.Ropt,
        row.gamma,
        row.threshold,
        metrics.n_atn,
        metrics.n_metaweb,
        metrics.n_shared,
        metrics.n_forbidden,
        metrics.n_missed,
        metrics.forbidden_pct,
        metrics.precision,
        metrics.recall,
        metrics.jaccard,
    ))
end

# --------------------------------------------------
# Export
# --------------------------------------------------

CSV.write("data/interaction_agreement.csv", results)
