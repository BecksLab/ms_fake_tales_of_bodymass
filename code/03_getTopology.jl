using CSV
using DataFrames
using FoodWebTools
using JLD2
using ProgressMeter
using SpeciesInteractionNetworks
using Statistics

# --------------------------------------------------
# Network topology summary for binary networks
# --------------------------------------------------

function _network_summary(
    N::SpeciesInteractionNetwork{<:Partiteness, <:Binary}
)
    gen = collect(values(SpeciesInteractionNetworks.generality(N)))
    vul = collect(values(SpeciesInteractionNetworks.vulnerability(N)))

    L = links(N)
    S = richness(N)
    l_s = L / S

    tls = trophic_level(N.edges.edges)

    motif_list = motifs(Unipartite, 3)

    S1 = length(findmotif(motif_list[1], N))
    S2 = length(findmotif(motif_list[2], N))
    S4 = length(findmotif(motif_list[4], N))
    S5 = length(findmotif(motif_list[5], N))

    return Dict(
        :richness => richness(N),
        :connectance => SpeciesInteractionNetworks.connectance(N),
        :complexity => complexity(N),
        :trophic_level => maximum(values(tls)),
        :generality => std(gen ./ l_s),
        :vulnerability => std(vul ./ l_s),
        :NDTI => (S1 + S2) == 0 ? NaN : (S1 - S2) / (S1 + S2),
        :NDCI => (S4 + S5) == 0 ? NaN : (S4 - S5) / (S4 + S5),
    )
end

# --------------------------------------------------
# Load all six groups
# --------------------------------------------------

saved = JLD2.load("data/networks.jld2")

group_names = [
    "CAN_metaweb",
    "GB_metaweb",
    "CAN_binary_atn",
    "GB_binary_atn",
    "CAN_probabilistic_atn",
    "GB_probabilistic_atn",
]

# --------------------------------------------------
# Process every network in every group
# --------------------------------------------------

all_topology = DataFrame[]

for group_name in group_names
    println("Processing $group_name")

    df = saved[group_name]
    group_results = DataFrame()

    @showprogress "Calculating network topology" for i in 1:nrow(df)
        network = df.network[i]

        # Retain all metadata except the network object.
        row = select(df[i:i, :], Not(:network))

        # Record the source group.
        row[!, :saved_group] = [group_name]

        # Calculate topology for binary networks.
        if network isa SpeciesInteractionNetwork{<:Partiteness, <:Binary}
            metrics = _network_summary(network)

            for (metric, value) in metrics
                row[!, String(metric)] = [value]
            end
        else
            # Retain probabilistic networks, with binary-only metrics missing.
            for metric in [
                :richness,
                :connectance,
                :complexity,
                :trophic_level,
                :generality,
                :vulnerability,
                :NDTI,
                :NDCI,
            ]
                row[!, String(metric)] = [missing]
            end
        end

        # Append the one-row DataFrame correctly.
        append!(group_results, row; cols = :union)
    end

    push!(all_topology, group_results)
end

# --------------------------------------------------
# Combine and export
# --------------------------------------------------

topology = vcat(all_topology...; cols = :union)

CSV.write("data/topology.csv", topology)