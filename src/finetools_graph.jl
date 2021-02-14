module finetools_graph
using IterTools

using FinEtools
using JuMP
using Cbc

function isconnected(el1::T, el2::T) where T <: NTuple{3, Int}
    # Two triangles share an edge if they have two nodes in common
    length(el1 ∩ el2) == 2
end

export get_edge_info
function get_edge_info(elems::FESetT3)

    # Index i = element number (index of `elems` and `elem_edges` (outer vector))
    # Index j = element-relative edge number (index of `elem_edges` (inner vectors))
    # Index k = global edge number (index of `edge_elems` and `edge_nodes`)

    nelems = count(elems)
    elem_edges::Vector{Vector{Int}} = fill([], nelems)
    edge_elems::Vector{NTuple{2, NTuple{2, Int}}} = []
    edge_nodes::Vector{NTuple{2, Int}} = []

    k = 0
    for (i1, i2) in subsets(1:nelems, 2)
        el1 = elems.conn[i1]
        el2 = elems.conn[i2]
        if isconnected(el1, el2)
            k += 1

            # Elems share the edge
            j1 = push!(elem_edges[i1], k) |> length
            j2 = push!(elem_edges[i2], k) |> length

            # Edge connects the elems
            push!(edge_elems, (
                (i1, j1),
                (i2, j2),
            ))

            # Store the nodes connected by the edge
            n1_ind, n2_ind = elems.conn[i1] ∩ elems.conn[i2]
            push!(edge_nodes, (n1_ind, n2_ind))
        end
    end
    return edge_elems, edge_nodes, elem_edges
end


function orient_edges(nelems, edge_elems, optimizer=Cbc.Optimizer)
    model = Model(optimizer)
    set_silent(model)

    # TODO: Initial guess heuristic?

    @variable(model, edge_orients[i=1:nelems, j=1:3], Bin)

    # Each element must contain edges of both orientations
    @constraint(
        model, 
        cons_both_orients[i=1:nelems], 
        1 <= sum(edge_orients[i, :]) <= 2
    )

    for ((i1, j1), (i2, j2)) in edge_elems
        el1_edge_orient = edge_orients[i1, j1]
        el2_edge_orient = edge_orients[i2, j2]
        # Adjacent elements must have opposite orientations for their shared edge
        @constraint(
            model,
            el1_edge_orient + el2_edge_orient == 1
        )
    end

    optimize!(model)
    return convert(BitArray, round.(Int, value.(edge_orients)))
end


"""
For each element, determine:
1. The "variant" of the tetrahedron to use
    - Variant 0 = [0, 1, 0]
    - Variant 1 = [1, 0, 1]
2. The "orientation", which is defined as the index 
in the element's node list of the unique node ("start node") 
for which both connected edges have the same orientation 
(which is equal to the variant)
"""
function orient_elems(edge_orients, edge_nodes, elem_nodes, elem_edges)
    nelems = length(elem_nodes)
    nedges = length(edge_orients)
    
    elem_variants::Vector{Int} = []

    # `elem_orients[i]` is the index of the start node
    # in `elem_nodes[i]`
    elem_orients::Vector{Int} = []
    sizehint!(elem_variants, nelems)
    sizehint!(elem_orients, nelems)

    for i = 1:nelems
        elem_node_inds = elem_nodes[i]
        elem_edge_inds = elem_edges[i]
        elem_edge_orients = edge_orients[i,:]

        elem_variant = sum(elem_edge_orients) == 1 ? 0 : 1
        push!(elem_variants, elem_variant)

        # We only have indices for edges between elements,
        # so the number of edges present is the number of adjacent elements
        nadj = length(elem_edge_inds)

        # The start node is opposite this edge ("center edge")
        j = findfirst(≠(elem_variant), elem_edge_orients)

        if j <= nadj
            # The center edge is a known edge

            # Get the nodes in the center edge
            k = elem_edge_inds[j]
            edge_node_inds = edge_nodes[k]
            # Get the local index of the node not in the center edge
            elem_orient = findfirst(
                ∌(edge_node_inds),
                elem_node_inds
            )
            push!(elem_orients, elem_orient)
        else
            if nadj == 2
                # Non-corner boundary element, orientation determined

                # We don't have info for the center edge,
                # but that's okay because 
                # we know that the start node is connected
                # to both other edges, which we do know about.
                ed1_node_inds, ed2_node_inds = [
                    edge_nodes[k]
                    for k in elem_edge_inds
                ]
                elem_orient = findfirst(
                    ∈(ed1_node_inds ∩ ed2_node_inds),
                    elem_node_inds
                )
                push!(elem_orients, elem_orient)
            elseif nadj == 1
                # Corner element, orientation underdetermined
                
                # We know that the one edge shared with another element
                # is not the center edge, (since j > nadj)
                # and therefore contains the start node.
                # We can choose either one, 
                # so pick the first node in the known edge.
                k = elem_edge_inds[1]
                start_node_ind, _ = edge_nodes[k]
                elem_orient = findfirst(
                    isequal(start_node_ind),
                    elem_node_inds
                )
                push!(elem_orients, elem_orient)
            else # nadj == 0
                # Element has no neighbors, so it can have any orientation.
                elem_orient = 1
                push!(elem_orients, elem_orient)
            end
        end
    end

    return elem_variants, elem_orients
end

"""
Create three tetrahedral elements to join two triangles

`lower` = node indices of lower triangle
`upper` = node indices of upper triangle
`variant` ∈ [0, 1] = "vertical" mirroring of the tetrahedrons
`orient` ∈ [1, 2, 3] = rotation of the tetrahedrons around the "vertical" axis
"""
function extrude_one(lower::AbstractVector{FInt}, upper::AbstractVector{FInt}, variant=1, orient=1)
    if variant == 0
        # Swap the triangles to switch variants
        tmp = upper
        upper = lower
        lower = tmp
    end

    if orient != 1
        # Rotate the nodes to change the orientation
        lower = circshift(lower, 1-orient)
        upper = circshift(upper, 1-orient)
    end

    return [
        lower[1] lower[2] lower[3] upper[2]
        lower[1] upper[2] lower[3] upper[3]
        lower[1] upper[2] upper[3] upper[1]
    ]

end


function doextrude(fens, fes::FESetT3, nLayers, extrusionh, naive=true, optimizer=nothing)
    nfes = count(fes)
    nn1 = count(fens);
    nnt = nn1*(nLayers+1);
    ngc = 3 * nfes * nLayers;
    hconn = zeros(FInt, ngc, 4);
    conn = connasarray(fes)
    xyz = zeros(FFlt, nnt, 3);

    # Compute new node coordinates
    for k=0:nLayers
        for j=1:nn1
            f=j+k*nn1;
            xyz[f, :] = extrusionh(fens.xyz[j, :], k);
        end
    end

    # Compute the correct variants and orientations for each 
    elem_variants = fill(1, nfes)
    elem_orients = fill(1, nfes)
    if !naive
        println("OH BOY")
        optimizer = something(optimizer, Cbc.Optimizer)
        edge_elems, edge_nodes, elem_edges = get_edge_info(fes)
        edge_orients = orient_edges(nfes, edge_elems)
        elem_variants, elem_orients = orient_elems(edge_orients, edge_nodes, fes.conn, elem_edges)
    end

    # Determine which nodes comprise the new elements
    gc=1;
    for k=1:nLayers
        for i=1:nfes
            lower = conn[i, :] .+ (k-1)*nn1
            upper = conn[i, :] .+ k*nn1
            variant = elem_variants[i]
            orient = elem_orients[i]
            hconn[gc:gc+2, :] = extrude_one(lower, upper, variant, orient)
            gc += 3;
        end
    end

    efes = FESetT4(hconn);
    efens = FENodeSet(xyz);
    return efens, efes
end


export T4extrudeT3
"""
    T4extrudeT3(fens::FENodeSet,  fes::FESetT3, nLayers::FInt, extrusionh::Function)

Extrude a mesh of triangles (T3) into a mesh of tetrahedra (T4).

`nLayers` = the number of newly created layers (in addition to the original, `k=0`)
`extrusionh(x::Vector{Float64}, k::FInt)::Vector{Float64}` computes the coordinates
of the the new node given the coordinates of the original node `x`, and the extruion level `k`.

If `naive = true`, all triangles in the mesh will be extruded into tetrahedrons of the
same orientation, which will likely generate non-matching edges 
(i.e. two elements parially have overlapping boundaries with different node connectivities).
Set `naive = false` to generate an extruded mesh with correctly matched edges.
Edge-matching requires JuMP and Cbc (or another mixed-integer optimizer).

`optimizer` = the optimizer to use if `naive = false`.
Defaults to `Cbc.Optimizer` if Cbc is installed
"""
function T4extrudeT3(fens::FENodeSet, fes::FESetT3, nLayers::FInt, extrusionh::F; naive::Bool=true, optimizer=nothing) where {F<:Function}
    println("Override extrude!")
    id = vec([i for i in 1:count(fens)])
    cn = connectednodes(fes);
    id[cn[:]] = vec([i for i in 1:length(cn)]);
    t3fes = deepcopy(fes);
    updateconn!(t3fes, id);
    t3fens = FENodeSet(fens.xyz[cn[:], :]);

    return doextrude(t3fens, t3fes, nLayers, extrusionh, naive, optimizer);
end

end # module
