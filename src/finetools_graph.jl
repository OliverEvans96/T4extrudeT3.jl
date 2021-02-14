module finetools_graph

using FinEtools
using SymRCM: adjgraph
using IterTools

using JuMP
using Cbc

export create_t3_mesh
function create_t3_mesh()
     # Start from tetrahedral mesh
    xmin, xmax = -6.0, 6.0
    ymin, ymax = -2.5, 2.5
    lf = xmax - xmin
    wf = ymax - ymin
    nxf = 4
    nyf = 5
    t4nodes, t4elems = T4block(lf, wf, 1.0, nxf, nyf, 1)
    t4nodes.xyz[:,1] .+= xmin
    t4nodes.xyz[:,2] .+= ymin
    # Select boundary
    t3elems = meshboundary(t4elems)
    # Exclude downward-ish faces (the floor)
    bottom_inds = selectelem(t4nodes, t3elems, facing=true, direction=[0.0, 0.0, -1.0])
    t3elems = subset(t3elems, bottom_inds)
    # Create new node set with only selected nodes
    connected = findunconnnodes(t4nodes, t3elems);
    t3nodes, new_numbering = compactnodes(t4nodes, connected);
    renumberconn!(t3elems, new_numbering);
    # Modify original mesh a bit
    f((x, y, z)) = [x y 0.4*(sin(2*x) + cos(2*y))]
    #t3nodes.xyz .= reduce(vcat, f.(eachrow(t3nodes.xyz)))

    return t3nodes, t3elems
end

export adjgraph
function adjgraph(nodes::FENodeSet, elems::AbstractFESet2Manifold)
    num_nodes = count(nodes)
    conn = connasarray(elems)
    return adjgraph(conn, num_nodes)
end

function isconnected(el1::T, el2::T) where T <: NTuple{3, FInt}
    # Two triangles share an edge if they have two nodes in common
    length(el1 ∩ el2) == 2
end

export dualadj
function dualadj(elems::AbstractFESet2Manifold)
    """
    Construct the dual adjacency graph for a triangular mesh.
    In other words, determine which elements share an edge.
    """

    # Index i = element number (index of `elems` and `elem_edges` (outer vector))
    # Index j = element-relative edge number (index of `elem_edges` (inner vectors))
    # Index k = global edge number (index of `edge_elems` and `edge_nodes`)

    nelems = count(elems)
    #elemadj::Vector{Vector{Int}} = fill([], nelems)
    elem_edges::Vector{Vector{Int}} = fill([], nelems)
    edge_elems::Vector{NTuple{2, NTuple{2, Int}}} = []
    edge_nodes::Vector{NTuple{2, Int}} = []

    k = 0
    for (i1, i2) in subsets(1:nelems, 2)
        el1 = elems.conn[i1]
        el2 = elems.conn[i2]
        if isconnected(el1, el2)
            k += 1

            # Elems are adjacent
            # el1_edge_ind = push!(elemadj[i1], i2) |> length
            # el2_edge_ind = push!(elemadj[i2], i1) |> length

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

            # # Save the index of this edge in each elem's list
            # push!()
        end
    end
    return edge_elems, edge_nodes, elem_edges
end


export orient_edges
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

export orient_elems
function orient_elems(edge_orients, edge_nodes, elem_nodes, elem_edges)
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

        println()
        @show elem_variant
        @show elem_edge_orients
        @show j
        @show nadj

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
        @show elem_orients[i]
    end

    return elem_variants, elem_orients
end

export greet
greet() = print("Hello World!")

end # module
