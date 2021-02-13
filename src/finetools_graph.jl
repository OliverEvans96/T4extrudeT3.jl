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
    nelems = count(elems)
    conn::Vector{Vector{Int}} = fill([], nelems)
    node_edges::Vector{Vector{Int}} = fill([], nelems)
    edges::Vector{NTuple{2, Int}} = []

    nedges = 0
    for (i1, i2) in subsets(1:nelems, 2)
        e1 = elems.conn[i1]
        e2 = elems.conn[i2]
        if isconnected(e1, e2)
            nedges += 1

            # Nodes are connected to each other
            push!(conn[i1], i2)
            push!(conn[i2], i1)

            # Edge contains the nodes
            push!(edges, (i1, i2))

            # Nodes are on the edge
            push!(node_edges[i1], nedges)
            push!(node_edges[i2], nedges)
        end
    end
    return conn, edges, node_edges
end


function orient_edges(edges, node_edges; optimizer=Cbc.Optimizer)
    nedges = length(edges)
    nnodes = length(node_edges)
    
    model = Model(optimizer)
    set_silent(model)
    # # Which node should be used as the first in the pattern
    # @variable(model, 1 <= orientation <= 3, Int)
    # # Which variant of the pattern to use
    # @variable(model, 1 <= variant <= 3, Bin)
    @variable(model, orients[i=1:nedges], Bin)
    @constraint(model, nodesum[i=1:nnodes], 1 <= sum(orients[node_edges[i]]) <= 2)
    @objective(model, Min, sum(orients))
    optimize!(model)

    edge_orients = convert(BitArray, round.(Int, value.(labels)))
    node_orients = [edge_orients[edges] for edges in node_edges]
end

export greet
greet() = print("Hello World!")

end # module
