module finetools_graph

using FinEtools
using SymRCM: adjgraph
using IterTools

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
    length(el1 ∩ el2) > 1
end

export dualadj
function dualadj(elems::AbstractFESet2Manifold)
    """
    Construct the dual adjacency matrix for a triangular mesh.
    In other words, determine which elements share an edge
    """
    num_elems = count(elems)
    A::Vector{Vector{FInt}} = fill([], num_elems)

    for (i1, i2) in subsets(1:num_elems, 2)
        e1 = elems.conn[i1]
        e2 = elems.conn[i2]
        if isconnected(e1, e2)
            push!(A[i1], i2)
            push!(A[i2], i1)
        end
    end
    return A
end


export greet
greet() = print("Hello World!")

end # module
