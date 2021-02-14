using finetools_graph
using FinEtools

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

nodes, elems = create_t3_mesh()
elem_nodes = elems.conn
edge_elems, edge_nodes, elem_edges = get_edge_info(elems)

edge_orients = orient_edges(count(elems), edge_elems)
elem_variants, elem_orients = orient_elems(edge_orients, edge_nodes, elem_nodes, elem_edges)