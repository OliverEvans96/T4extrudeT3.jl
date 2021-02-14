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

t3nodes, t3elems = create_t3_mesh()
nLayers = 3
function extrusionh(x, k)
    # coordinates?
    # extrusion level index
    # return k == 0 ? x 
    # return x .+ [0.3k*sin(x[2]),0.2k*cos(x[1]),k]
    return x .+ [0,0,k]
end
t4nodes, t4elems = T4extrudeT3(t3nodes, t3elems, nLayers, extrusionh, naive=false)

# Write vtks
vtkexportmesh("t3.vtk", t3nodes, t3elems)
vtkexportmesh("t4.vtk", t4nodes, t4elems)
