using finetools_graph


nodes, elems = create_t3_mesh()
nodeadj = adjgraph(nodes, elems)
elemadj, edges, node_edges = dualadj(elems)

orientations = finetools_graph.orient_edges(edges, node_edges)