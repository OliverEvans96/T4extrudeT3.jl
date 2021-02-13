using finetools_graph

nodes, elems = create_t3_mesh()
nodeadj = adjgraph(nodes, elems)
edgeadj = dualadj(elems)