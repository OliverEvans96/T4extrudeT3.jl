using finetools_graph


nodes, elems = create_t3_mesh()
nodeadj = adjgraph(nodes, elems)
elem_nodes = elems.conn
edge_elems, edge_nodes, elem_edges = dualadj(elems)

edge_orients = orient_edges(count(elems), edge_elems)
elem_variants, elem_orients = orient_elems(edge_orients, edge_nodes, elem_nodes, elem_edges)