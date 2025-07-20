latt = SquaLatt(4, 4)

@show latt

@show getAllSites(latt)

@show getAllNNPairs(latt; boundary=:PBC)
@assert all(x->x in getAllNNPairs(latt; boundary=:PBC), getAllNNPairs(latt; boundary=:OBC))
@show setdiff(getAllNNPairs(latt; boundary=:PBC), getAllNNPairs(latt; boundary=:OBC))

@show getAllNNNPairs(latt; boundary=:OBC)
@assert all(x->x in getAllNNNPairs(latt; boundary=:PBC), getAllNNNPairs(latt; boundary=:OBC))
@show setdiff(getAllNNNPairs(latt; boundary=:PBC), getAllNNNPairs(latt; boundary=:OBC))
