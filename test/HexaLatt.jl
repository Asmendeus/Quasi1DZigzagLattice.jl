latt = HexaLatt(4, 4)

@show latt

@show getAllSites(latt)

@show getAllNNPairs(latt; boundary=:PBC)
@assert length(getAllNNPairs(latt; boundary=:PBC)) == length(Set(getAllNNPairs(latt; boundary=:PBC)))
@assert length(getAllNNPairs(latt; boundary=:OBC)) == length(Set(getAllNNPairs(latt; boundary=:OBC)))
@assert all(x->x in getAllNNPairs(latt; boundary=:PBC), getAllNNPairs(latt; boundary=:OBC))
@show setdiff(getAllNNPairs(latt; boundary=:PBC), getAllNNPairs(latt; boundary=:OBC))

@show getAllNNNPairs(latt; boundary=:PBC)
@assert length(getAllNNNPairs(latt; boundary=:PBC)) == length(Set(getAllNNNPairs(latt; boundary=:PBC)))
@assert length(getAllNNNPairs(latt; boundary=:OBC)) == length(Set(getAllNNNPairs(latt; boundary=:OBC)))
@assert all(x->x in getAllNNNPairs(latt; boundary=:PBC), getAllNNNPairs(latt; boundary=:OBC))
@show setdiff(getAllNNNPairs(latt; boundary=:PBC), getAllNNNPairs(latt; boundary=:OBC))
