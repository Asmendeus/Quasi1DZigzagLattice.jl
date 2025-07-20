latt = DiagKagoLatt(4, 4)

@show latt

@show getAllSites(latt)

@show getAllNNPairs(latt; boundary=:OBC)
@show getAllNNPairs(latt; boundary=:PBC)

@show getAllNNNPairs(latt; boundary=:OBC)
@show getAllNNNPairs(latt; boundary=:PBC)
