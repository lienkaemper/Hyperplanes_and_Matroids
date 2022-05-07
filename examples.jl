h_rand  = rand(10,5)
χ_rand = chirotope(h_rand)
cocircuits_3(h_rand)

h_degen = [1 0; -1 0; 0 1; -1 1; 1 1]
χ_degen = chirotope(h_degen)
σ = [5]
complete_flat(h_degen, σ)
coatoms(h_degen)
topes(h_degen)

C_star = cocircuits_3(h_degen)
C_star_2 = cocircuits_2(h_degen)

n = 5
cocircuit_supports = [(σ, setdiff(collect(1:n), σ)) for σ in coatoms(h_degen)]

basic_cocircuit(h_degen, [1; 4],1 )
