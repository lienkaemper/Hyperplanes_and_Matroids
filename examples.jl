using Random

include("sign_vector.jl")
include("oriented_matroids.jl")
include("conversions.jl")
include("hyperplane_to_matroid.jl")
include("check_axioms.jl")
include("generate_arrangement.jl")



Random.seed!(1)
# generate a random hyperplane arrangement
# each row is the normal vector to a hyperplane
h_rand  = rand(10,3)

#= initialize an oriented matroid M
assign n, rank, chirotope, topes, and cocircuits from the hyperplane
arrangement  ha
=#
M = fromHyperplanes(h_rand)

#now, lets look at stuff!

M.topes
M.cocircuits
M.chirotope

#right now, axiom checking only works for uniform matroids
checkCoCircuitAxioms(M)

#lets try with an affine oriented matroid!

#this generates a rank 3 affine oriented matroid, corresponding to a line
#arrangement in the plane
h_aff = randomAffine(4, 3)
M = fromHyperplanes(h_aff)
affineTopes(M)




# a non-uniform example

h_degen = [1 0; -1 0; 0 1; -1 1; 1 1]

M_non_uniform = fromHyperplanes(h_degen)

M_non_uniform.topes
M_non_uniform.cocircuits
M_non_uniform.chirotope

# an non-uniform affine example: three lines meeting at a point, viewed as affine
h_nuaf = [1 0 0; 0 1 0; 1 1 0; 0 0 1]
M_nuaf = fromHyperplanes(h_nuaf)
affineTopes(M_nuaf)
