using LinearAlgebra

#input: hyperplane arrangement, rows are hyperplane normals
#output: dictionary, σ -> sign(det(ha_sigma))

#this will work, even if matroid is not uniform
function chirotope(ha)
    (n, r) = size(ha)
    chirotope =Dict()
    for σ in powerset(collect(1:n), r, r) #all subsets of 1:n of size r
        σ_2 = sort(σ)
        chirotope[σ_2] = sign(det(ha[σ_2, :]))
    end
    return chirotope
end

function fromHyperplanes(ha)
    M = OrientedMatroid()
    M.chirotope = chirotope(ha)
    M.n = size(ha, 1)
    M.r = rank(ha)
    return M

end
