#rows: hyperplane normals

# genrate a random central hyperplane arrangement
# rank r and with n planes
function randomCentral(n, r)
    return rand(n, r) - (1/2)*ones(n,r)
end

#note: this does ~not~ modify A
function make_affine(A)
    n, r = size(A)
    A = vcat(A, zeros(r)')
    A = hcat(A, ones(n+1))
    return A
end
#rows: hyperplane normals
function randomAffine(n, r)
    A =  rand(n, r-1) - (1/2)*ones(n,r-1)
    return make_affine(A)
end
