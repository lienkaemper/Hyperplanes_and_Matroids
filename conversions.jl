include("oriented_matroids.jl")


function comp_sign(σ, e)
    r = length(σ)
    current = ((r + 1) %2)*2 -1
    original = (minimum([findall(σ .> e); r+1]) % 2)*2 -1
    return current * original
end

#returns the unique signed cocircuit which is disjoint from χ∖e, positve on e
function basic_cocircuit(χ, σ, e,n)
    signature =  SignVector(zeros(n))
    δ = setdiff(σ , [e])
    a = [δ ;[e]]
    signature[e] = χ[σ] * comp_sign(δ,e)
    for f in setdiff(collect(1:n) ,σ)
        b = [δ ;[f]]
        signature[f] = χ[sort(b)] * comp_sign(δ, f)
    end
    return signature
end

#given chirotope, return cocircuits as list of sign vectors
function cocircuits(χ, n)
    cocircs = Vector{SignVector}()
    for (σ, sgn) in χ
        if sgn != 0
            for e in σ
                signature = basic_cocircuit(χ, σ, e, n)
                if signature ∉ cocircs
                    append!(cocircs, [signature, -signature])
                end
            end
        end
    end
    return cocircs
end

#return topes of rank r matroid on ground set n with chirotope χ.
function topes(χ, n, r)
    topes = []
    for (σ, sgn) ∈ χ
        if sgn != 0
            for α ∈ allPM(r)
                T_σα = composition([α[i] * basic_cocircuit(χ, σ, σ[i], n) for i = 1:r]...)
                if T_σα ∉ topes
                    push!(topes, T_σα)
                end
            end
        end
    end
    return topes
end

# M an oriented matroid with circuits filled in
# fills in topes of M, also returns them
function circuitsToTopes!(M::OrientedMatroid)
    # only need to take one circuit from each opposite pair,
    # since orthogonality X⟂Y iff and only if X⟂-Y.
    circuits = [circ[1] for circ in values(M.circuits)]
    n = length(circuits[1])
    topes = SignVector[]
    # something is a tope if it's orthogonal to all circuits
    for T in allPM(n)
        orth = true
        for C in circuits
            if !orthogonal(C, T)
                orth = false
                break
            end
        end
        if orth == true
            append!(topes, [T])
        end
    end
    M.topes = topes
    return topes
end

# find a list of circuits on support orthogonal to topes
# if topes actually satisfies tope axioms, this will be one pair of opposites
# otherwise, will have more pairs
function findCircuit(topes, support)
    n = length(support)
    m = length(topes[1])
    restriction = [T[support] for T in topes]

    patterns_found = zeros(2^(n))
    count = 0
    i = 1
    while count < 2^(n) && i <= length(restriction)
        label = binaryLabel(restriction[i])
        if patterns_found[label+1]==0
            count += 1
            patterns_found[label+1] = 1
        end
        i+=1
    end
    circuits = SignVector[]
    if count < 2^n
        circuit_labels = findall(iszero, patterns_found)
        for label in circuit_labels
            circuit = SignVector(zeros(m))
            circuit[support] = labelToSignVec(label-1,n)
            append!(circuits,[circuit])
        end
    end
    return circuits
end


# M an oriented matroid with topes filled in
# need to specify a rank, currently only works w/ uniform matroids
# computes circuits from the topes and fills them in
function topesToCircuits!(M; uniform = true)
    rank = M.r
    topes = M.topes
    n = length(topes[1])
    circuits = Dict{Array{Int}, Array{SignVector}}()
    if uniform
        for sigma in powerset(collect(1:n), rank+1, rank+1)
            circuits[sigma] = findCircuit(topes, sigma)
        end
    end
    M.circuits = circuits
    return circuits
end

function chirotopeToTopes!(M)
    χ = M.chirotope
    n = M.n
    r = M.r
    M.topes = topes(χ, n, r)
end

function chirotopeToCocircuits!(M)
    χ = M.chirotope
    n = M.n
    r = M.r
    setCoCircuits!(M, cocircuits(χ, n))
end
