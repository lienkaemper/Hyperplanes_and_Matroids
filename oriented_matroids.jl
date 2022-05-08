using Combinatorics #we need this for powerset() function
include("sign_vector.jl")


#---------------Oriented Matroid Basics----------------------#

mutable struct OrientedMatroid
    topes::Array{SignVector}
    circuits::Dict{Array{Int}, Array{SignVector}}
    cocircuits::Dict{Array{Int}, Array{SignVector}}
    chirotope::Dict{Array{Int}, Int}
    n::Int
    r::Int
    OrientedMatroid() = new()
end

function setGround!(M::OrientedMatroid, n::Int)
    M.n::n
end


function setRank!(M::OrientedMatroid, r::Int)
    M.r::r
end

function setTopes!(M,topes::Vector{SignVector})
    M.topes = topes
    return M
end

function setTopes!(M,topes::Vector{Vector{Int}})
    M.topes = SignVector.(topes)
    return M
end


function setTopes!(M,topes::String)
    M.topes = SignVector.(String.(strip.(split(topes, ','))))
    return M
end

function setCircuits!(M,circuits::Vector{Vector{Int}})
    values = SignVector.(circuits)
    M.circuits = Dict{Array{Int}, Array{SignVector}}()
    for circuit in values
        supp = support(circuit)
        if haskey(M.circuits, supp)
            append!(M.circuits[supp], [circuit])
        else
            M.circuits[supp] = [circuit]
        end
    end
    return M
end

function setCircuits!(M,circuits::Vector{SignVector})
    values = circuits
    M.circuits = Dict{Array{Int}, Array{SignVector}}()
    for circuit in values
        supp = support(circuit)
        if haskey(M.circuits, supp)
            append!(M.circuits[supp], [circuit])
        else
            M.circuits[supp] = [circuit]
        end
    end
    return M

end

function setCircuits!(M,circuits::String)
    values = SignVector.(String.(strip.(split(circuits, ','))))
    M.circuits = Dict{Array{Int}, Array{SignVector}}()
    for circuit in values
        supp = support(circuit)
        if haskey(M.circuits, supp)
            append!(M.circuits[supp], [circuit])
        else
            M.circuits[supp] = [circuit]
        end
    end
    return M
end



#-----------Translating Between Different Descriptions-----------------#

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
    cocircs = []
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
function topesToCircuits!(M; uniform = true, rank = 0)
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
    M.cocircuits = cocircuits(χ, n)
end


#--------------------Checking Axioms--------------------------#
###NEED TO WRITE A VERSION OF THIS THAT WORKS FOR DICTIONARIES TOO
## SO WE CAN CHECK IT FOR CIRCUITS

# checks symmetry axiom
# returns true if X ∈ sign_vectors ⟹ -X ∈ sign_vectors
function checkSymmetry(sign_vectors::Array{SignVector})
    violated = false
    for X in sign_vectors
        if !(-X in sign_vectors)
            violated = true
            println("first violation: ",X)
            break
        end
    end
    return !violated
end

function modularSupports(support, n)
    supports = []
    complement = setdiff(collect(1:n), support)
    supports = [sort!(union(setdiff(support, i), j)) for i in support for j in complement]
    return supports
end

# checks circuit elimination axiom
#= returns true if:
    for all X, Y ∈ circuits
    if X ≂̸ -Y, then for all e ∈ sep(X,Y) there exists Z ∈ circuits
    such that pos(Z)⊆ pos(X) ∪ pos(Y) ∖ e
              neg(Z)⊆ neg(X) ∪ neg(Y) ∖ e
=#
function checkModularCircuitElimination(circuits, n; verbose = true)
    supports = keys(circuits)
    violated = false
    for σ in supports
        for τ in modularSupports(σ, n)
            for X in circuits[σ], Y in circuits[τ]
                for e in intersect(positivePart(X), negativePart(Y))
                    elim = circuitElimSetFunction(X, Y, e)
                    new_support = sort!(setdiff(union(σ,τ), e))
                    found = false
                    for Z in circuits[new_support]
                        if elim(Z)
                            found = true
                            break
                        end
                    end
                    if !found
                        violated = true
                        if verbose
                            println("first violation: ", X, ", ", Y, ", ", e)
                        end
                        return !violated
                    end
                end
            end
        end
    end
    return !violated
end

#helper for checkCircuitElimination
function circuitElimSetFunction(X::SignVector,Y::SignVector, e::Int )
    new_X = copy(X)
    new_Y = copy(Y)
    new_X[e] = 0
    new_Y[e] = 0
    positives = union(positivePart(new_X), positivePart(new_Y))
    negatives = union(negativePart(new_X), negativePart(new_Y))
    function elim(Z::SignVector)
        return ⊆(positivePart(Z), positives) && ⊆(negativePart(Z), negatives)
    end
    return elim
end

#returns true if circuits of M satsify all circuit axions
function checkCircuitAxioms(M::OrientedMatroid)
    circuits = M.circuits
    return checkSymmetry(circuits)&&checkIncomparability(circuits)&&checkCircuitElimination(circuits)
end
