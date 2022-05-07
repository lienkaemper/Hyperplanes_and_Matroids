#you need to add this package if you don't have it
using Combinatorics #we need this for powerset() function

include("sign_vector.jl")


#---------------Oriented Matroid Basics----------------------#

mutable struct OrientedMatroid
    topes::Array{SignVector}
    circuits::Dict{Array{Int}, Array{SignVector}}
    OrientedMatroid() = new()
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

#---------------------Matroid Completion----------------------------#
# makes sure topes actually satisfy symmetry
# works on both oriented matroids and partial oriented matroids
function symmetrizeTopes!(M::OrientedMatroid)
    topes = M.topes
    new_topes = SignVector[]
    for tope in topes
        append!(new_topes, [tope, -tope])
    end
    M.topes = unique(x -> x.vec, new_topes)
end

function symmetrizeTopes(topes::Array{SignVector})
    new_topes = SignVector[]
    for tope in topes
        append!(new_topes, [tope, -tope])
    end
    return unique(x -> x.vec, new_topes)
end

# looks for completions of a partial oriented matroid M
# need to compute circuits first
function findCompletions(M::OrientedMatroid; lazy = true)
    C = M.potential_circuits
    found_matroid = false
    d = []
    for supp in C
        append!(d, Int(length(supp)/2))
    end
    indices = collect(Iterators.product([collect(1:d_i) for d_i in d]...))
    n_i = length(indices[1])
    valid_sets = []
    for index in indices #choose a pair of opposite circuits for each support
        test_set = SignVector[]
        for i = 1:n_i
            c_i = C[i][index[i]]
            append!(test_set, [c_i, -c_i])
        end
        if checkCircuitElimination(test_set, verbose = false)
            append!(valid_sets, [test_set])
            if lazy
                return test_set
            end
        end
    end
    return valid_sets
end



function vcDim(M)
    topes = M.topes
    n = length(topes[1])
    sh_C = [[[i] for i = 1:n],]
    for d = 1:n
        sh_C_d = []
        for sigma in sh_C[d]
            j = maximum(sigma)
            for i = (j+1):n
                tau = vcat(sigma, i)
                if findCircuit(topes, tau) == []
                    append!(sh_C_d, [tau])
                end
            end
        end
        if sh_C_d == []
            break
        end
        append!(sh_C, [sh_C_d])
    end
    return length(sh_C)
end
