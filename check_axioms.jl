
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

function checkSymmetry(circuit_dict::Dict)
    violated = false
    for (key,value) in circuit_dict
        if !checkSymmetry(value)
            violated = true
            println("first violation: ",key, value)
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

# checks incomparability axiom
# returns true if support(X) ⊆ support(Y) ⟹ X = ±Y ∀ X, Y ∈ circuits
function checkIncomparability(circuits; verbose = true)
    violated = false
    for X in keys(circuits), Y in  keys(circuits)
        if X ⊊ Y
            violated = true
            if verbose
                println("first violation: ", X, Y)
            end
            break
        end
    end
    return !violated
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
    return checkSymmetry(circuits)&&checkIncomparability(circuits)&& checkModularCircuitElimination(circuits)
end

#returns true if circuits of M satsify all circuit axions
function checkCoCircuitAxioms(M::OrientedMatroid)
    cocircuits = M.cocircuits
    n = M.n
    return checkSymmetry(cocircuits)&&checkIncomparability(cocircuits)&&checkModularCircuitElimination(cocircuits, n)
end
