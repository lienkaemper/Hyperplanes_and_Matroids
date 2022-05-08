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
