
using Combinatorics
using LinearAlgebra
using StatsBase

#rows: hyperplane normals
function randomHyperplanes(n, r)
    return rand(n, r) - (1/2)*ones(n,r)
end

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

# compute the cocircuits of a hyperplane arrangement
# function cocircuits(ha)
#     (n, r) = size(ha)
#     m = n - r + 1
#     covectors = []
#     for sigma in powerset(collect(1:n), m, m)
#         covector = zeros(n)
#         tau = setdiff(collect(1:n), sigma)
#         e = sigma[1]
#         covector[e] = 1
#         X = vcat([e], tau)
#         s1 = sign(det((ha[X, :])))
#         for f in sigma
#             Y = vcat([f], tau)
#             s2 = sign(det((ha[Y, :])))
#             covector[f] = s1*s2
#         end
#         append!(covectors, SignVector.([covector, -covector]))
#     end
#     return covectors
# end

# compute the cocircuits of a hyperplane arrangement
# function cocircuits_2(ha)
#     n = size(ha,1)
#     cocircuit_supports = [(σ["gen"],σ["flat"], setdiff(collect(1:n), σ["flat"])) for σ in coatoms(ha)]
#     cocircuits = []
#     for (δ, σ, τ) in cocircuit_supports
#         println( (δ, σ, τ))
#         signature = SignVector(zeros(n))
#         i = pop!(τ)
#         println(τ)
#         signature[i] = 1
#         for j in τ
#             a = [δ ;[j]]
#             b = [δ ;[i]]
#             signature[j] = sign(Int((det(ha[a , :]) * det(ha[b, :]))))
#             println("prod", Int(det(ha[a , :]) * det(ha[b, :])))
#             println(j, signature[j])
#             if signature[j]== 0
#                 print("this shouldn't happen")
#             end
#         end
#         println(signature)
#         append!(cocircuits, [signature, -signature])
#     end
#     return cocircuits
# end

function basic_cocircuit(χ, σ, e,n)
    #n = size(ha,1)
    signature =  SignVector(zeros(n))
    δ = setdiff(σ , [e])
    a = [δ ;[e]]
    signature[e] = χ[σ] * comp_sign(δ,e)
    for f in setdiff(collect(1:n) ,σ)
        b = [δ ;[f]]
        signature[f] = χ[sort(b)] * comp_sign(δ, f)
        #signature[f] = det(ha[b , :])
    end
    return signature
end

function cocircuits_3(χ, n)
    #n = size(ha,1)
    #χ = chirotope(ha)
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

function topes(ha)
    topes = []
    χ = chirotope(ha)
    r = rank(ha)
    for (σ, sgn) ∈ χ
        if sgn != 0
            for α ∈ allPM(r)
                for i = 1:r
                end
                T_σα = composition([α[i] * basic_cocircuit(ha, σ, σ[i]) for i = 1:r]...)
                if T_σα ∉ topes
                    push!(topes, T_σα)
                end
            end
        end
    end
    return topes
end


#need to remove duplicates, modify to take a dictionary
# or maybe just do something completely different, idk
function topesFromCocircuits(cocircuits)
    topes = Vector{SignVector}()
    for cocircuit in cocircuits
        zero_indices = zeroPart(cocircuit)
        for fill in allPM(length(zero_indices))
            tope = copy(cocircuit)
            tope[zero_indices] = fill
            if !(tope in topes)
                push!(topes, tope)
            end
        end
    end
    return unique(topes)
end

function orientedMatroid(ha)
    M = OrientedMatroid()
    c = cocircuits(ha)
    T = topesFromCocircuits(c)
    setTopes!(M, T)
    return M
end

function coatoms(ha)
    χ = chirotope(ha)
    result = []
    for σ in keys(χ)
        if χ[σ] != 0
            for i in σ
                this_flat =  complete_flat(ha, setdiff(σ, [i]))
                if this_flat["flat"] ∉ [entry["flat"] for entry in result]
                    push!(result, this_flat)
                end
            end
        end
    end
    return result
end

function complete_flat(ha, σ)
    n = size(ha, 1)
    r = rank(ha[σ, :])
    flat = deepcopy(σ)
    for i in setdiff(collect(1:n), σ)
        candidate = [σ; [i]]
        if rank(ha[candidate, :]) == r
            push!(flat, i)
        end
    end

    return Dict("gen" => σ, "flat" => sort!(flat))
end

function comp_sign(σ, e)
    r = length(σ)
    current = ((r + 1) %2)*2 -1
    original = (minimum([findall(σ .> e); r+1]) % 2)*2 -1
    return current * original
end
