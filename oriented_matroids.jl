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

function setTopes!(M,topes::Vector)
    M.topes = SignVector.(topes)
    return M
end

function setTopes!(M,topes::String)
    M.topes = SignVector.(String.(strip.(split(topes, ','))))
    return M
end

function setCircuits!(M,circuits::Vector)
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

function setCoCircuits!(M,cocircuits::Vector{SignVector})
    values = cocircuits
    M.cocircuits = Dict{Array{Int}, Array{SignVector}}()
    for cocircuit in values
        supp = support(cocircuit)
        if haskey(M.cocircuits, supp)
            append!(M.cocircuits[supp], [cocircuit])
        else
            M.cocircuits[supp] = [cocircuit]
        end
    end
    return M
end

function setCoCircuits!(M,cocircuits::String)
    values = SignVector.(String.(strip.(split(cocircuits, ','))))
    M.cocircuits = Dict{Array{Int}, Array{SignVector}}()
    for cocircuit in values
        supp = support(cocircuit)
        if haskey(M.cocircuits, supp)
            append!(M.cocircuits[supp], [cocircuit])
        else
            M.cocircuits[supp] = [cocircuit]
        end
    end
    return M
end


#---------------Affine oriented matroids----------------------#

function affineTopes(M; k = M.n)
    aff_topes = []
    actual_planes = setdiff(collect(1:M.n), k)
    for tope in M.topes
        if tope[k] > 0
            push!(aff_topes, tope[actual_planes])
        end
    end
    return aff_topes
end
