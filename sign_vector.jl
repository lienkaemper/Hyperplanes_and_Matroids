#---------------Sign Vector Stuff------------------------#

struct SignVector
    vec::Array{Int64}
end

function SignVector(str::String)
    vec = zeros(length(str))
    for i in 1:length(str)
        if str[i] == '+'
            vec[i] = 1
        elseif str[i] == '-'
            vec[i] = -1
        end
    end
    return SignVector(vec)
end

#print sign vectors good
function Base.show(io::IO, X::SignVector)
    str = ""
    V = X.vec
    for i in V
        if i == -1
            str = str*"-"
        elseif i == 1
            str = str*"+"
        else i == 0
            str = str*"0"
        end
    end
    print(io, str)
end

Base.getindex(X::SignVector, i::Int)= X.vec[i]
Base.getindex(X::SignVector, I::Array{Int})= SignVector(X.vec[I])
Base.getindex(X::SignVector, u::UnitRange{Int})= SignVector(X.vec[u])
Base.length(X::SignVector)= length(X.vec)
Base.:-(X::SignVector)= SignVector(-X.vec)
Base.copy(X::SignVector) = SignVector(copy(X.vec))
Base.:(==)(X::SignVector, Y::SignVector) = Int.(sign.(X.vec)) == Int.(sign.(Y.vec))
Base.isequal(X::SignVector, Y::SignVector) = X.vec == Y.vec
Base.last(X::SignVector) = last(X.vec)

function Base.setindex!(X::SignVector, val::Real, i::Int)
    X.vec[i] = Int(sign(val))
end

function Base.setindex!(X::SignVector, val::Array{Real}, i::Array{Int})
    X.vec[i] = Int(sign(val))
end

function Base.setindex!(X::SignVector, val::SignVector, i::Array{Int})
    X.vec[i] = val.vec
end


#----------------Sign Vector Tools--------------------------#

function positivePart(X::SignVector)
    return(findall(x-> x> 0, X.vec))
end

function negativePart(X::SignVector)
    return(findall(x-> x< 0, X.vec))
end

function support(X::SignVector)
    return(findall(x-> x!= 0, X.vec))
end

function zeroPart(X::SignVector)
    return(findall(x-> x== 0, X.vec))
end

function orthogonal(X::SignVector, Y ::SignVector)
    x = X.vec
    y = Y.vec
    n = length(x)

    #first, check if supports are disjoint
    if sum([abs(x[i]*y[i]) for i = 1:n]) == 0
        return true
    end

    #otherwise, check that vectors are neither equal nor opposite
    diff = 0
    same = 0
    for i = 1:n
        if x[i]*y[i] == 1
            same = 1
        elseif x[i]*y[i] == -1
            diff = 1
        end
    end
    orth = diff*same
    return Bool(orth)
end


function composition(X::SignVector, Y::SignVector)
    result = copy(X)
    for i = 1:length(X)
        if X[i]==0
            result[i] = Y[i]
        end
    end
    return result

end

function composition(X::SignVector, Y::SignVector)
    result = copy(X)
    for i = 1:length(X)
        if X[i]==0
            result[i] = Y[i]
        end
    end
    return result
end

function composition(X::SignVector...)
    result = X[1]
    n_inputs = length(X)
    for i = 2:n_inputs
        result = composition(result, X[i])
    end
    return result
end

function sep(X::SignVector,Y::SignVector)
    result = []
    for i =1:length(X)
        if X[i]*Y[i] == -1
            append!(result,i)
        end
    end
    return result
end

#labels a sign vector of length n w/ a number between 0, 2^n-1
function binaryLabel(X::SignVector)
    n = length(X)
    return sum([2^(n-i)*Int64((X[i]+1)/2) for i=1:length(X)])
end

#inverse of binaryLabel
#k: binary label of sign vec
#n: desired length of sign vec
function labelToSignVec(k,n)
    v = last(bitstring(k), n)
    n = length(v)
    result = zeros(Int64, n)
    for i = 1:n
        if v[i] == '0'
            result[i] = -1
        elseif v[i]=='1'
            result[i] = 1
        end
    end
    return SignVector(result)
end



# returns a vector containing all 2^n vectors of length n w/ entries ±1
function allPM(n)
    vecs = [labelToSignVec(i,n) for i =0:2^n-1]
    return vecs
end
