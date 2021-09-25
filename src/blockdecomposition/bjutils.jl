mutable struct BlockDecompositionError <: Exception
  info
  message
end

bjerror(message) = throw(BlockDecompositionError(nothing, message))
bjerror(info,message) = throw(BlockDecompositionError(info, message))
function Base.showerror(io::IO, e::BlockDecompositionError)
  println(io, " ")
  println(io, "BlockDecompositionError : ", e.message)
  if e.info != nothing
    println(io, " ")
    println(io, e.info)
  end
  println(io, " ")
  println(io, "-------------------------------")
end

isa_jumparray(contnr) = isa(contnr, JuMP.JuMPArray)
isa_jumpdict(contnr) = isa(contnr, JuMP.JuMPDict)
isa_jumpcontnr(contnr) = (isa_jumpdict(contnr) || isa_jumparray(contnr))
isa_array(contnr) = isa(contnr, Array)
isa_jumpcstr(contnr) = isa(contnr, JuMP.ConstraintRef)
isa_jumpvar(contnr) = isa(contnr, JuMP.Variable)
contains_in_type(contnr, t1, t2) = contains("$(typeof(contnr))", t1) && contains("$(typeof(contnr))", t2)
contains_jumpcstr(contnr) = contains_in_type(contnr, "JuMP.", "ConstraintRef")
contains_jumpvar(contnr) = contains_in_type(contnr, "JuMP.", "Variable")

name(x::JuMP.JuMPContainer) = Symbol(x.meta[:model].varData[x].name)
name(x::JuMP.Variable) = Symbol(x.m.colNames[x.col])

jumpmodel(x::JuMP.JuMPContainer) = x.meta[:model]
jumpmodel(c::JuMP.JuMPArray) = c.innerArray[end].m
jumpmodel(x::JuMP.Variable) = x.m
jumpmodel(c::JuMP.ConstraintRef) = c.m

# Get the keys of an array of dimension dim
#=mutable struct KeyArrayIterator
  a::Array
  len::Int
  dim::Int
  function KeyArrayIterator(a)
    new(a, length(a), ndims(a))
  end
end

Base.keys(a::Array) = KeyArrayIterator(a)
Base.start(k::KeyArrayIterator) = 1
Base.done(k::KeyArrayIterator, state) = (state > k.len)

function Base.next(key::KeyArrayIterator, state)
  index = fill(0, key.dim)
  size_prod = mapreduce(i -> size(key.a, i), *, [1:key.dim...])
  for k in key.dim:-1:1
    size_prod /= size(key.a, k)
    mod = 0
    size_mod = 1
    for i in k+1:key.dim
      size_mod *= size(key.a, i-1)
      mod += (index[i] - 1) * size_mod
    end
    index[k] = ceil(state / size_prod) - mod
  end
  return (tuple(index...), state+1)
end =#

# Copy from affexpr.jl:88
function assert_isfinite(a::JuMP.AffExpr)
    coeffs = a.coeffs
    for i in 1:length(a.vars)
        isfinite(coeffs[i]) || error("Invalid coefficient $(coeffs[i]) on variable $(a.vars[i])")
    end
end

# Copy from JuMP.jl:457
function verify_ownership(m::JuMP.Model, vec::Vector{JuMP.Variable})
    n = length(vec)
    @inbounds for i in 1:n
        vec[i].m !== m && return false
    end
    return true
end

# Copy from JuMP (solvers.jl:535)
# Convert all the affine constraints into a sparse column-wise
# matrix of coefficients.
function prepConstrMatrix(m::JuMP.Model)

    linconstr = m.linconstr::Vector{JuMP.LinearConstraint}
    numRows = length(linconstr)
    # Calculate the maximum number of nonzeros
    # The actual number may be less because of cancelling or
    # zero-coefficient terms
    nnz = 0
    for c in 1:numRows
        nnz += length(linconstr[c].terms.coeffs)
    end
    # Non-zero row indices
    I = Array{Int}(undef, nnz)
    # Non-zero column indices
    J = Array{Int}(undef, nnz)
    # Non-zero values
    V = Array{Float64}(undef, nnz)

    # Fill it up!
    # Number of nonzeros seen so far
    nnz = 0
    for c in 1:numRows
        # Check that no coefficients are NaN/Inf
        assert_isfinite(linconstr[c].terms)
        coeffs = linconstr[c].terms.coeffs
        vars   = linconstr[c].terms.vars
        # Check that variables belong to this model
        if !verify_ownership(m, vars)
            error("Variable not owned by model present in a constraint")
        end
        # Record all (i,j,v) triplets
        @inbounds for ind in 1:length(coeffs)
            nnz += 1
            I[nnz] = c
            J[nnz] = vars[ind].col
            V[nnz] = coeffs[ind]
        end
    end

    # sparse() handles merging duplicate terms and removing zeros
    return sparse(I,J,V,numRows,m.numCols)
end
