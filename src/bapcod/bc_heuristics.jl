############# Heuristics ################
function ifloor(a)
  r = floor(a + a*1e-10 + 1e-6)
  (r < a - 1 + 1e-6) && (r+=1)
  return r
end

function delete_negative_profits!(profits, weights)
  l = length(profits)
  corresp = Vector{Integer}()
  j = 1
  for i in 1:l
    if profits[j] <= 0
      deleteat!(profits, j)
      deleteat!(weights, j)
    else
      push!(corresp, i)
      j += 1
    end
  end
  return corresp
end

function pisinger_minknap(profits::Vector{Cint}, weights::Vector{Cint},
  capacity::Cint)
  nbitems = length(profits)
  solution = fill(Cint(0), nbitems)
  obj = @pisinger_knp_ccall("callminknap", Clong, (Cint, Ptr{Cint},
                                                   Ptr{Cint}, Ptr{Cint}, Cint),
                                                   nbitems, profits, weights,
                                                   solution, capacity)
  return (obj, solution)
end

pisinger_minknap(profits::Vector{Int64}, weights::Vector{Int64},
  capacity::Integer) = pisinger_minknap(convert(Vector{Cint}, profits),
  convert(Vector{Cint}, weights), convert(Cint, capcity))

function pisinger_minknap(profits::Vector{Cdouble}, weights::Vector{Cint},
  capacity::Cint)
  nbitems = length(profits)
  solution = fill(Cint(0), nbitems)
  @pisinger_knp_ccall("calldminknap", Cdouble, (Cint, Ptr{Cdouble},
                                                   Ptr{Cint}, Ptr{Cint}, Cint),
                                                   nbitems, profits, weights,
                                                   solution, capacity)
  obj = mapreduce(i -> solution[i]*profits[i], +, 1:nbitems)
  return (obj, solution)
end

pisinger_minknap(profits::Vector{Cdouble}, weights::Vector{Int64},
  capacity::Integer) = pisinger_minknap(profits, convert(Vector{Cint}, weights),
  convert(Cint, capacity))

function pisinger_bouknap(profits::Vector{Cint}, weights::Vector{Cint},
  ub::Vector{Cint}, capacity::Cint)
  obj = @pisinger_knp_ccall("callbouknap", Clong, (Cint, Ptr{Cint}, Ptr{Cint},
                                                   Ptr{Cint}, Ptr{Cint}, Cint),
                                                   nbitems, profits, weights,
                                                   ub, solution, capacity)
  return (obj, solution)
end

pisinger_bouknap(profits::Vector{Int64}, weights::Vector{Int64},
  ub::Vector{Int64}, capacity::Integer) = pisinger_bouknap(
  convert(Vector{Cint}, profits), convert(Vector{Cint}, weights),
  convert(Vector{Cint}, ub), convert(Cint, capacity))

function pisinger_bouknap(profits::Vector{Cdouble}, weights::Vector{Cint},
  ub::Vector{Cint}, capacity::Cint)
  @pisinger_knp_ccall("calldbouknap", Cdouble, (Cint, Ptr{Cdouble}, Ptr{Cint},
                                                Ptr{Cint}, Ptr{Cint}, Cint),
                                                nbitems, profits, weights,
                                                ub, solution, capacity)
  obj = mapreduce(i -> solution[i]*profits[i], +, 1:nbitems)
  return (obj, solution)
end

pisinger_bouknap(profits::Vector{Cdouble}, weights::Vector{Int64},
  ub::Vector{Int64}, capacity::Integer) = pisinger_bouknap(profits,
  convert(Vector{Cint}, weights), convert(Vector{Cint}, ub),
  convert(Cint, capacity))

function solvewith_pisinger_minknap(profits::Union{Vector{Int32}, Vector{Int64}, Vector{Cdouble}},
  weights::Union{Vector{Int32}, Vector{Int64}}, capacity::Integer)
  nbitems = length(profits)
  solution = fill(Cint(0), nbitems)
  corresp = delete_negative_profits!(profits, weights)
  if !isempty(profits)
    (obj, part_solution) = pisinger_minknap(profits, weights, capacity)
    for (i,j) in enumerate(corresp)
      solution[j] = part_solution[i]
    end
    return (obj, solution)
  end
  return (0, solution)
end

function solvewith_pisinger_bouknap(profits::Union{Vector{Int32}, Vector{Int64}, Vector{Cdouble}},
  weights::Union{Vector{Int32}, Vector{Int64}}, ub::Union{Vector{Int32}, Vector{Int64}}, capacity::Integer)
  nbitems = length(profits)
  solution = fill(Cint(0), nbitems)
  corresp = delete_negative_profits!(profits, weights)
  if !isempty(profits)
    (obj, part_solution) = pisinger_bouknap(profits, weights, capacity)
    for (i,j) in enumerate(corresp)
      solution[j] = part_solution[i]
    end
    return (obj, solution)
  end
  return (0, solution)
end
