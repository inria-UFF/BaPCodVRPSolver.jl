

function vec(v::Vector{Float64})
  outstr = ""
  for value in v
    outstr *= "$value  "
  end
  "[  " * outstr * "]"
end

# copied from JuMP print
# Modified to print a JuMPContainer of Vector{Float64}
function val_str(j::JuMP.JuMPArray{Vector{Float64},N}) where N
    m = j.meta[:model] # _getmodel(j)
    data = m.varData[j] # printdata(j)
    out_str = "$(data.name): $N dimensions:\n"
    if isempty(j)
        return out_str * "  (no entries)"
    end

    function val_str_rec(depth, parent_index::Vector{Any}, parent_str::AbstractString)
        # Turn index set into strings
        indexset = data.indexsets[depth]
        index_strs = map(string, indexset)

        # Determine longest index so we can align columns
        max_index_len = 0
        for index_str in index_strs
            max_index_len = max(max_index_len, strwidth(index_str))
        end

        # If have recursed, we need to prepend the parent's index strings
        # accumulated, as well as white space so the alignment works.
        for i = 1:length(index_strs)
            index_strs[i] = parent_str * lpad(index_strs[i],max_index_len," ")
        end

        # Create a string for the number of spaces we need to indent
        indent = " "^(2*(depth-1))

        # Determine the need to recurse
        if depth == N
            # Deepest level
            for i = 1:length(indexset)
                value = length(parent_index) == 0 ?
                            j[indexset[i]] :
                            j[parent_index...,indexset[i]]
                out_str *= indent * "[" * index_strs[i] * "] = $(vec(value))\n"
            end
        else
            # At least one more layer to go
            for i = 1:length(indexset)
                index = indexset[i]
                # Print the ":" version of indices we will recurse over
                out_str *= indent * "[" * index_strs[i] * ",:"^(N-depth) * "]\n"
                val_str_rec(depth+1,
                     length(parent_index) == 0 ? Any[index] : Any[parent_index...,index],
                    index_strs[i] * ",")
            end
        end
    end
    val_str_rec(1,Any[],"")
    return out_str
end
Base.show(io::IO, j::JuMP.JuMPContainer{Vector{Float64}}) = print(io, val_str(j))
