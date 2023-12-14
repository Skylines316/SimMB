using Combinatorics

function generate_combinations(N, L)
  M = L^2
  array = zeros(Int, M)

  index_combinations = combinations(1:M, N)

  combinations_array = [[i in indices ? 1 : 0 for i in 1:M] for indices in index_combinations]

  matrix2D = zeros(Int64, length(combinations_array), L, L)
  
  for i in 1:length(combinations_array)
    matrix2D[i,:,:] = transpose(hcat([combinations_array[i][n*L+1:(n+1)*L] for n in 0:L-1]...))
  end

  return matrix2D
end

function cyclic_permutations(x::Matrix)
  L = size(x)[1]
  cyclic = Vector[]
  for i in 0:L-1
    for j in 0:L-1
      push!(cyclic,[parse(Int, join(circshift(x,(i,j))[k,:]), base=2) for k in 1:L])
    end
  end
  return cyclic
end

function Translation(matrix)
  decimal = Vector[]
  for j in 1:size(matrix)[1]
    push!(decimal,[parse(Int, join(matrix[j,k,:]), base=2) for k in 1:L])
  end
  i = 1
  nu = Int[]
  while i < length(decimal)
    cyclicpermut = cyclic_permutations(matrix[i,:,:])
    push!(nu,length(Set(cyclicpermut)))
    filtered_permutations = setdiff(Set(decimal[i+1:end]),Set(cyclicpermut))
    decimal = vcat(decimal[begin:i],collect(filtered_permutations))
    i += 1
    
    # print(filtered_permutations)
  end
  return decimal, nu
end

L=3
Nup = 3
x = generate_combinations(Nup,L)

base, nu = Translation(x)
print("base T invariant= ",length(base))
# print(cyclic_permutations(x[1,:,:]))
print("nu= ",nu)
display(size(x)[1])