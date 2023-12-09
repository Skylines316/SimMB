using Combinatorics
function generate_combinations(N, M)
  array = zeros(Int, M)

  index_combinations = combinations(1:M, N)

  combinations_array = [[i in indices ? 1 : 0 for i in 1:M] for indices in index_combinations]

  return combinations_array
end

function cyclic_permutations(arr)
  n = length(arr)
  return [circshift(arr, i) for i in 0:n-1]
end

function Translation(array)
  i = 1
  while i < length(array)
    cyclicpermut = cyclic_permutations(array[i])
    filtered_permutations = setdiff(Set(array[i+1:end]),Set(cyclicpermut))
    array = vcat(array[begin:i],collect(filtered_permutations))
    i += 1
    # print(filtered_permutations)
  end
  return array
end

function genNconsv(Nu,Nd)
  N = Nu+Nd
  baseNoT = generate_combinations(Nu,N)
  # for i in baseNoT
  #   print(i)
  # end

  baseNoT = collect(baseNoT) 
  base = Translation(baseNoT)
  baseDec = zeros(Int,length(base))
  for j in 1:length(base)
    baseDec[j] = parse(Int, join(collect(base)[j]), base=2)
  end
  return baseDec
end

base = genNconsv(10,10)
print(base)
# print("Hola Mundo")
