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
  nu = Int[]
  while i < length(array)
    cyclicpermut = cyclic_permutations(array[i])
    push!(nu,length(Set(cyclicpermut)))
    filtered_permutations = setdiff(Set(array[i+1:end]),Set(cyclicpermut))
    array = vcat(array[begin:i],collect(filtered_permutations))
    i += 1
    
    # print(filtered_permutations)
  end
  return array, nu
end

function kAv(rep, nu)
  # k = zeros(Int,length(rep))
  k = Vector{Int64}[]
  n = length(rep[1])
  k_base = [i for i in 0:n-1]
  # print(k)
  for i in 1:length(rep)
    if nu[i] != length(rep[i])
      #compute k available
      k_n = Int[]
      aux = length(rep[i])/nu[i] - 1
      for j in 0:n-1
        if length(rep[i]) % 2 == 0 && nu[i] % 2 == 1
          val = abs(sum([exp(im*2*pi*k*nu[i]*j/n)*(-1)^k for k in 0:aux]))
        else
          val = abs(sum([exp(im*2*pi*k*nu[i]*j/n) for k in 0:aux]))
        end
          # print(val)
        if val < 0.1
          x=4
        else
          push!(k_n, j)
        end
      end 
      # print(k_n,"aqui")
      push!(k,k_n)
    else
      push!(k,k_base)
    end
  end
  return k
end

function genNconsv(Nu,Nd)
  N = Nu+Nd
  baseNoT = generate_combinations(Nu,N)
  # for i in baseNoT
  #   print(i)
  # end

  baseNoT = collect(baseNoT) 
  base, nu = Translation(baseNoT)
  k = kAv(base,nu)
  baseDec = zeros(Int,length(base))
  for j in 1:length(base)
    baseDec[j] = parse(Int, join(base[j]), base=2)
  end
  return baseDec, nu, k
end

function decimal_to_binary(number,n)
  binary_array = Int[]

  while number > 0
      pushfirst!(binary_array, number % 2)
      number ÷= 2
  end

  binary_array = vcat(zeros(Int,n-length(binary_array)),binary_array)

  return binary_array
end

function dotProduct(interm, base,n)
  intermBin = decimal_to_binary(interm,n)
  intermCyc = cyclic_permutations(intermBin)
  nu = length(Set(intermCyc))
  intermCyc = intermCyc[1:nu]
  # print(intermCyc )
  r = 0
  basjav = 0
  distance = 0
  for i in intermCyc
    for j in 1:length(base)
      if base[j] == parse(Int, join(i), base=2)
        basjav = j 
        distance = r%nu
      end
    end
    r+=1
  end 
  return basjav, distance 
end

base, nu, k = genNconsv(2,2)
# print(decimal_to_binary(27,6))
q,r = dotProduct(3, base, 4)
print(q,r)
print(base, nu, k)
# print("Hola Mundo")
