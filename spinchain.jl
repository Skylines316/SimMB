using Combinatorics
using LinearAlgebra
using Plots
using Statistics
# using Distributed
# using Base.Threads

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

function dotProduct(interm, base)
  intermBin = interm
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


function Hamil(base_i, base, nu, k, n)
  initial = decimal_to_binary(base[base_i],n)
  Ham_pm_k = zeros(Complex{Float64}, length(base), n)
  for i in 1:length(initial)
    base_f, d = dotProduct(initial, base)
    for j in k
      Ham_pm_k[base_f, j+1] += (-1)^(initial[i]+1)*(-1)^(initial[mod1(i+1,n)]+1)*1/2 * exp(im*2*pi*j*d/n) * sqrt(nu[base_i]) / sqrt(nu[base_f])
    end
    interm = copy(initial)
    if initial[i] == 1 && initial[mod1((i+1),n)] == 0
      interm[i] = 0
      interm[mod1((i+1),n)] = 1
      base_f, d = dotProduct(interm, base)
      for j in k
        Ham_pm_k[base_f, j+1] +=  exp(im*2*pi*j*d/n) * sqrt(nu[base_i]) / sqrt(nu[base_f])
      end
    end  
    if initial[i] == 0 && initial[mod1((i+1),n)] == 1
      interm[i] = 1
      interm[mod1((i+1),n)] = 0
      base_f, d = dotProduct(interm, base)
      for j in k
        Ham_pm_k[base_f, j+1] += exp(im*2*pi*j*d/n) * sqrt(nu[base_i]) / sqrt(nu[base_f])
        # print("llego", base_f)
      end
    end
  end
  return Ham_pm_k
end

function Hamiltonian(base, nu, k, n)
  Ham = zeros(Complex{Float64}, length(base), length(base), n)
  for i in 1:length(base)
    Ham[i,:,:] = Hamil(i, base, nu, k[i], n)
  end
  return Ham
end

Nup = 8
Ndown = 6   
Ntot = Ndown + Nup

base, nu, k = genNconsv(Nup,Ndown)
# print(decimal_to_binary(27,6))
# q,r = dotProduct(3, base, Ntot)
# print(hcat(k...)[1,:])

Ham = Hamiltonian(base, nu, k, Ntot)
# print(base)
# display(real(Ham[:,:,1]))

# display(Ham[:,:,1])
# # print(base)
# # display(Ham[:,:,1])
if sum(imag(eigvals(Ham[:,:,1]))) == 0.0
  En = sort(-real(eigvals(Ham[:,:,1])))
  spacings = [En[i+1]-En[i] for i in 1:length(En)-1]
  # spacings_normalized = spacings/mean(spacings)
  # print(En)
  print(sum(En))
  # # println(spacings_normalized)
  # plot(En, 1:length(En)) 
  r = [min(spacings[i+1],spacings[i])/max(spacings[i+1],spacings[i]) for i in 1:length(spacings)-1]
  # # print(r)
  histogram(spacings,normalize=:true)
  print("good")
else 
  print("error")
end
png("spacings")

# print(q,r)
# print(base, nu, k)
# print("Hola Mundo")
