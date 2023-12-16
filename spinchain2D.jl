using Combinatorics
using LinearAlgebra: eigvals
using Plots
# using Statistics

using Base.Threads


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
    cyclicpermut = Set(cyclic_permutations(matrix[i,:,:]))
    push!(nu,length(Set(cyclicpermut)))
    # println(decimal[i+1:end])
    # println(cyclicpermut)
    filtered_permutations = Vector[x for x in decimal[i+1:end] if !(x in cyclicpermut)]
    # filtered_permutations = setdiff(Set(decimal[i+1:end]),Set(cyclicpermut))
    decimal = vcat(decimal[begin:i],filtered_permutations)
    i += 1
    
    # print(filtered_permutations)
  end
  return decimal, nu
end

function genNconsv(Nu,L)
#  N = Nu+Nd
  baseNoT = generate_combinations(Nu,L)
  # for i in baseNoT
  #   print(i)
  # end

  base, nu = Translation(baseNoT)
  # k = kAv(base,nu)
  # baseDec = zeros(Int,length(base))
  return base, nu #, k
end

################################################################################
function decimal_to_binary(base_x,L)
  array = copy(base_x)
  matrix = zeros(Int, L, L)
  
  for i in 1:length(array)
    binary_array = Int[]
    while array[i] > 0
        pushfirst!(binary_array, array[i] % 2)
        array[i] ÷= 2
    end
    matrix[i,:] = vcat(zeros(Int,L-length(binary_array)),binary_array)
  end

  return matrix
end

function dotProduct(interm, base)
  intermBin = interm
  intermCyc = cyclic_permutations(intermBin)
  nu = length(Set(intermCyc))
  # intermCyc = intermCyc[1:nu]
  # print(intermCyc )
  r = 0
  basjav = 0
  distance = 0
  for i in intermCyc
    for j in 1:length(base)
      # println(i, j)
      if base[j] == i
        basjav = j
        distance = r%nu
        return basjav, distance
      end
    end
    r+=1
  end
end

function Hamil(base_i, base, nu, L)
  initial = decimal_to_binary(base[base_i],L)
  Ham_pm_k = zeros(Float64, length(base))
  for i in 0:L^2-1
    base_f, d = dotProduct(initial, base)
    Ham_pm_k[base_f] += ((-1)^(initial[mod1(div(i,L)+1,L),mod1(i%L+1,L)]+1)*(-1)^(initial[mod1(div(i,L)+1,L),mod1(i%L+2,L)]+1) + (-1)^(initial[mod1(div(i,L)+1,L),mod1(i%L+1,L)]+1)*(-1)^(initial[mod1(div(i,L)+2,L),mod1(i%L+1,L)]+1)) * 1/2 * sqrt(nu[base_i]) / sqrt(nu[base_f])
    interm = copy(initial)
    if initial[mod1(div(i,L)+1,L),mod1(i%L+1,L)] == 1 && initial[mod1(div(i,L)+1,L),mod1(i%L+2,L)] == 0
      interm[mod1(div(i,L)+1,L),mod1(i%L+1,L)] = 0
      interm[mod1(div(i,L)+1,L),mod1(i%L+2,L)] = 1
      base_f, d = dotProduct(interm, base)
      Ham_pm_k[base_f] +=  sqrt(nu[base_i]) / sqrt(nu[base_f])
    end  
    interm = copy(initial)
    if initial[mod1(div(i,L)+1,L),mod1(i%L+1,L)] == 1 && initial[mod1(div(i,L)+2,L),mod1(i%L+1,L)] == 0
      interm[mod1(div(i,L)+1,L),mod1(i%L+1,L)] = 0
      interm[mod1(div(i,L)+2,L),mod1(i%L+1,L)] = 1
      # print("qui",initial, base_f)
      base_f, d = dotProduct(interm, base)
      Ham_pm_k[base_f] +=  sqrt(nu[base_i]) / sqrt(nu[base_f])
    end  
    interm = copy(initial)
    if initial[mod1(div(i,L)+1,L),mod1(i%L+1,L)] == 0 && initial[mod1(div(i,L)+1,L),mod1(i%L+2,L)] == 1
      interm[mod1(div(i,L)+1,L),mod1(i%L+1,L)] = 1
      interm[mod1(div(i,L)+1,L),mod1(i%L+2,L)] = 0
      base_f, d = dotProduct(interm, base)
      Ham_pm_k[base_f] +=  sqrt(nu[base_i]) / sqrt(nu[base_f])
    end  
    interm = copy(initial)
    if initial[mod1(div(i,L)+1,L),mod1(i%L+1,L)] == 0 && initial[mod1(div(i,L)+2,L),mod1(i%L+1,L)] == 1
      interm[mod1(div(i,L)+1,L),mod1(i%L+1,L)] = 1
      interm[mod1(div(i,L)+2,L),mod1(i%L+1,L)] = 0
      base_f, d = dotProduct(interm, base)
      Ham_pm_k[base_f] +=  sqrt(nu[base_i]) / sqrt(nu[base_f])
      # print("llego", base_f)
    end
  end
  return Ham_pm_k
end





function Hamiltonian(base, nu, L)
  Ham = zeros(Float64, length(base), length(base))
  file_path = "./Hamiltonians/H2D_u14_d11/L_5.txt"
  file = open(file_path, "w")
  @threads for i in 1:length(base)
    Ham[i,:] = Hamil(i, base, nu, L)
  end
  for i in 1:length(base)
    println(file, join(real(Ham[i,:]), ","))
  end
  close(file)
  return Ham
end



L = 5   #I'm considering L^2 lattice
Nup = 14
base, nu = genNconsv(Nup,L)
Ham = Hamiltonian(base, nu, L)

println("computing eigenvalues...")

En = sort(-eigvals(Ham))
spacings = [En[i+1]-En[i] for i in 1:length(En)-1]
r = [min(spacings[i+1],spacings[i])/max(spacings[i+1],spacings[i]) for i in 1:length(spacings)-1]
# # print(r)
histogram(En,normalize=:true)
png("energy2DL_25Nup_14")
histogram(spacings,normalize=:true)
png("spacings2DL_25Nup_14")
histogram(r,normalize=:true)
png("ratio2DL_25Nup_14")
println("good")



# display(Ham)
# print("base T invariant= ",length(base))
# println(base)
# baseM = decimal_to_binary(base[4],L)
# println(base)
# base_f,d = dotProduct(baseM,base)
# println(base_f)
# display(baseM)
# for k in 0:L^2-1
#   println([mod1(div(k,L)+1,L),mod1(k%L+1,L)],[mod1(div(k,L)+2,L),mod1(k%L+1,L)])
# end
# println(cyclic_permutations(baseM))
# print(cyclic_permutations(x[1,:,:]))
# println("nu= ",nu)
# display(size(base)[1])