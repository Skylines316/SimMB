#import "@preview/physica:0.8.1": *

We have the following Hamiltonian in XXX-Heisenberg model in one dimension with periodic boundary conditions.

$
hat(H) = -J sum_(j=1)^(N)(hat(S)_j^x hat(S)_(j+1)^x + hat(S)_j^y hat(S)_(j+1)^y + hat(S)_j^z hat(S)_(j+1)^z)
$

if we compute the commutator of $hat(H)$ and spin vector operator $arrow(S) = (hat(S)^x , hat(S)^y , hat(S)^z)$, this is

$
[hat(H),arrow(S)] = [-J sum_(j=1)^(N)(hat(S)_j^x hat(S)_(j+1)^x + hat(S)_j^y hat(S)_(j+1)^y + hat(S)_j^z hat(S)_(j+1)^z),(hat(S)^x , hat(S)^y , hat(S)^z)] = -J sum_(j=1)^(N) [hat(S)_j^x hat(S)_(j+1)^x + hat(S)_j^y hat(S)_(j+1)^y + hat(S)_j^z hat(S)_(j+1)^z,(hat(S)^x , hat(S)^y , hat(S)^z)] \
= -J sum_(j=1)^(N) ([hat(S)_j^y hat(S)_(j+1)^y + hat(S)_j^z hat(S)_(j+1)^z, hat(S)^x], [hat(S)_j^x hat(S)_(j+1)^x + hat(S)_j^z hat(S)_(j+1)^z, hat(S)^y], [hat(S)_j^x hat(S)_(j+1)^x + hat(S)_j^y hat(S)_(j+1)^y, hat(S)^z]) = arrow(0)
$

so we can use eigenvectors of any component of $arrow(S)$. 

For preference let use $S_z$ then the we can use the basis ${ket(arrow.t),ket(arrow.b)}$. Indeed the number of $N_(arrow.t)$ and $N_(arrow.b)$ is fixed. Also the total number $N = N_(arrow.t) + N_(arrow.b)$ is fixed.

In this sense it's useful to express $hat(S)_j^x hat(S)_(j+1)^x + hat(S)_j^y hat(S)_(j+1)^y$ in terms of operators $hat(S)_(+)$ and $hat(S)_(-)$ defined by 

$
hat(S)^+ = hat(S)^x + i hat(S)^y \ hat(S)^- = hat(S)^x - i hat(S)^y 
$

and the operators $hat(S)^x$ and $hat(S)^y$ expressed in terms of these operators

$
hat(S)^x = 1/2(hat(S)^+ + hat(S)^-) \ hat(S)^y = i/2(hat(S)^- - hat(S)^+) 
$

then

$
hat(S)_j^x hat(S)_(j+1)^x = (1/2(hat(S)^+_j + hat(S)_j^-))(1/2(hat(S)_(j+1)^+ + hat(S)_(j+1)^-))=1/4(hat(S)^+_j hat(S)_(j+1)^+ + hat(S)^+_j hat(S)^-_(j+1) + hat(S)^-_j hat(S)^+_(j+1) + hat(S)^-_j hat(S)_(j+1)^-) \ 
hat(S)_j^y hat(S)_(j+1)^y = (i/2(hat(S)^-_j - hat(S)_j^+))(i/2(hat(S)_(j+1)^- - hat(S)_(j+1)^+))=-1/4(hat(S)^-_j hat(S)_(j+1)^- - hat(S)^-_j hat(S)^+_(j+1) - hat(S)^+_j hat(S)^-_(j+1) + hat(S)^+_j hat(S)_(j+1)^+)
$

so

$
hat(S)_j^x hat(S)_(j+1)^x + hat(S)_j^y hat(S)_(j+1)^y = 1/2(hat(S)^+_j hat(S)^-_(j+1) + hat(S)^-_j hat(S)^+_(j+1))
$

and the hamiltonian now is

$
hat(H) = -J/2 sum_(j=1)^(N)(hat(S)^+_j hat(S)^-_(j+1) + hat(S)^-_j hat(S)^+_(j+1) + 2hat(S)_j^z hat(S)_(j+1)^z)
$

we know that

$
hat(S)^+ ket(arrow.t) = 0 #h(1cm) hat(S)^+ ket(arrow.b) = ket(arrow.t) #h(1cm) hat(S)^- ket(arrow.t) = ket(arrow.b) #h(1cm) hat(S)^- ket(arrow.b) = 0 #h(1cm) hat(S)^z ket(arrow.t) = 1/2 ket(arrow.t) #h(1cm) hat(S)^z ket(arrow.b) = -1/2 ket(arrow.b)
$

= Construction of the basis
Now we will explote that $N_(arrow.t)$ and $N_(arrow.b)$ are constants. Also we explote translational invariance and periodic boundary conditions.

So we will have a set of representative states 
$
ket(r_k) = hat(P)_k ket(r)/sqrt(mel(r,hat(P)_k,r)) 
$

where a is the lattice parameter. Since we are using this basis each representative state we are using is characterized by its momenta, in other words the momenta is a quantum number well defined. Also we should notice

$
hat(P)_k hat(P)_k = hat(P)_k \
[hat(P)_k , hat(P)_k'] = 0
$

And we will apply $hat(H)$ to these states

$
bra(r_q)hat(H)ket(r_k) = (bra(r_q)hat(H)hat(P)_k ket(r))/sqrt(mel(r,hat(P)_k,r)) = (bra(r_q)hat(P)_k hat(H) ket(r))/sqrt(mel(r,hat(P)_k,r)) = bra(r_q)hat(P)_k (-J/2 sum_(j=1)^(N)(hat(S)^+_j hat(S)^-_(j+1) + hat(S)^-_j hat(S)^+_(j+1) + 2hat(S)_j^z hat(S)_(j+1)^z)) hat(c)^dagger_(N, sigma) dots.h.c hat(c)^dagger_(j, sigma) dots.h.c hat(c)^dagger_(1, sigma) ket(0)/sqrt(mel(r,hat(P)_k,r)) \ 
= -J/2 sum_(j=1)^(N)(bra(r_q)hat(P)_k (hat(S)^+_j hat(S)^-_(j+1) + hat(S)^-_j hat(S)^+_(j+1) + 2hat(S)_j^z hat(S)_(j+1)^z) hat(c)^dagger_(N, sigma) dots.h.c hat(c)^dagger_(j, sigma) dots.h.c hat(c)^dagger_(1, sigma) ket(0))/sqrt(mel(r,hat(P)_k,r))
$

Since we know that

$
[hat(S)^z_l ,hat(c)^dagger_(m,sigma)] = (plus.minus)1/2 delta_(l m) hat(c)^dagger_(m,sigma) #h(1cm) +(-) : arrow.t (arrow.b)\
[hat(S)^+_l ,hat(c)^dagger_(m,sigma)] = 1/2 delta_(l m) hat(c)^dagger_(m,macron(sigma)) delta_(macron(sigma),arrow.t) \
[hat(S)^-_l ,hat(c)^dagger_(m,sigma)] = 1/2 delta_(l m) hat(c)^dagger_(m,macron(sigma)) delta_(macron(sigma),arrow.b) 
$

so for each term in the hamiltonian we will have

$
hat(S)_j^z hat(S)_(j+1)^z ket(r)= hat(c)^dagger_(N, sigma) dots.h.c hat(S)_j^z hat(S)_(j+1)^z hat(c)^dagger_(j, sigma)hat(c)^dagger_(j+1, sigma) dots.h.c hat(c)^dagger_(1, sigma) ket(0) = hat(c)^dagger_(N, sigma) dots.h.c hat(S)_j^z hat(c)^dagger_(j, sigma) hat(S)_(j+1)^z   hat(c)^dagger_(j+1, sigma) dots.h.c hat(c)^dagger_(1, sigma) ket(0) \
= hat(c)^dagger_(N, sigma) dots.h.c (1/2  hat(c)^dagger_(j,sigma) + hat(c)^dagger_(j,sigma) hat(S)^z_j ) (1/2  hat(c)^dagger_(j+1,sigma) + hat(c)^dagger_(j+1,sigma) hat(S)^z_(j+1) ) dots.h.c hat(c)^dagger_(1, sigma) ket(0) = hat(c)^dagger_(N, sigma) dots.h.c (1/2  hat(c)^dagger_(j,sigma) + hat(c)^dagger_(j,sigma) hat(S)^z_j ) (1/2  hat(c)^dagger_(j+1,sigma)) dots.h.c hat(c)^dagger_(1, sigma) ket(0) \
= 1/4 hat(c)^dagger_(N, sigma) dots.h.c hat(c)^dagger_(j,sigma) hat(c)^dagger_(j+1,sigma) dots.h.c hat(c)^dagger_(1, sigma) ket(0) = 1/4 ket(r)
$

$
hat(S)_j^+ hat(S)_(j+1)^- ket(r)= hat(c)^dagger_(N, sigma) dots.h.c hat(S)_j^+ hat(S)_(j+1)^- hat(c)^dagger_(j, sigma)hat(c)^dagger_(j+1, sigma) dots.h.c hat(c)^dagger_(1, sigma) ket(0) = hat(c)^dagger_(N, sigma) dots.h.c hat(S)_j^+ hat(c)^dagger_(j, sigma) hat(S)_(j+1)^- hat(c)^dagger_(j+1, sigma) dots.h.c hat(c)^dagger_(1, sigma) ket(0) \
= hat(c)^dagger_(N, sigma) dots.h.c (1/2  hat(c)^dagger_(j,macron(sigma)) delta_(macron(sigma),arrow.t) + hat(c)^dagger_(j,sigma) hat(S)^+_j ) (1/2  hat(c)^dagger_(j+1,macron(sigma)) delta_(macron(sigma),arrow.b) + hat(c)^dagger_(j+1,sigma) hat(S)^-_(j+1) ) dots.h.c hat(c)^dagger_(1, sigma) ket(0) = hat(c)^dagger_(N, sigma) dots.h.c (1/2  hat(c)^dagger_(j,macron(sigma)) delta_(macron(sigma),arrow.t) + hat(c)^dagger_(j,sigma) hat(S)^+_j ) (1/2  hat(c)^dagger_(j+1,macron(sigma)) delta_(macron(sigma),arrow.b)) dots.h.c hat(c)^dagger_(1, sigma) ket(0) \
= 1/4 delta_(macron(sigma)_(j+1),arrow.b) delta_(macron(sigma)_j,arrow.t) hat(c)^dagger_(N, sigma) dots.h.c hat(c)^dagger_(j,macron(sigma)) hat(c)^dagger_(j+1,macron(sigma)) dots.h.c hat(c)^dagger_(1, sigma) ket(0) 
$

$
hat(S)_j^- hat(S)_(j+1)^+ ket(r) = hat(c)^dagger_(N, sigma) dots.h.c hat(S)_j^- hat(S)_(j+1)^+ hat(c)^dagger_(j, sigma)hat(c)^dagger_(j+1, sigma) dots.h.c hat(c)^dagger_(1, sigma) ket(0) = hat(c)^dagger_(N, sigma) dots.h.c hat(S)_j^- hat(c)^dagger_(j, sigma) hat(S)_(j+1)^+ hat(c)^dagger_(j+1, sigma) dots.h.c hat(c)^dagger_(1, sigma) ket(0) \
= hat(c)^dagger_(N, sigma) dots.h.c (1/2  hat(c)^dagger_(j,macron(sigma)) delta_(macron(sigma),arrow.b) + hat(c)^dagger_(j,sigma) hat(S)^-_j ) (1/2  hat(c)^dagger_(j+1,macron(sigma)) delta_(macron(sigma),arrow.t) + hat(c)^dagger_(j+1,sigma) hat(S)^+_(j+1) ) dots.h.c hat(c)^dagger_(1, sigma) ket(0) = hat(c)^dagger_(N, sigma) dots.h.c (1/2  hat(c)^dagger_(j,macron(sigma)) delta_(macron(sigma),arrow.b) + hat(c)^dagger_(j,sigma) hat(S)^-_j ) (1/2  hat(c)^dagger_(j+1,macron(sigma)) delta_(macron(sigma),arrow.t)) dots.h.c hat(c)^dagger_(1, sigma) ket(0) \
= 1/4 delta_(macron(sigma)_(j+1),arrow.t) delta_(macron(sigma)_j,arrow.b) hat(c)^dagger_(N, sigma) dots.h.c hat(c)^dagger_(j,macron(sigma)) hat(c)^dagger_(j+1,macron(sigma)) dots.h.c hat(c)^dagger_(1, sigma) ket(0) 
$

then we can express the Hamiltonian as

$
bra(r_q)hat(H)ket(r_k) = -J/2 sum_(j=1)^(N)[1/2(bra(r_q)hat(P)_k ket(r))/sqrt(mel(r,hat(P)_k,r)) + 1/4 (bra(r_q)hat(P)_k delta_(macron(sigma)_(j+1),arrow.b) delta_(macron(sigma)_j,arrow.t) hat(c)^dagger_(N, sigma) dots.h.c hat(c)^dagger_(j,macron(sigma)) hat(c)^dagger_(j+1,macron(sigma)) dots.h.c hat(c)^dagger_(1, sigma) ket(0))/sqrt(mel(r,hat(P)_k,r)) \
#h(10cm) + 1/4 (bra(r_q)hat(P)_k delta_(macron(sigma)_(j+1),arrow.t) delta_(macron(sigma)_j,arrow.b) hat(c)^dagger_(N, sigma) dots.h.c hat(c)^dagger_(j,macron(sigma)) hat(c)^dagger_(j+1,macron(sigma)) dots.h.c hat(c)^dagger_(1, sigma) ket(0))/sqrt(mel(r,hat(P)_k,r))]
$

= Some aspects about the basis
In general we choose only representative states, with the projection operator

$
hat(P)_k = 1/L sum_(j=0)^(L-1)exp(i 2 pi k/L j)hat(T)^(j)
$

and if we apply $hat(T)^(j)$ to some state we will obtain a minus sign if the total number of lattice elements in the chain is even.

== Example
$
hat(T)ket(0101) = hat(T) hat(c)^dagger_(4, arrow.b) hat(c)^dagger_(3, arrow.t) hat(c)^dagger_(2, arrow.b) hat(c)^dagger_(1, arrow.t) ket(0) = hat(c)^dagger_(1, arrow.b) hat(c)^dagger_(4, arrow.t) hat(c)^dagger_(3, arrow.b) hat(c)^dagger_(2, arrow.t) ket(0) \
hat(T)ket(0101) = - hat(c)^dagger_(4, arrow.t) hat(c)^dagger_(3, arrow.b) hat(c)^dagger_(2, arrow.t) hat(c)^dagger_(1, arrow.b) ket(0) = -ket(1010)
$

$
hat(T)ket(101101) = hat(T) hat(c)^dagger_(6, arrow.t) hat(c)^dagger_(5, arrow.b) hat(c)^dagger_(4, arrow.t) hat(c)^dagger_(3, arrow.t) hat(c)^dagger_(2, arrow.b) hat(c)^dagger_(1, arrow.t) ket(0) = hat(c)^dagger_(1, arrow.t) hat(c)^dagger_(6, arrow.t) hat(c)^dagger_(5, arrow.b) hat(c)^dagger_(4, arrow.t) hat(c)^dagger_(3, arrow.t) hat(c)^dagger_(2, arrow.b) ket(0) \
hat(T)ket(101101) = - hat(c)^dagger_(6, arrow.t) hat(c)^dagger_(5, arrow.b) hat(c)^dagger_(4, arrow.t) hat(c)^dagger_(3, arrow.t) hat(c)^dagger_(2, arrow.b) hat(c)^dagger_(1, arrow.t) ket(0) = -ket(110110)
$

$
hat(T)ket(101101101) = hat(T) hat(c)^dagger_(9, arrow.t) hat(c)^dagger_(8, arrow.b) hat(c)^dagger_(7, arrow.t) hat(c)^dagger_(6, arrow.t) hat(c)^dagger_(5, arrow.b) hat(c)^dagger_(4, arrow.t) hat(c)^dagger_(3, arrow.t) hat(c)^dagger_(2, arrow.b) hat(c)^dagger_(1, arrow.t) ket(0) = hat(c)^dagger_(1, arrow.t) hat(c)^dagger_(9, arrow.t) hat(c)^dagger_(8, arrow.b) hat(c)^dagger_(7, arrow.t) hat(c)^dagger_(6, arrow.t) hat(c)^dagger_(5, arrow.b) hat(c)^dagger_(4, arrow.t) hat(c)^dagger_(3, arrow.t) hat(c)^dagger_(2, arrow.b) ket(0) \
hat(T)ket(101101101) = hat(c)^dagger_(9, arrow.t) hat(c)^dagger_(8, arrow.b) hat(c)^dagger_(7, arrow.t) hat(c)^dagger_(6, arrow.t) hat(c)^dagger_(5, arrow.b) hat(c)^dagger_(4, arrow.t) hat(c)^dagger_(3, arrow.t) hat(c)^dagger_(2, arrow.b) hat(c)^dagger_(1, arrow.t) ket(0) = ket(110110110)
$

when we apply now this $hat(P)_k$ over the representative vector we obtatin cases where the chain repeat itself whe the number of translations is less than the number of lattice sites. From the previous example we see the three possible cases, when the subset of the chain the lattice set that repeat itself is even and the total lattice is even (1), when subset is odd and total lattice site is even (2) and when the subset is even and the total is odd (3). We can see that it's not possible the case where subset is even and the total is odd.

From this relation also we can compute which elements of $k$ will give is a non vanishing wavefunction. Let see case by case considering the length of the subset is $nu$

+ *Subset even and chain even* \
  The number that the element will repeat itself is $r=L/nu$ and since $nu$ is even the sign of this elements is the same. So if the value of the following sum is enoughly close to zero then k will give us a vanishing wavefunction.
  ```jl
  S = abs(sum([exp(im*2*pi*i*nu*k/n) for i in 0:r-1]))
  ```

+ *Subset odd and chain even* \
  The number that the element will repeat itself is $r=L/nu$ and since $nu$ is odd the sign of this elements is different. So if the value of the following sum is enoughly close to zero then k will give us a vanishing wavefunction.
  ```jl
  S = abs(sum([exp(im*2*pi*i*nu*k/n)*(-1)^r for i in 0:r-1]))
  ```

+ *Subset odd and chain odd* \
  The number that the element will repeat itself is $r=L/nu$ and since $L$ is even the sign of this elements is the same. So if the value of the following sum is enoughly close to zero then k will give us a vanishing wavefunction.
  ```jl
  S = abs(sum([exp(im*2*pi*i*nu*k/n) for i in 0:r-1]))
  ```

= Some aspects about the dot product
We have the following dot value $mel(r_q,hat(H),r_k)$, let focus for a moment in $hat(H)ket(r_k)$

$
hat(H)ket(r_k) = (hat(H) hat(P)_k ket(r_0))/(sqrt(mel(r_0,hat(P)_k,r_0))) =(hat(P)_k hat(H) ket(r_0))/(sqrt(mel(r_0,hat(P)_k,r_0))) = (hat(P)_k (lambda_0) ket(r_(i')))/(sqrt(mel(r_0,hat(P)_k,r_0))) = (hat(P)_k (lambda_0) e^(i 2 pi k/L d) ket(r_i))/(sqrt(mel(r_0,hat(P)_k,r_0))) 
$

where \
#h(1cm) $r_0:$ is the representative of the initial state \
#h(1cm) $r_(i'):$ is a state that belongs to same cycle as the $r_i$ \
#h(1cm) $d:$ is the #text(style: "italic")[distance] from $r_i$ to $r_(i')$ \
#h(1cm) $r_i:$ is the representative of the intermediate state \
#h(1cm) $r_f:$ is the representative of the final state

$
mel(r_q,hat(H),r_k) = (bra(r_f) hat(P)_k hat(P)_k (lambda_0) e^(i 2 pi k/L d) ket(r_i))/(sqrt(mel(r_0,hat(P)_k,r_0))sqrt(mel(r_f,hat(P)_k,r_f))) = lambda_0 e^(i 2 pi k/L d) sqrt(mel(r_f,hat(P)_k,r_f))/sqrt(mel(r_0,hat(P)_k,r_0))
$

in the last step $r_i$ has to be equal to $r_f$ to have non zero values. And the value of $mel(r_l,hat(P)_k,r_l)$ is equal to $1/nu_l$. So we can express 

$
mel(r_q,hat(H),r_k) = lambda_0 e^(i 2 pi k/L d) sqrt(nu_0)/sqrt(nu_f)
$

and since we have a matrix for each k we can represent this as

$ H^(k=0) = mat(
  H_(11), H_(12), ..., H_(1 b);
  H_(21), H_(22), ..., H_(2 b);
  dots.v, dots.v, dots.down, dots.v;
  H_(b 1), H_(b 2), ..., H_(bb);
) $

where $b$ is the total number of representative states we have.
