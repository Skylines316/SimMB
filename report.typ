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

So we will have a set of representative states $ket(r_k) = hat(P)_k ket(r)/sqrt(mel(r,hat(P)_k,r))$. And we will apply $hat(H)$ to these states

$
bra(r_q)hat(H)ket(r_k) = (bra(r_q)hat(H)hat(P)_k ket(r))/sqrt(mel(r,hat(P)_k,r)) = (bra(r_q)hat(P)_k hat(H) ket(r))/sqrt(mel(r,hat(P)_k,r)) = bra(r_q)hat(P)_k (-J/2 sum_(j=1)^(N)(hat(S)^+_j hat(S)^-_(j+1) + hat(S)^-_j hat(S)^+_(j+1) + 2hat(S)_j^z hat(S)_(j+1)^z)) hat(c)^dagger_(N, sigma) dots.h.c hat(c)^dagger_(j, sigma) dots.h.c hat(c)^dagger_(1, sigma) ket(0)/sqrt(mel(r,hat(P)_k,r)) \ 
= -J/2 sum_(j=1)^(N)(bra(r_q)hat(P)_k (hat(S)^+_j hat(S)^-_(j+1) + hat(S)^-_j hat(S)^+_(j+1) + 2hat(S)_j^z hat(S)_(j+1)^z) hat(c)^dagger_(N, sigma) dots.h.c hat(c)^dagger_(j, sigma) dots.h.c hat(c)^dagger_(1, sigma) ket(0))/sqrt(mel(r,hat(P)_k,r))
$

Since we know that

$
[hat(S)^z_l ,hat(c)^dagger_(m,sigma)] = 1/2 delta_(l m) hat(c)^dagger_(m,sigma) \
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
