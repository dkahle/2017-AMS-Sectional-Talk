##### installing and loading what you need in R
#####################################################

# these are the exact versions of the software i used
# for the presentation
if(!requireNamespace("devtools")) 
  install.packages("devtools")

library("devtools")
install_github("dkahle/mpoly", ref = "670645f")
install_github("coneill-math/m2r", ref = "3cf5e94d")




 
##### demo
#####################################################

library("m2r")

m2("1 + 1")
m2("1.2") # output not parsed
m2_parse(m2("1.2"))

# persistence
m2("a = 1")
m2("a")

# helpers
m2_ls()
m2_exists(c("a", "b"))
m2_getwd()



##### basic commutative algebra
########################################

# rings, ideals, grobner bases
ring("t", "x", "y", "z", coefring = "QQ")
(I <- ideal("t^4 - x", "t^3 - y", "t^2 - z"))
gb(I)
gb("t^4 - x", "t^3 - y", "t^2 - z")

gb(I, code = TRUE)

# radical
ring("x", coefring = "QQ")
I <- ideal("x^2")
radical(I)

# saturation
I <- ideal("(x-1) x (x+1)")
J <- ideal("x")
saturate(I, J) # = (x - 1) (x + 1)

# primary decomposition
QQxyz <- ring("x", "y", "z", coefring = "QQ")
I <- ideal("x z", "y z")
(ideal_list <- primary_decomposition(I))

# dimension
dimension(ideal_list)

# ideal arithmetic
I <- ideal("x", "y")
J <- ideal("z")
I + J
I * J
I == J

# predicate functions
is.radical <- function (ideal) ideal == radical(ideal)
is.radical(I)

# using R workflows
# x %>% f(y)   means   f(x, y)
library("magrittr")
use_ring(QQxyz)
ideal("x z", "y z") %>% 
  primary_decomposition %>% 
  dimension



##### other algorithms
########################################

# factoring
(x <- 2^5 * 3^4 * 5^3 * 7^2 * 11^1)
factor_n(x)

ring("x", "y", coefring = "QQ")
factor_poly("x^4 - y^4") # $factor is mpolyList


# smith normal form
(M <- matrix(c(
   2,  4,   4,
  -6,  6,  12,
  10, -4, -16
), nrow = 3, byrow = TRUE))

(mats <- snf(M)); P <- mats$P; D <- mats$D; Q <- mats$Q

P
str(P)

P %*% M %*% Q                
solve(P) %*% D %*% solve(Q)  

det(P)


##### solving systems via gb
########################################

use_ring(QQxyz)
I <- ideal("x + y + z", "x^2 + y^2 + z^2 - 9", "x^2 + y^2 - z^2")
(grobner_basis <- gb(I))

extract_unipoly <- function(mpolyList) Filter(is.unipoly, mpolyList)[[1]]
which_unipoly <- function(mpolyList) which(sapply(mpolyList, is.unipoly))
solve_gb <- function(gb) {
  poly <- extract_unipoly(gb)
  elim_var <- vars(poly)
  solns <- solve_unipoly(poly, real_only = TRUE)
  
  if(length(gb) == 1) 
    return(structure(t(t(solns)), .Dimnames = list(NULL, elim_var)))
  
  gb <- structure(
    gb[-which_unipoly(gb)[1], drop = FALSE], 
    class = "mpolyList"
  )
  new_systems <- lapply(solns, function(soln) plug(gb, elim_var, soln))
  
  low_solns_list <- lapply(new_systems, solve_gb)
  lower_var_names <- colnames(low_solns_list[[1]])
  
  Map(cbind, solns, low_solns_list) %>% do.call("rbind", .) %>% 
    structure(.Dimnames = list(NULL, c(elim_var, lower_var_names)))
}

# actual solutions are ±(3/sqrt(2))*(1,0,-1) and ±(3/sqrt(2))*(0,1,-1)

# gb solutions
# correct to 14 digits
(solns <- solve_gb(grobner_basis)) %>% 
  structure(.Dimnames = list(paste("Soln", 1:4, ":"), c("z","y","x")))

# naive nonlinear solver solutions
# correct to only 3 digits
resid <- function(v) {
  x <- v[1]; y <- v[2]; z <- v[3]
  (x + y + z)^2 + (x^2 + y^2 + z^2 - 9)^2 + (x^2 + y^2 - z^2)^2
}
optim(c(x = 0, y = 0, z = 0), resid)$par











# what's it doing?
gb(I, code = TRUE)



# standard evaluation
polys <- c("t^4 - x", "t^3 - y", "t^2 - z")
gb_(polys, ring = R)












##### algebraic statistical application
##### independence on a 2x2 table
########################################

ring("p00", "p01", "p10", "p11", coefring = "QQ")
indep_ideal <- ideal(
 "p00 - (p00 + p01) (p00 + p10)",  "p01 - (p00 + p01) (p01 + p11)",    
 "p10 - (p10 + p11) (p00 + p10)",  "p11 - (p10 + p11) (p01 + p11)",    
 "p00 + p01 + p10 + p11 - 1"         
)
gb(indep_ideal)

# first equation is the sum-to-one condition
# second equation equal to p00 p11 - p01 p10 == 0


# the dimension of the corresponding variety is
# important to the practical application of hypothesis tests
dimension(indep_ideal)



##### demo
#####################################################

# in system solution application, we had:
ring("x", "y", "z", coefring = "QQ")
(I <- ideal("x + y + z", "x^2 + y^2 + z^2 - 9", "x^2 + y^2 - z^2"))
gb(I)

# it's faster to do it like this:
ring.("x", "y", "z", coefring = "QQ")
(I. <- ideal.("x + y + z", "x^2 + y^2 + z^2 - 9", "x^2 + y^2 - z^2"))
gb.(I.)




##### macaulay2 in the cloud
#####################################################

stop_m2()
  start_m2(cloud = TRUE)
m2("1+1")







