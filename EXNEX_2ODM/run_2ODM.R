# compare naive exnex and logodds exnex:

# common setup:
J           = 5
n           = rep(30, J)
n.exch      = 2
nex.mean    = 0.2
nex.prec    = 0.01
ex.mu.mean  = c(0, 0.5)
ex.mu.prec  = c(0.01, 0.01)
ex.tau.prec = c(1, 1)
tau         = 0.8
num.run     = 1000

# # uncomment these 4 lines:
# H0 = rep(0.2, J)
# H1 = rep(0.4, J)
# p0 = H0
# scenario = "1-GN"

source("src_compare_exnex_2ODM.R")

