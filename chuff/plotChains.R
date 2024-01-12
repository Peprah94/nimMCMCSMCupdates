library(ggmcmc)

S <- ggs(example1UpdatedModel$samples)
ggs_traceplot(S, family = "c")

S <- ggs( example1ReducedModel$samples)
ggs_traceplot(S, family = "c")


S <- ggs(example2UpdatedModelTrue$samples)
ggs_traceplot(S, family = "alphaPhi")

S <- ggs(example2ReducedModelTrue$samples)
ggs_traceplot(S, family = "alphaPhi")
