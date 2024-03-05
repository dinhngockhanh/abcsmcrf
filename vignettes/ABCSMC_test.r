library(SimBIID)

# parameters_truth <- data.frame(
#     lambda = 10,
#     mu = 2
# )
# statistics_target <- model(parameters = parameters_truth)[-c(1:ncol(parameters_truth))]
# #   Prior distributions for each parameter
# lambda <- runif(10000, 0, 15)
# mu <- runif(10000, 0, lambda)

# priors <- data.frame(
#     parnames = c("lambda","mu"),
#     dist = c("unif","unif"),
#     stringsAsFactors = FALSE
# )
# priors$p1 <- c(0,0)
# priors$p2 <- c(15,lambda)



# ==================================AFS
theta <- runif(1, 1, 20)
parameters_truth <- data.frame(
    theta = theta
)
statistics_target <- model(parameters = parameters_truth)[-c(1:ncol(parameters_truth))]
data <- c(K = statistics_target$K)
# function_ABC
func_ABC <- function(pars, data, tols, u) {
    K <- AFS_model(theta = pars, beta = u["beta"], model_type = u["model_type"], n = u["n"])[-c(1:length(pars))]
    if (all(abs(K - data) <= tols)) {
        return(K)
    } else {
        return(NA)
    }
}
u <- c(beta = 0, model_type = 1, n = 100)
# =====priors
priors <- data.frame(
    parnames = c("theta"),
    dist = c("unif"),
    stringsAsFactors = FALSE
)
priors$p1 <- c(1)
priors$p2 <- c(20)

tols <- c(K = 20)

#   Run ABC-SMC
post <- ABCSMC(
    x = data,
    priors = priors,
    func = func_ABC,
    u = u,
    tols = tols,
    ptol = 0.8,
    ngen = 7,
    npart = 5000,
    parallel = TRUE
)

filename <- "post.jpeg"
jpeg(filename, width = 2000, height = 1000)
plot(post)
dev.off()

filename <- "post_output.jpeg"
jpeg(filename, width = 2000, height = 1000)
plot(post, "output")
dev.off()

