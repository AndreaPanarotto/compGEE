softmax <- function(eta)
{
    pi <- exp(eta - max(eta)) + 1e-10 # correction
    return(pi / sum(pi))
}