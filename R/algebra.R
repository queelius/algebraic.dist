join <- function(x,y) make_joint(list(unlist(x),unlist(y)))

make_joint <- function(x)
{
    if (!is.list(x)) x <- list(x)
    structure(x,class=list("joint_dist","dist"))
}

marginal.joint_dist <- function(x,indices)
{
    ifelse(length(indices) == 1,
           x[[indices[1]]],
           make_joint(x[indices]))
}

conditional.dist <- function(x,cond,n,...)
{
    make_conditional_dist(x,cond,...)
}


pdf.joint_dist <- function(x,...)
{

}



make_conditional_dist <- function(x,cond,n,...)
{
    structure(list(
        x=x,
        cond=cond,
        n=n),
        class=unique(c("conditional_dist",class(x))))
}


sampler.conditional_dist <- function(x,...)
{
    samp <- sampler(x$x)
    function(n=1)
    {
        data <- samp(x$n)
        # use x$cond so that fileter only allows cond(data[i]) == T into
        # sample
    }
}




print.dist <- function(x,...)
{
    cat("distribution\n")
    print.default(unclass(x),...)
}
