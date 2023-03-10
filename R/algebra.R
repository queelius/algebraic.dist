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

conditional.dist <- function(x,cond,n,...)
{
  make_conditional_dist(x,cond,...)
}



print.dist <- function(x,...)
{
  cat("distribution\n")
  print.default(unclass(x),...)
}







#' \code{x1} and \code{x2} are \code{dist} (distribution) objects.
#' If \code{x1} and \code{x2} are both \code{exp_dist} objects, then
#' the result is an \code{exp_dist} object. Otherwise, the result is something
#' else which depends on how \code{min} and \code{dispatch_dist} are defined for
#' objects of type .
#' 
#' @param x1 \code{dist} object
#' @param x2 \code{dist} object
#' @export
minimum <- function(x1,x2)
{
  if (is_exp_dist(x1) && is_exp_dist(x2))
    exp_dist(rate(x1)+rate(x2))
  else
    dispatch(x1,x2,min)
}

%+%.dist <- function(x1,x2)
{
  if (is_exp_dist(x1) && is_exp_dist(x2))
    exp_dist(rate(x1)+rate(x2))
  else
    dispatch(x1,x2,add)
}
{
  if (is_normal_dist(x1) && is_normal_dist(x2))
    normal_dist(mean(x1)+mean(x2),
                var(x1)+var(x2))
}

sub <- function(x1,x2)
{
  if (is_normal_dist(x1) && is_normal_dist(x2))
    normal_dist(mean(x1)-mean(x2),
                var(x1)+var(x2))
}
