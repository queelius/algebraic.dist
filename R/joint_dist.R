joint_dist <- function(x,y) joint_dist(list(unlist(x),unlist(y)))
joint_dist <- function(x,y,z) joint_dist(list(unlist(x),unlist(y),unlist(z)))
joint_dist <- function(x)
{
  if (!is.list(x)) x <- list(x)
  ifelse(length(x) == 1, x, structure(x,class=unique(c("joint_dist","dist"))))
}

marginal.joint_dist <- function(x,indices)
{
  n <- length(indices)
  stopifnot(n >= 1 && n <= length(x))
  ifelse(n==1,x[[indices[1]]],joint_dist(x[indices]))
}

sampler.joint_dist <- function(x,...)
{
  function(n=1)
  {
    data <- NULL
    for (d in x)
      data <- cbind(data,sampler(d)(n,...))
    data
  }
}

dim.joint_dist <- function(x)
{
  n <- 0L
  for (d in x)
    n <- n + dim(d)
  n
}

pdf.joint_dist <- function(x,logp=F,...)
{
  if (logp)
  {
    function(v)
    {
      j <- 1L
      p <- 0
      for (d in x)
      {
        n <- dim(d)
        p <- p + pdf(d,logp=T,...)(v[j:(j+n)])
        j <- j + n
      }
      p
    }
  }
  else
  {
    function(v)
    {
      j <- 1L
      p <- 1
      for (d in x)
      {
        n <- dim(d)
        p <- p * pdf(d,logp=F,...)(v[j:(j+n)])
        j <- j + n
      }
      p
    }
  }
}
