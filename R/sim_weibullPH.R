#' Function to simulate Weibull with proportional hazard parametrisation
#' @param lp `<numeric()>` a vector a linear predictor i.e. sum of log hazard-ratio
#' @param shape `<integer(1)>`Weibull shape parameter, positive
#' @param scale `<integer(1)>` Weibull scale parameter,  positive
#' @param t_max `<numeric(1)>` maximium follow-up time
#' @param t_entry `<numeric()>`vector of entry time
#' @param lp_cure `<numeric()>` vector of log-odds of being cured/immortal. If a scalar, everyone shares the same immortal rate
#' @param lp_cen `<numeric()>` vector of log-hazard of dropping out. If a scalar, everyone shares the same immortal rate
#' @return a `data.frame` with these columns: tt_ev, ev, lp
sim_weibullPH <-
  function(lp,
           shape, scale,
           t_max, t_entry = rep(0, length(lp)),
           lp_cure = NULL, lp_cen = NULL){
    # Make the original dataset
    backend <- Vectorize(sim_weibullPH_one, c('lp', 't_entry'))
    tt_ev <- backend(lp, shape, scale, t_max, t_entry)

    N <- length(lp)

    # Make censoring   ----
    if (length(lp_cen)){
      t_cen <- if (length(lp_cen)==1) {
        rexp(N, p=exp(lp_cen))
      } else {
        sapply(exp(lp_cen), rexp, n=1)
      }

      t_max <- pmin(t_cen, t_max)
    }
    ev <- tt_ev <= t_max
    tt_ev <- pmin(tt_ev, t_max)

    # Make immortal list ----
    if (length(lp_cure)){
      cure <- if (length(lp_cure)==1) {
        rbinom(N, size=1, p=plogis(lp_cure))
      } else {
        sapply(plogis(lp_cure), rbinom, n=1, size=1)
      }
      ev[cure==1] <- 0
      tt_ev[cure==1] <- t_max
    }


    data.frame(
      tt_ev = tt_ev,
      ev = ev,
      lp = lp
    )
  }

sim_weibullPH_one <- function(lp, shape, scale, t_max, t_entry = 0){
  if (t_entry == 0){
    scale <- scale * exp(-lp/shape)
    return(rweibull(1, shape=shape, scale=scale))
  }
  u <- runif(1)
  inv <- function(t, lp, shape, scale, u, t_entry){
    # if (t_entry == 0) return(surv(t, lp, shape, scale) - u)
    surv(t, lp, shape, scale) - u*surv(t_entry, lp, shape, scale)
  }

  # Blhaz <- blhaz(shape, scale)
  tt_ev <-
    tryCatch({
      uniroot(
        inv,
        interval = c(0, t_max),
        lp=lp[[1]], shape=shape, scale=scale, u=u,
        tol=1e-9,
        maxiter=1e4,
        t_entry=t_entry[[1]]
      )$root
    }, error = function(e){
      if (grepl('values at end points not of opposite sign', e$message, fixed=TRUE))
        Inf
      else (stop(e$message))
    })
  tt_ev
}

#' Baseline hazard
# blhaz <- function(shape, scale){
#   function(t){
#     shape <- abs(shape)
#     scale <- abs(scale)
#     shape * (1/sqrt(scale)) * t^(shape - 1)
#   }
# }

#' Real hazard at t
# haz <- function(t, blhaz, lp){
#   # blhaz(t) * exp(lp)
# }

#' Survival at t = exp(-H) = exp(int_t(h(t)))
surv <- function(t, lp, shape, scale){
  if (t==0) return(1)
  H <- scale^(-shape) * t^shape * exp(lp)
  # H <- tryCatch(

  #   integrate(haz, lower = 0, upper = t,  blhaz=blhaz, lp=lp)$value,
  #   abs.tol = 1e-9, rel.tol=1e-9, subdivisions = 1000L,
  #   error = \(e){
  #     if (e$message == 'non-finite function value') return(Inf)
  #     stop(e$message)
  #   }
  # )

  S <- exp(-H)
  S
}
