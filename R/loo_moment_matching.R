#' Moment matching for efficient approximate leave-one-out cross-validation (LOO)
#'
#' Moment matching algorithm for updating a loo object when Pareto k estimates are high.
#'
#'
#' @export mmloo_manual
#'
#'
#' @param x A fitted model object.
#' @param loo A loo object that is modified.
#' @param post_draws A function the takes \code{x} as the first argument and
#'   returns a matrix of posterior draws of the model parameters.
#' @param log_lik A function that takes \code{x} and \code{i} and returns a
#'   vector of log-likeliood draws of the \code{i}th observation
#'   based on the model \code{x}.
#' @param unconstrain_pars A function that takes arguments \code{x}, and
#'   \code{pars} and returns posterior draws on the unconstrained space based on
#'   the posterior draws on the constrained space passed via \code{pars}.
#' @param log_prob_upars A function that takes arguments \code{x} and
#'   \code{upars} and returns a matrix of log-posterior density values of the
#'   unconstrained posterior draws passed via \code{upars}.
#' @param log_lik_upars A function that takes arguments \code{x}, \code{upars},
#'   and \code{i} and returns a vector of log-likeliood draws of the \code{i}th
#'   observation based on the unconstrained posterior draws passed via
#'   \code{upars}.
#' @param max_iters Maximum number of moment matching iterations. Usually this does not
#' need to be modified. If the maximum number of iterations is reached, there will be a warning, and
#' increasing \code{max_iters} may improve accuracy.
#' @param k_thres Threshold value for Pareto k values above which the moment
#'   matching algorithm is used. The default value is 0.5.
#' @param split Logical; Indicate whether to do the split transformation or not
#'   at the end of moment matching for each LOO fold.
#' @param cov Logical; Indicate whether to match the covariance matrix of the samples or not.
#'   If \code{FALSE}, only the mean and marginal variances are matched.
#' @template cores
#' @param ... Further arguments passed to the custom functions documented above.
#'
#' @return An updated \code{loo} object.
#'
#'
#' @seealso [loo()], [psis()], [loo_compare()]
#' @template moment-matching-references
#'
#'
#'
#'
#'
mmloo_manual <- function(x, loo, post_draws, log_lik,
                         unconstrain_pars, log_prob_upars,
                         log_lik_upars, max_iters = 30L,
                         k_thres = 0.5, split = TRUE,
                         cov = TRUE, cores = getOption("mc.cores", 1),
                         ...) {

  # input checks
  checkmate::assertClass(loo,classes = "loo")
  checkmate::assertFunction(post_draws)
  checkmate::assertFunction(log_lik)
  checkmate::assertFunction(unconstrain_pars)
  checkmate::assertFunction(log_prob_upars)
  checkmate::assertFunction(log_lik_upars)
  checkmate::assertNumber(max_iters)
  checkmate::assertNumber(k_thres)
  checkmate::assertLogical(split)
  checkmate::assertLogical(cov)
  checkmate::assertNumber(cores)


  S <- dim(loo)[1]
  N <- dim(loo)[2]
  pars <- post_draws(x, ...)
  # transform the model parameters to unconstrained space
  upars <- unconstrain_pars(x, pars = pars, ...)
  # number of parameters in the **parameters** block only
  npars <- dim(upars)[2]
  # if more parameters than samples, do not do Cholesky transformation
  cov <- cov && S >= npars
  # compute log-probabilities of the original parameter values
  orig_log_prob <- log_prob_upars(x, upars = upars, ...)

  # loop over all observations whose Pareto k is high
  ks <- loo$diagnostics$pareto_k
  kfs <- rep(0,N)
  r_effs <- loo$diagnostics$n_eff
  I <- which(ks > k_thres)
  for (i in I) {
    message("Moment matching observation ", i)
    # initialize values for this LOO-fold
    uparsi <- upars
    ki <- ks[i]
    kfi <- 0
    log_liki <- log_lik(x, i, ...)
    r_effi <- r_effs[i]

    # compute log-weights per draw
    lwi <- -log_liki
    lwi <- lwi - max(lwi)
    psis_i <- suppressWarnings(loo::psis(lwi, r_eff = r_effi, cores = cores))
    lwi <- as.vector(weights(psis_i))

    # initialize objects that keep track of the total transformation
    total_shift <- rep(0, npars)
    total_scaling <- rep(1, npars)
    total_mapping <- diag(npars)

    # try several transformations one by one
    # if one does not work, do not apply it and try another one
    # to accept the transformation, Pareto k needs to improve
    # when transformation succeeds, start again from the first one
    iterind <- 1
    while (iterind <= max_iters && ki > k_thres) {

      if (iterind == max_iters) {
        throw_moment_match_max_iters_warning()
      }

      # 1. match means
      trans <- shift(x, uparsi, lwi)
      # gather updated quantities
      quantities_i <- update_quantities_i(x, trans$upars,  i = i,
                                          orig_log_prob = orig_log_prob,
                                          log_prob_upars = log_prob_upars,
                                          log_lik_upars = log_lik_upars,
                                          r_effi = r_effi,
                                          cores = cores,
                                          ...)
      if (quantities_i$ki < ki) {
        uparsi <- trans$upars
        total_shift <- total_shift + trans$shift

        lwi <- quantities_i$lwi
        ki <- quantities_i$ki
        kfi <- quantities_i$kfi
        log_liki <- quantities_i$log_liki
        iterind <- iterind + 1
        next
      }

      # 2. match means and marginal variances
      trans <- shift_and_scale(x, uparsi, lwi)
      # gather updated quantities
      quantities_i <- update_quantities_i(x, trans$upars,  i = i,
                                          orig_log_prob = orig_log_prob,
                                          log_prob_upars = log_prob_upars,
                                          log_lik_upars = log_lik_upars,
                                          r_effi = r_effi,
                                          cores = cores,
                                          ...)
      if (quantities_i$ki < ki) {
        uparsi <- trans$upars
        total_shift <- total_shift + trans$shift
        total_scaling <- total_scaling * trans$scaling

        lwi <- quantities_i$lwi
        ki <- quantities_i$ki
        kfi <- quantities_i$kfi
        log_liki <- quantities_i$log_liki
        iterind <- iterind + 1
        next
      }

      # 3. match means and covariances
      if (cov) {
        trans <- shift_and_cov(x, uparsi, lwi)
        # gather updated quantities
        quantities_i <- update_quantities_i(x, trans$upars,  i = i,
                                            orig_log_prob = orig_log_prob,
                                            log_prob_upars = log_prob_upars,
                                            log_lik_upars = log_lik_upars,
                                            r_effi = r_effi,
                                            cores = cores,
                                            ...)

        if (quantities_i$ki < ki) {
          uparsi <- trans$upars
          total_shift <- total_shift + trans$shift
          total_mapping <- trans$mapping %*% total_mapping

          lwi <- quantities_i$lwi
          ki <- quantities_i$ki
          kfi <- quantities_i$kfi
          log_liki <- quantities_i$log_liki
          iterind <- iterind + 1
          next
        }
      }
      # none of the transformations improved khat
      # so there is no need to try further
      break
    }

    # transformations are now done
    # if we don't do split transform, or
    # if no transformations were successful
    # stop and collect values
    if (!split || (iterind == 1)) {
      elpd_loo_i <- matrixStats::logSumExp(log_liki + lwi)
    } else {
      # compute split transformation
      split_mm <- split_mm(
        x, upars, cov, total_shift, total_scaling, total_mapping, i,
        log_prob_upars = log_prob_upars, log_lik_upars = log_lik_upars,
        cores = cores, r_effi = r_effi
      )
      elpd_loo_i <- matrixStats::logSumExp(
        split_mm$log_liki + split_mm$lwi
      )
      log_liki <- split_mm$log_liki
      lwi <- split_mm$lwi
    }


    # pointwise estimates
    # p_loo: use original p_loo, add original elpd, subtract new elpd
    loo$pointwise[i, 3] <- loo$pointwise[i, 3] +
      loo$pointwise[i, 1] - elpd_loo_i
    # elpd_loo
    loo$pointwise[i, 1] <- elpd_loo_i
    # mcse_elpd_loo
    E_epd_i <- exp(elpd_loo_i)
    var_epd_i <- sum((exp(lwi))^2 * (exp(log_liki) - E_epd_i)^2)
    z <- rnorm(1000, mean = E_epd_i, sd = sqrt(var_epd_i))
    loo$pointwise[i, 2] <- sqrt(var(log(z[z > 0])) / r_effi)
    # looic
    loo$pointwise[i, 4] <- -2 * elpd_loo_i

    # diagnostics
    loo$diagnostics$pareto_k[i] <- ki
    loo$diagnostics$n_eff[i] <- 1.0 / sum(exp(2 * lwi)) * r_effi
    kfs[i] <- kfi

    # update psis object
    if (!is.null(loo$psis_object)) {
      loo$psis_object$log_weights[, i] <- lwi
    }

    if (!split) {
      throw_large_kf_warning(kfs)
    }

  }
  # combined estimates
  loo$estimates[1, 1] <- sum(loo$pointwise[, 1])
  loo$estimates[2, 1] <- sum(loo$pointwise[, 3])
  loo$estimates[3, 1] <- sum(loo$pointwise[, 4])
  loo$estimates[1, 2] <- sqrt(N * var(loo$pointwise[, 1]))
  loo$estimates[2, 2] <- sqrt(N * var(loo$pointwise[, 3]))
  loo$estimates[3, 2] <- sqrt(N * var(loo$pointwise[, 4]))
  # loo$estimates <- loo:::table_of_estimates(loo$pointwise)

  # update psis object
  if (!is.null(loo$psis_object)) {
    loo$psis_object$diagnostics <- loo$diagnostics
  }


  # these will be deprecated at some point
  loo$elpd_loo <- loo$estimates[1, 1]
  loo$p_loo <- loo$estimates[2, 1]
  loo$looic <- loo$estimates[3, 1]
  loo$se_elpd_loo <- loo$estimates[1, 2]
  loo$se_p_loo <- loo$estimates[2, 2]
  loo$se_looic <- loo$estimates[3, 2]

  # Warn if some Pareto ks are still high
  psislw_warnings(loo$diagnostics$pareto_k)

  loo
}










#' Update the importance weights, Pareto diagnostic and log-likelihood
#' for observation \code{i} based on model \code{x}.
#'
#' @noRd
#' @param x A fitted model object.
#' @param upars A matrix representing a sample of vector-valued parameters in the unconstrained space.
#' @param i observation number.
#' @param orig_log_prob log probability densities of the original draws from the model \code{x}.
#' @param log_prob_upars A function that takes arguments \code{x} and
#'   \code{upars} and returns a matrix of log-posterior density values of the
#'   unconstrained posterior draws passed via \code{upars}.
#' @param log_lik_upars A function that takes arguments \code{x}, \code{upars},
#'   and \code{i} and returns a vector of log-likeliood draws of the \code{i}th
#'   observation based on the unconstrained posterior draws passed via
#'   \code{upars}.
#' @param r_effi MCMC effective sample size divided by the total sample size for 1/exp(log_ratios) for observation i.
#' @template cores
#' @return List with the updated importance weights, Pareto diagnostics and log-likelihood values.
#'
update_quantities_i <- function(x, upars, i, orig_log_prob,
                                log_prob_upars, log_lik_upars,
                                r_effi, cores, ...) {
  log_prob_new <- log_prob_upars(x, upars = upars, ...)
  log_liki_new <- log_lik_upars(x, upars = upars, i = i, ...)
  # compute new raw log importance weights, unnormalized
  lwi_new <- -log_liki_new + log_prob_new - orig_log_prob
  lwi_new <- lwi_new - max(lwi_new)
  psis_new <- suppressWarnings(loo::psis(lwi_new, r_eff = r_effi, cores = cores))
  ki_new <- psis_new$diagnostics$pareto_k

  lwfi_new <- log_prob_new - orig_log_prob
  lwfi_new <- lwfi_new - max(lwfi_new)
  psisf_new <- suppressWarnings(loo::psis(lwfi_new, r_eff = r_effi, cores = cores))
  kfi_new <- psisf_new$diagnostics$pareto_k
  # pareto smoothed weights
  lwi_new <- as.vector(weights(psis_new))
  # gather results
  list(
    lwi = lwi_new,
    ki = ki_new,
    kfi = kfi_new,
    log_liki = log_liki_new
  )
}



#' Shift a matrix of parameters to their weighted mean.
#' Also calls update_quantities_i which updates the importance weights based on the supplied model object.
#'
#' @noRd
#' @param x A fitted model object.
#' @param upars A matrix representing a sample of vector-valued parameters in the unconstrained space
#' @param lwi A vector representing the log-weight of each parameter
#' @return List with the shift that was performed, and the new parameter matrix.
#'
shift <- function(x, upars, lwi) {
  # compute moments using log weights
  mean_original <- colMeans(upars)
  mean_weighted <- colSums(exp(lwi) * upars)
  shift <- mean_weighted - mean_original
  # transform posterior draws
  upars_new <- sweep(upars, 2, shift, "+")
  list(
    upars = upars_new,
    shift = shift
  )
}




#' Shift a matrix of parameters to their weighted mean and scale the marginal variances
#' to match the weighted marginal variances.
#' Also calls update_quantities_i which updates the importance weights based on the supplied model object.
#'
#' @noRd
#' @param x A fitted model object.
#' @param upars A matrix representing a sample of vector-valued parameters in the unconstrained space
#' @param lwi A vector representing the log-weight of each parameter
#' @return List with the shift and scaling that were performed, and the new parameter matrix.
#'
#'
shift_and_scale <- function(x, upars, lwi) {
  # compute moments using log weights
  S <- dim(upars)[1]
  mean_original <- colMeans(upars)
  mean_weighted <- colSums(exp(lwi) * upars)
  shift <- mean_weighted - mean_original
  mii <- exp(lwi)* upars^2
  mii <- colSums(mii) - mean_weighted^2
  mii <- mii*S/(S-1)
  scaling <- sqrt(mii / matrixStats::colVars(upars))
  # transform posterior draws
  upars_new <- sweep(upars, 2, mean_original, "-")
  upars_new <- sweep(upars_new, 2, scaling, "*")
  upars_new <- sweep(upars_new, 2, mean_weighted, "+")

  list(
    upars = upars_new,
    shift = shift,
    scaling = scaling
  )
}

#' Shift a matrix of parameters to their weighted mean and scale the covariance
#' to match the weighted covariance.
#' Also calls update_quantities_i which updates the importance weights based on the supplied model object.
#'
#' @noRd
#' @param x A fitted model object.
#' @param upars A matrix representing a sample of vector-valued parameters in the unconstrained space
#' @param lwi A vector representing the log-weight of each parameter
#' @return List with the shift and mapping that were performed, and the new parameter matrix.
#'
shift_and_cov <- function(x, upars, lwi, ...) {
  # compute moments using log weights
  mean_original <- colMeans(upars)
  mean_weighted <- colSums(exp(lwi) * upars)
  shift <- mean_weighted - mean_original
  covv <- stats::cov(upars)
  wcovv <- stats::cov.wt(upars, wt = exp(lwi))$cov
  chol1 <- chol(wcovv)
  chol2 <- chol(covv)
  mapping <- t(chol1) %*% solve(t(chol2))
  # transform posterior draws
  upars_new <- sweep(upars, 2, mean_original, "-")
  upars_new <- tcrossprod(upars_new, mapping)
  upars_new <- sweep(upars_new, 2, mean_weighted, "+")
  colnames(upars_new) <- colnames(upars)

  list(
    upars = upars_new,
    shift = shift,
    mapping = mapping
  )
}


#' Warning message if max_iters is reached
#' @noRd
throw_moment_match_max_iters_warning <- function() {
  warning(
    "The maximum number of moment matching iterations ('max_iters' argument) was reached.\n",
    "Increasing the value may improve accuracy.",
    call. = FALSE
  )
}

#' Warning message if not using split transformation and accuracy is compromised
#' @noRd
throw_large_kf_warning <- function(kf) {
  if (any(kf > 0.5)) {
    warning(
      "The accuracy of self-normalized importance sampling may be bad.\n",
      "Setting the argument 'split' to 'TRUE' will likely improve accuracy.",
      call. = FALSE
    )
  }

}

#' warnings about pareto k values ------------------------------------------
#' @noRd
psislw_warnings <- function(k) {
  if (any(k > 0.7)) {
    .warn(
      "Some Pareto k diagnostic values are too high. ",
      .k_help()
    )
  } else if (any(k > 0.5)) {
    .warn(
      "Some Pareto k diagnostic values are slightly high. ",
      .k_help()
    )
  }
}
