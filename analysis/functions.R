# functions.R
# DATA PROCESSING --------------------------------------------------------------
# extract covariates at a set of a points from raster file
get_cov <- function(file, points, assign_na = T) {
  cov <- rast(file)
  # get covariate value at points
  cov_pts <- extract(cov, points)
  
  if (assign_na) {
    coarse_cov <- aggregate(cov, fact = 10)
    nearest_points <- 
      st_join(points[which(is.na(cov_pts[, 2])), ], 
              st_as_sf(as.points(coarse_cov)), 
              join = st_nearest_feature) |>
      st_set_geometry(NULL)
    
    cov_pts[which(is.na(cov_pts[, 2])), ] <- nearest_points[, ncol(points)]
  }
  
  
  return(cov_pts[, 2])
}
# UNIT LEVEL MODELS ------------------------------------------------------------
fitLonobinomialBYM2 <- function(formula,
                                domain, 
                                design,
                                X.pop = NULL,
                                adj.mat = NULL,
                                domain.size = NULL,
                                pc.u = 1,
                                pc.alpha = 0.01,
                                pc.u.phi = 0.5,
                                pc.alpha.phi = 2/3,
                                pc.prior.clust = 
                                  list(prec = list(prior = "pc.prec",
                                                   param = c(1, 0.05))),
                                overdisp.mean = 0, overdisp.prec = 0.4,
                                level = .95, n.sample = 250,
                                return.samples = F,
                                X.pop.weights = NULL) {
  
  # SETUP ----------------------------------------------------------------------
  
  # set up return value
  out <- list()
  attr(out, "inla.fitted") <- c()
  
  # get domain variable
  domain.var <- all.vars(domain)
  if (length(domain.var) != 1) {
    stop("Domain labels must be contained in a single column. Spatio-temporal not supported yet.")
  } 
  
  # set up formulas
  resp.frm <- as.formula(paste0("~", all.vars(formula)[1]))
  cov.frm <- update(formula, NULL ~ .)
  
  # DIRECT ESTIMATES -----------------------------------------------------------
  # compute direct estimates (Hajek estimates if domain size unknown)
  if (!is.null(domain.size)) {
    direct.est <- svyby(resp.frm, domain, design = design, svytotal, na.rm = T)
    direct.est$domain.size <- 
      domain.size$size[match(direct.est[, 1], domain.size[[domain.var]])]
    direct.est[, 2] = direct.est[, 2] / direct.est$domain.size
    direct.est[, 3] = (direct.est[, 3] / direct.est$domain.size) ^ 2
    direct.est <- direct.est[, 1:3]
  } else {
    direct.est <- svyby(resp.frm, domain, design = design, svymean, na.rm = T)
    direct.est[, 3] <- direct.est[, 3] ^ 2
  }
  rownames(direct.est) <- NULL
  colnames(direct.est) <- c("domain", "mean", "var")
  direct.est <- 
    data.frame(domain = direct.est$domain,
               mean = direct.est$mean,
               median = direct.est$mean,
               var = direct.est$var,
               lower = direct.est$mean + qnorm((1-level)/2) * sqrt(direct.est$var),
               upper = direct.est$mean + qnorm(1 - (1-level)/2) * sqrt(direct.est$var),
               method = paste0("Direct"))
  
  
  if (!is.null(adj.mat) & !setequal(X.pop[[domain.var]], rownames(adj.mat))) {
    stop("Domains in X.pop do not match domains in adj.mat.")
  }
  if (any(is.na(match(X.pop[[domain.var]], direct.est$domain)))) {
    warning(cat("There are domains in X.pop not in design/direct estimates.",
                "\nGenerating estimates for all domains in X.pop."))
  }
  # if no adjacency matrix matches, take domain names from X.pop
  domain.table <- data.frame(domain = unique(as.character(X.pop[[domain.var]])))
  # if adjacency matrix provided, take domain names from row names
  if (!is.null(adj.mat)) {
    domain.table <- data.frame(domain = rownames(adj.mat))
  }
  
  direct.est <- 
    merge(direct.est, data.frame(domain = domain.table$domain), 
          by = "domain", all.y = T)
  direct.est$method = "Direct"
  
  out$direct.est <- direct.est
  attr(out, "domain.names") <- sort(direct.est$domain)
  attr(out, "method.names") <- c("direct.est")
  
  # UNIT LEVEL MODEL -----------------------------------------------------------

  mf <- model.frame(formula, design$variables)
  resp <- model.response(mf, "numeric")
  
  pop.dat <- data.frame(
    domain = X.pop[[domain.var]]
  )
  mm.pop <- model.matrix(cov.frm, X.pop)
  mod.dat <- data.frame(
    resp = as.vector(resp),
    domain = design$variables[[domain.var]],
    cluster = design$cluster[, 1]
  )
  mod.dat <- cbind(mod.dat, model.matrix(cov.frm, design$variables))

  # domain labels as indexes for use in INLA
  domain.table$domain.struct <- seq_len(nrow(domain.table))
  mod.dat$domain.struct <- match(mod.dat$domain, domain.table$domain)
  pop.dat$domain.struct <- match(pop.dat$domain, domain.table$domain)
  
  # model formula
  ftxt <- paste("resp ~ 1")
  if (length(all.vars(cov.frm)) > 0) {
    ftxt <- 
      paste(ftxt, as.character(cov.frm)[-1], sep = " + ")
  }
  
  if (is.null(adj.mat)) {
    # set priors
    hyperpc.iid.int <- list(
      prec = list(prior = "pc.prec",
                  param = c(pc.u , pc.alpha))
    )
    warning("No spatial information provided, using iid domain effects")
    model.method <- "lonobinomial.iid.model"
    ftxt <- paste0(ftxt, " + f(domain.struct, model = 'iid', hyper = hyperpc.iid.int)")
  } else {
    hyperpc.bym.int <- list(
      prec = list(prior = "pc.prec", param = c(pc.u , pc.alpha)),  
      phi = list(prior = 'pc', param = c(pc.u.phi , pc.alpha.phi))
    )
    model.method <- "lonobinomial.bym2.model"
    ftxt <- paste0(ftxt, " + f(domain.struct, model = 'bym2', graph=adj.mat, hyper = hyperpc.bym.int)")
  }
  
  ftxt <- paste0(ftxt, "+ f(cluster, model = 'iid', hyper = pc.prior.clust)")
  
  # fit model
  mod.frm <- as.formula(ftxt)
  fit <- INLA::inla(mod.frm, family = "binomial", Ntrials = mod.dat$n, 
                    data = mod.dat,
                    control.compute = list(dic = T, mlik = T,
                                           cpo = T, config = TRUE),
                    control.predictor = list(compute = TRUE),
                    lincomb = NULL, quantiles = c((1-level)/2, 0.5, 1-(1-level)/2))
  
  # generate draws
  samp.all <- INLA::inla.posterior.sample(n = n.sample, result = fit, intern = TRUE)

  # identify indices of fixed effects and random effects
  fe.idx <- grep(colnames(mm.pop)[1], rownames(samp.all[[1]]$latent))
  fe.idx <- fe.idx:(fe.idx + ncol(mm.pop) - 1)
  re.idx <- grep("domain.struct", x = rownames(samp.all[[1]]$latent))
  # aggregate sample predictions
  summary.sample <- function(x) {
    clust.sd <- sqrt(1 / exp(x$hyperpar[3]))
    pop.unit.ests <-
      x$latent[re.idx][pop.dat$domain.struct] + mm.pop %*% x$latent[fe.idx] 
    pop.unit.ests <- 
      pop.unit.ests + 
      rnorm(length(pop.unit.ests), sd = clust.sd)
    pop.unit.ests <- expit(pop.unit.ests)
    if (!is.null(X.pop.weights)) {
      area.ests <- 
        aggregate(pop.unit.ests * X.pop.weights,  list(domain = pop.dat$domain.struct), sum)
    } else {
      area.ests <- 
        aggregate(pop.unit.ests, list(domain = pop.dat$domain.struct), mean)
    }
    
    return(area.ests[match(1:nrow(domain.table), area.ests[, 1]), 2])
  }
  est.mat <- do.call(cbind, lapply(samp.all, summary.sample))
  out[[paste0(model.method, ".fit")]] <- fit
  out[[paste0(model.method, ".est")]]  <-
    data.frame(domain = domain.table$domain,
               mean = rowMeans(est.mat),
               median = apply(est.mat, 1,
                              function(x) median(x, na.rm = T)),
               var = apply(est.mat, 1, var),
               lower = apply(est.mat, 1,
                             function(x) quantile(x, (1-level)/2, na.rm = T)),
               upper = apply(est.mat, 1,
                             function(x) quantile(x, 1-(1-level)/2, na.rm = T)),
               method = paste0("Lonobinomial: ", 
                               ifelse(is.null(adj.mat), "IID", "BYM2")))
  if (return.samples) {
    out[[paste0(model.method, ".sample")]]  <- est.mat
  }
  out[[paste0(model.method, ".est")]] <- 
    out[[paste0(model.method, ".est")]][match(out$direct.est$domain, 
                                              out[[paste0(model.method, ".est")]]$domain),]
  attr(out, "method.names") <- c(attr(out, "method.names"), paste0(model.method, ".est"))
  attr(out, "inla.fitted") <- c(attr(out, "inla.fitted"), model.method)
  
  # finish building return value
  out$call <- match.call()
  class(out) <- c("svysae")
  return(out)
  
}

fitBetabinomialBYM2 <- function(formula,
                                domain, 
                                design,
                                X.pop = NULL,
                                adj.mat = NULL,
                                domain.size = NULL,
                                pc.u = 1,
                                pc.alpha = 0.01,
                                pc.u.phi = 0.5,
                                pc.alpha.phi = 2/3,
                                overdisp.mean = 0, overdisp.prec = 0.4,
                                level = .95, n.sample = 250,
                                return.samples = F,
                                X.pop.weights = NULL) {
  
  # SETUP ----------------------------------------------------------------------
  
  # set up return value
  out <- list()
  attr(out, "inla.fitted") <- c()
  
  # get domain variable
  domain.var <- all.vars(domain)
  if (length(domain.var) != 1) {
    stop("Domain labels must be contained in a single column. Spatio-temporal not supported yet.")
  } 
  
  # set up formulas
  resp.frm <- as.formula(paste0("~", all.vars(formula)[1]))
  cov.frm <- update(formula, NULL ~ .)
  
  # DIRECT ESTIMATES -----------------------------------------------------------
  # compute direct estimates (Hajek estimates if domain size unknown)
  if (!is.null(domain.size)) {
    direct.est <- svyby(resp.frm, domain, design = design, svytotal, na.rm = T)
    direct.est$domain.size <- 
      domain.size$size[match(direct.est[, 1], domain.size[[domain.var]])]
    direct.est[, 2] = direct.est[, 2] / direct.est$domain.size
    direct.est[, 3] = (direct.est[, 3] / direct.est$domain.size) ^ 2
    direct.est <- direct.est[, 1:3]
  } else {
    direct.est <- svyby(resp.frm, domain, design = design, svymean, na.rm = T)
    direct.est[, 3] <- direct.est[, 3] ^ 2
  }
  rownames(direct.est) <- NULL
  colnames(direct.est) <- c("domain", "mean", "var")
  direct.est <- 
    data.frame(domain = direct.est$domain,
               mean = direct.est$mean,
               median = direct.est$mean,
               var = direct.est$var,
               lower = direct.est$mean + qnorm((1-level)/2) * sqrt(direct.est$var),
               upper = direct.est$mean + qnorm(1 - (1-level)/2) * sqrt(direct.est$var),
               method = paste0("Direct"))
  
  
  if (!is.null(adj.mat) & !setequal(X.pop[[domain.var]], rownames(adj.mat))) {
    stop("Domains in X.pop do not match domains in adj.mat.")
  }
  if (any(is.na(match(X.pop[[domain.var]], direct.est$domain)))) {
    warning(cat("There are domains in X.pop not in design/direct estimates.",
                "\nGenerating estimates for all domains in X.pop."))
  }
  # if no adjacency matrix matches, take domain names from X.pop
  domain.table <- data.frame(domain = unique(as.character(X.pop[[domain.var]])))
  # if adjacency matrix provided, take domain names from row names
  if (!is.null(adj.mat)) {
    domain.table <- data.frame(domain = rownames(adj.mat))
  }
  
  direct.est <- 
    merge(direct.est, data.frame(domain = domain.table$domain), 
          by = "domain", all.y = T)
  direct.est$method = "Direct"
  
  out$direct.est <- direct.est
  attr(out, "domain.names") <- sort(direct.est$domain)
  attr(out, "method.names") <- c("direct.est")
  
  # UNIT LEVEL MODEL -----------------------------------------------------------

  mf <- model.frame(formula, design$variables)
  resp <- model.response(mf, "numeric")
  
  pop.dat <- data.frame(
    domain = X.pop[[domain.var]]
  )
  mm.pop <- model.matrix(cov.frm, X.pop)
  mod.dat <- data.frame(
    resp = as.vector(resp),
    domain = design$variables[[domain.var]],
    cluster = design$cluster[, 1]
  )
  mod.dat <- cbind(mod.dat, model.matrix(cov.frm, design$variables))
  mod.dat <- mod.dat |>
    group_by(cluster) |>
    mutate(resp = sum(resp),
           n = n()) |>
    slice(which.min(resp)) |>
    ungroup() 

  
  # model formula
  ftxt <- paste("resp ~ 1")
  if (length(all.vars(cov.frm)) > 0) {
    ftxt <- 
      paste(ftxt, paste0(colnames(model.matrix(cov.frm, design$variables))[-1], collapse = "+"), sep = " + ")
  }
  
  # domain labels as indexes for use in INLA
  domain.table$domain.struct <- seq_len(nrow(domain.table))
  mod.dat$domain.struct <- match(mod.dat$domain, domain.table$domain)
  pop.dat$domain.struct <- match(pop.dat$domain, domain.table$domain)
  
  
  if (is.null(adj.mat)) {
    # set priors
    hyperpc.iid.int <- list(
      prec = list(prior = "pc.prec",
                  param = c(pc.u , pc.alpha))
    )
    warning("No spatial information provided, using iid domain effects")
    model.method <- "betabinomial.iid.model"
    ftxt <- paste0(ftxt, " + f(domain.struct, model = 'iid', hyper = hyperpc.iid.int)")
  } else {
    hyperpc.bym.int <- list(
      prec = list(prior = "pc.prec", param = c(pc.u , pc.alpha)),  
      phi = list(prior = 'pc', param = c(pc.u.phi , pc.alpha.phi))
    )
    model.method <- "betabinomial.bym2.model"
    
    ftxt <- paste0(ftxt, " + f(domain.struct, model = 'bym2', graph=adj.mat, hyper = hyperpc.bym.int, scale.model=T, constr=T, adjust.for.con.comp=T)")
  }
  control.family <- list(hyper = list(rho = list(param = c(overdisp.mean, overdisp.prec))))
  
 
  
  # fit model
  mod.frm <- as.formula(ftxt)
  fit <- INLA::inla(mod.frm, family = "betabinomial", Ntrials = mod.dat$n, 
                    data = mod.dat,
                    control.family = control.family,
                    control.compute = list(dic = T, mlik = T,
                                           cpo = T, config = TRUE),
                    control.predictor = list(compute = TRUE),
                    lincomb = NULL, quantiles = c((1-level)/2, 0.5, 1-(1-level)/2))
  
  
  # generate draws
  samp.all <- INLA::inla.posterior.sample(n = n.sample, result = fit, intern = TRUE)

  # identify indices of fixed effects and random effects
  fe.idx <- grep(colnames(mm.pop)[1], rownames(samp.all[[1]]$latent))
  fe.idx <- fe.idx:(fe.idx + ncol(mm.pop) - 1)
  re.idx <- grep("domain.struct", x = rownames(samp.all[[1]]$latent))

  # aggregate sample predictions
  summary.sample <- function(x) {
    pop.unit.ests <-
      x$latent[re.idx][pop.dat$domain.struct] + mm.pop %*% x$latent[fe.idx] 
    pop.unit.ests <- expit(pop.unit.ests)
    if (!is.null(X.pop.weights)) {
      area.ests <- 
        aggregate(pop.unit.ests * X.pop.weights,  list(domain = pop.dat$domain.struct), sum)
    } else {
      area.ests <- 
        aggregate(pop.unit.ests, list(domain = pop.dat$domain.struct), mean)
    }
    
    return(area.ests[match(1:nrow(domain.table), area.ests[, 1]), 2])
  }
  est.mat <- do.call(cbind, lapply(samp.all, summary.sample))
  out[[paste0(model.method, ".fit")]] <- fit
  out[[paste0(model.method, ".est")]]  <-
    data.frame(domain = domain.table$domain,
               mean = rowMeans(est.mat),
               median = apply(est.mat, 1,
                              function(x) median(x, na.rm = T)),
               var = apply(est.mat, 1, var),
               lower = apply(est.mat, 1,
                             function(x) quantile(x, (1-level)/2, na.rm = T)),
               upper = apply(est.mat, 1,
                             function(x) quantile(x, 1-(1-level)/2, na.rm = T)),
               method = paste0("Betabinomial: ", 
                               ifelse(is.null(adj.mat), "IID", "BYM2")))
  if (return.samples) {
    out[[paste0(model.method, ".sample")]]  <- est.mat
  }
  out[[paste0(model.method, ".est")]] <- 
    out[[paste0(model.method, ".est")]][match(out$direct.est$domain, 
                                              out[[paste0(model.method, ".est")]]$domain),]
  attr(out, "method.names") <- c(attr(out, "method.names"), paste0(model.method, ".est"))
  attr(out, "inla.fitted") <- c(attr(out, "inla.fitted"), model.method)
  
  # finish building return value
  out$call <- match.call()
  class(out) <- c("svysae")
  return(out)
  
}

#### NESTED MODEL ####

smoothUnitNested <- function(formula,
                       domain, 
                       design,
                       family = c("gaussian", "binomial")[1],
                       X.pop = NULL,
                       adj.mat = NULL,
                       domain.size = NULL,
                       pc.u = 1,
                       pc.alpha = 0.01, 
                       pc.u.phi = 0.5,
                       pc.alpha.phi = 2/3,
                       level = .95, n.sample = 250,
                       return.samples = F,
                       X.pop.weights = NULL) {
  
  if (design$has.strata) {
    warning("This model does not account for stratification yet.")
  }
  if (ncol(design$cluster) > 2) {
    stop("This function is not ready for > 2 stages of sampling")
  }
  # RESPONSE FAMILY
  if (family != "binomial" & family != "gaussian" ) {
    stop(paste0("family = ", family, " not supported"))
  }
  
  # SETUP ----------------------------------------------------------------------
  
  # set up return value
  out <- list()
  attr(out, "inla.fitted") <- c()
  
  # get domain variable
  domain.var <- all.vars(domain)
  if (length(domain.var) != 1) {
    stop("Domain labels must be contained in a single column. Spatio-temporal not supported yet.")
  } 
  
  # set up formulas
  resp.frm <- as.formula(paste0("~", all.vars(formula)[1]))
  cov.frm <- update(formula, NULL ~ .)
  
  # DIRECT ESTIMATES -----------------------------------------------------------
  # compute direct estimates (Hajek estimates if domain size unknown)
  if (!is.null(domain.size)) {
    direct.est <- svyby(resp.frm, domain, design = design, svytotal, na.rm = T)
    direct.est$domain.size <- 
      domain.size$size[match(direct.est[, 1], domain.size[[domain.var]])]
    direct.est[, 2] = direct.est[, 2] / direct.est$domain.size
    direct.est[, 3] = (direct.est[, 3] / direct.est$domain.size) ^ 2
    direct.est <- direct.est[, 1:3]
  } else {
    direct.est <- svyby(resp.frm, domain, design = design, svymean, na.rm = T)
    direct.est[, 3] <- direct.est[, 3] ^ 2
  }
  rownames(direct.est) <- NULL
  colnames(direct.est) <- c("domain", "mean", "var")
  direct.est <- 
    data.frame(domain = direct.est$domain,
               mean = direct.est$mean,
               median = direct.est$mean,
               var = direct.est$var,
               lower = direct.est$mean + qnorm((1-level)/2) * sqrt(direct.est$var),
               upper = direct.est$mean + qnorm(1 - (1-level)/2) * sqrt(direct.est$var),
               method = paste0("Direct"))
  
  
  if (!is.null(adj.mat) & !setequal(X.pop[[domain.var]], rownames(adj.mat))) {
    stop("Domains in X.pop do not match domains in adj.mat.")
  }
  if (any(is.na(match(X.pop[[domain.var]], direct.est$domain)))) {
    warning(cat("There are domains in X.pop not in design/direct estimates.",
                "\nGenerating estimates for all domains in X.pop."))
  }
  # if no adjacency matrix matches, take domain names from X.pop
  domain.table <- data.frame(domain = unique(as.character(X.pop[[domain.var]])))
  # if adjacency matrix provided, take domain names from row names
  if (!is.null(adj.mat)) {
    domain.table <- data.frame(domain = rownames(adj.mat))
  }
  
  direct.est <- 
    merge(direct.est, data.frame(domain = domain.table$domain), 
          by = "domain", all.y = T)
  direct.est$method = "Direct"
  
  out$direct.est <- direct.est
  attr(out, "domain.names") <- sort(direct.est$domain)
  attr(out, "method.names") <- c("direct.est")
  
  # UNIT LEVEL MODEL -----------------------------------------------------------
  mf <- model.frame(formula, design$variables)
  resp <- model.response(mf, "numeric")
  
  pop.dat <- data.frame(
    domain = X.pop[[domain.var]]
  )
  mm.pop <- model.matrix(cov.frm, X.pop)
  mod.dat <- data.frame(
    resp = as.vector(resp),
    domain = design$variables[[domain.var]]
  )
  mod.dat <- cbind(mod.dat, model.matrix(cov.frm, design$variables))
  
  
  
  # model formula
  ftxt <- paste("resp ~ 1")
  if (length(all.vars(cov.frm)) > 0) {
    ftxt <- 
      paste(ftxt, paste0(colnames(model.matrix(cov.frm, design$variables))[-1], collapse = "+"), sep = " + ")
  }
  
  # domain labels as indexes for use in INLA
  domain.table$domain.struct <- seq_len(nrow(domain.table))
  mod.dat$domain.struct <- match(mod.dat$domain, domain.table$domain)
  pop.dat$domain.struct <- match(pop.dat$domain, domain.table$domain)
  
  if (is.null(adj.mat)) {
    # set priors
    hyperpc.iid.int <- list(
      prec = list(prior = "pc.prec",
                  param = c(pc.u , pc.alpha))
    )
    warning("No spatial information provided, using iid domain effects")
    model.method <- "iid.model"
    ftxt <- paste0(ftxt, " + f(domain.struct, model = 'iid', hyper = hyperpc.iid.int)")
  } else {
    hyperpc.bym.int <- list(
      prec = list(prior = "pc.prec", param = c(pc.u , pc.alpha)),  
      phi = list(prior = 'pc', param = c(pc.u.phi , pc.alpha.phi))
    )
    model.method <- "bym2.model"
    ftxt <- paste0(ftxt, " + f(domain.struct, model = 'bym2', graph=adj.mat, hyper = hyperpc.bym.int, scale.model=T, constr=T, adjust.for.con.comp=T)")
  }
  # fit model
  mod.frm <- as.formula(ftxt)
  print(mod.frm)
  fit <- INLA::inla(mod.frm, family = family, data = mod.dat,
                    control.compute = list(dic = T, mlik = T,
                                           cpo = T, config = TRUE),
                    control.predictor = list(compute = TRUE),
                    lincomb = NULL, quantiles = c((1-level)/2, 0.5, 1-(1-level)/2))
  
  
  # generate draws
  samp.all <- INLA::inla.posterior.sample(n = n.sample, result = fit, intern = TRUE)
  
  # identify indices of fixed effects and random effects
  fe.idx <- grep(colnames(mm.pop)[1], rownames(samp.all[[1]]$latent))
  fe.idx <- fe.idx:(fe.idx + ncol(mm.pop) - 1)
  re.idx <- grep("domain.struct", x = rownames(samp.all[[1]]$latent))
  
  # aggregate sample predictions
  summary.sample <- function(x) {
    pop.unit.ests <-
      x$latent[re.idx][pop.dat$domain.struct] + mm.pop %*% x$latent[fe.idx] 
    if (family == "binomial") {
      pop.unit.ests <- expit(pop.unit.ests)
    }
    if (!is.null(X.pop.weights)) {
      area.ests <- 
        aggregate(pop.unit.ests * X.pop.weights,  list(domain = pop.dat$domain.struct), sum)
    } else {
      area.ests <- 
        aggregate(pop.unit.ests, list(domain = pop.dat$domain.struct), mean)
    }
    
    return(area.ests[match(1:nrow(domain.table), area.ests[, 1]), 2])
  }
  est.mat <- do.call(cbind, lapply(samp.all, summary.sample))
  out[[paste0(model.method, ".fit")]] <- fit
  out[[paste0(model.method, ".est")]]  <-
    data.frame(domain = domain.table$domain,
               mean = rowMeans(est.mat),
               median = apply(est.mat, 1,
                              function(x) median(x, na.rm = T)),
               var = apply(est.mat, 1, var),
               lower = apply(est.mat, 1,
                             function(x) quantile(x, (1-level)/2, na.rm = T)),
               upper = apply(est.mat, 1,
                             function(x) quantile(x, 1-(1-level)/2, na.rm = T)),
               method = paste0("Unit level model: ", 
                               ifelse(is.null(adj.mat), "IID", "BYM2")))
  if (return.samples) {
    out[[paste0(model.method, ".sample")]]  <- est.mat
  }
  out[[paste0(model.method, ".est")]] <- 
    out[[paste0(model.method, ".est")]][match(out$direct.est$domain, 
                                              out[[paste0(model.method, ".est")]]$domain),]
  attr(out, "method.names") <- c(attr(out, "method.names"), paste0(model.method, ".est"))
  attr(out, "inla.fitted") <- c(attr(out, "inla.fitted"), model.method)
  
  # finish building return value
  out$call <- match.call()
  class(out) <- c("svysae")
  return(out)
}
# GEOSTATISTICAL MODELS --------------------------------------------------------
fitContLGM <- function(formula, 
                       family,
                       cluster,
                       cluster.effect = F,
                       data,
                       mesh,
                       pc.prior.range,
                       pc.prior.sigma,
                       pc.prior.clust,
                       overdisp.mean = 0, 
                       overdisp.prec = 0.4) {
  spde = inla.spde2.pcmatern(mesh = mesh,
                             prior.range = pc.prior.range,
                             prior.sigma = pc.prior.sigma)
  

  
  # Stack
  cluster.var <- all.vars(cluster)
  data$cluster <- data[[cluster.var]]
  control.family <- list()
  if (family == "betabinomial") {
    mf <- model.frame(formula, data)
    data$resp <- model.response(mf, "numeric")
    data <- data |>
      group_by(cluster) |>
      mutate(resp = sum(resp),
             n = n()) |>
      slice(which.min(resp)) |>
      ungroup() 
    formula <- update(formula, resp ~ .)
    control.family <- list(hyper = list(rho = list(param = c(overdisp.mean, overdisp.prec))))
  } 
  if (cluster.effect) {
    cov.frm <- update(formula, NULL ~ . )
    ftxt <- 
      paste(format(update(formula, . ~ intercept)),
            as.character(cov.frm)[-1],
            "f(cluster, model = 'iid', hyper = pc.prior.clust)",
            "f(s, model = spde) - 1", 
            sep = " + ")
  } else {
    cov.frm <- update(formula, NULL ~ .)
    ftxt <- 
      paste(format(update(formula, . ~ intercept)),
            as.character(cov.frm)[-1],
            "f(s, model = spde) - 1", 
            sep = " + ")
  }
  # Mapping
  A = inla.spde.make.A(mesh = mesh,
                       loc = st_coordinates(data))
  covariate.list <- c(list(intercept = 1),
                      setNames(lapply(all.vars(cov.frm), 
                                      function(x) data[[x]]), 
                               all.vars(cov.frm)))
  if (cluster.effect) {
    covariate.list <- c(covariate.list,
                        list(cluster = data$cluster))
  }
  effect.list <- list(list(s = 1:spde$n.spde),
                      covariate.list)
  
  response.frm <- update(formula, . ~ NULL)
  response.list <- setNames(lapply(all.vars(response.frm), 
                                   function(x) data[[x]]), 
                            all.vars(response.frm))
  stk = inla.stack(data = response.list,
                   A = list(A, 1),
                   effects = effect.list)
  
  # Formula
  mod.frm <- as.formula(ftxt)
  # Fit
  if (family == "betabinomial") {
    res.inla.spde = inla(formula = mod.frm,
                         data = inla.stack.data(stk),
                         Ntrials = data$n,
                         family = family,
                         control.family = control.family,
                         control.inla = 
                           list(int.strategy = "ccd"),
                         control.fixed = 
                           list(prec = 1e-4,
                                prec.intercept = 1e-4),
                         control.predictor =
                           list(A = inla.stack.A(stk), compute = TRUE),
                         control.compute = 
                           list(config = TRUE, 
                                return.marginals.predictor = TRUE))
  } else {
    res.inla.spde = inla(formula = mod.frm,
                         data = inla.stack.data(stk),
                         family = family,
                         control.family = control.family,
                         control.inla = 
                           list(int.strategy = "ccd"),
                         control.fixed = 
                           list(prec = 1e-4,
                                prec.intercept = 1e-4),
                         control.predictor =
                           list(A = inla.stack.A(stk), compute = TRUE),
                         control.compute = 
                           list(config = TRUE, 
                                return.marginals.predictor = TRUE))
  }
  
  return(res.inla.spde)
}

smoothContLGM <- function(res.inla,
                          X.pop,
                          domain,
                          mesh,
                          n.sample = 1000,
                          level,
                          cluster.effect = TRUE,
                          X.pop.weights = NULL,
                          return.samples = F) {
  out <- list()
  # Sample from inla object
  post.sample.spde = inla.posterior.sample(n = n.sample, result = res.inla)
  
  
  # Extract different effects
  
  if (cluster.effect) {
    spde.idx <- res.inla$misc$configs$contents$start[4]:(res.inla$misc$configs$contents$start[5]-1)
    cov.idx <- res.inla$misc$configs$contents$start[5]
    clust.sd = vapply(post.sample.spde, function(x) sqrt(1 / x$hyperpar[1]), 
                      FUN.VALUE = numeric(1))
  } else {
    spde.idx <- res.inla$misc$configs$contents$start[3]:(res.inla$misc$configs$contents$start[4]-1)
    cov.idx <- res.inla$misc$configs$contents$start[4]
  }
  cov.idx <- cov.idx:(cov.idx + length(res.inla$names.fixed) - 1)
  
  
  
  spde.sample <-  vapply(post.sample.spde, function(x) x$latent[spde.idx], 
                         FUN.VALUE = numeric(length(spde.idx)))
  cov.sample <- vapply(post.sample.spde, function(x) x$latent[cov.idx], 
                       FUN.VALUE = numeric(length(cov.idx)))
  cov.ftxt <- paste("~", paste(colnames(res.inla$model.matrix), collapse= " + "), "-1")
  cov.frm <- as.formula(cov.ftxt)
  X.pop$intercept <- 1
  
  domain.var <- all.vars(domain)
  if (length(domain.var) > 2) {
    stop("Only suitable for two levels of domains and subdomains")
  } 
  X.pop$domain <- X.pop[[domain.var[1]]]
  if (length(domain.var) != 1) {
    warning("Assuming subdomains are nested.")
    X.pop$subdomain <- X.pop[[domain.var[2]]]
  } 
  # Get predictions for each domain separately
  all.domains <- sort(unique(X.pop$domain))
  predictDomain <- function(dom) {

    X.pop.dom <- X.pop[X.pop$domain == dom, ]
    mm.pop.dom <- model.matrix(cov.frm, X.pop.dom)
    Amap.dom <- inla.spde.make.A(mesh = mesh, loc = st_coordinates(X.pop.dom))
    pop.unit.ests.dom <- 
      Amap.dom %*% spde.sample + 
      mm.pop.dom %*% cov.sample
    if (cluster.effect) {
      pop.unit.ests.dom <- 
        pop.unit.ests.dom + 
        matrix(rnorm(length(pop.unit.ests.dom), sd = clust.sd),
               nrow = nrow(pop.unit.ests.dom))
    }
    pop.unit.ests.dom <- as.matrix(pop.unit.ests.dom)
    if (res.inla$.args$family == "binomial" | res.inla$.args$family == "betabinomial") {
      pop.unit.ests.dom <- expit(pop.unit.ests.dom)
    }
    out <- list()
    if (!is.null(X.pop.weights)) {
      X.pop.wts.dom <- X.pop.weights[X.pop$domain == dom]
      out$domain.ests <- 
        colSums(pop.unit.ests.dom * X.pop.wts.dom)
    } else {
      out$domain.ests <- 
        colMeans(pop.unit.ests.dom)
    }
    out$domain.name <- dom
    if (length(domain.var) != 1) {
      
      if (!is.null(X.pop.weights)) {
        total.wts.subdom <- 
          aggregate(X.pop.wts.dom, 
                    list(subdomain = X.pop.dom$subdomain), sum)
        X.pop.wts.subdom <- 
          X.pop.wts.dom / 
          total.wts.subdom[match(X.pop.dom$subdomain, total.wts.subdom[, 1]), 2]
        out$subdomain.ests <- 
          aggregate(as.data.frame(pop.unit.ests.dom * X.pop.wts.subdom),  
                    list(subdomain = X.pop.dom$subdomain), sum)
      } else {
        out$subdomain.ests <-
          aggregate(as.data.frame(pop.unit.ests.dom), 
                    list(subdomain = X.pop.dom$subdomain), mean)
      }
      out$subdomain.names <- out$subdomain.ests[, 1]
      out$subdomain.ests <- as.matrix(out$subdomain.ests[, -1])
    }
  
    return(out)
  }
  all.ests <- lapply(all.domains, predictDomain)
  est.mat <- do.call(rbind, lapply(all.ests, function(x) x$domain.ests))
  model.method <- paste0(res.inla$.args$family, ".spde.lgm")
  if (return.samples) {
    out[[paste0(model.method, ".sample")]]  <- est.mat
  }
  out[[paste0(model.method, ".est")]]  <-
      data.frame(domain = all.domains,
                 mean = rowMeans(est.mat),
                 median = apply(est.mat, 1,
                                function(x) median(x, na.rm = T)),
                 var = apply(est.mat, 1, var),
                 lower = apply(est.mat, 1,
                               function(x) quantile(x, (1-level)/2, na.rm = T)),
                 upper = apply(est.mat, 1,
                               function(x) quantile(x, 1-(1-level)/2, na.rm = T)),
                 method = paste0(res.inla$.args$family, " SPDE LGM"))
  
  
  est.mat.sub <- do.call(rbind, lapply(all.ests, function(x) x$subdomain.ests))
  if (return.samples) {
    out[[paste0(model.method, ".sample.subdomain")]]  <- est.mat.sub 
  }
  out[[paste0(model.method, ".est.subdomain")]]  <-
    data.frame(domain = do.call(c, lapply(all.ests, function(x) x$subdomain.names)),
               mean = rowMeans(est.mat.sub),
               median = apply(est.mat.sub, 1,
                              function(x) median(x, na.rm = T)),
               var = apply(est.mat.sub, 1, var),
               lower = apply(est.mat.sub, 1,
                             function(x) quantile(x, (1-level)/2, na.rm = T)),
               upper = apply(est.mat.sub, 1,
                             function(x) quantile(x, 1-(1-level)/2, na.rm = T)),
               method = paste0(res.inla$.args$family, " SPDE LGM"))
  
  return(out)
}


simulateFrame <- function(grid, easpa, easize) {
  grid <- grid |> mutate(cell_id = 1:nrow(grid))
  sample_frame <- grid |> 
    st_set_geometry(NULL) |>
    left_join(easpa, by = c("admin1_name", "urban")) |>
    left_join(easize, by = c("admin1_name", "urban")) |>
    group_by(admin1_name, urban) |>
    mutate(pi = cen_pop / sum(cen_pop) * n_ea,
           samp_ind = UPmultinomial(pi)) |>
    ungroup()
  sample_frame <- 
    sample_frame[rep(sample_frame$cell_id, times = sample_frame$samp_ind),]
  sample_frame <- sample_frame |>
    group_by(admin1_name) |>
    mutate(adm1_pop_weight = size / sum(size)) |>
    ungroup() |>
    group_by(admin2_name) |>
    mutate(adm2_pop_weight = size / sum(size)) |>
    ungroup() 
  sample_frame <- grid |>
    select(cell_id, geometry) |> right_join(sample_frame, by = "cell_id")
  return(sample_frame)
}
# PLOTTING FUNCTIONS --------------------------------------------------------
compareEstimates <- function(x,
                             posterior.sample = NULL,
                             title = NULL) {
  
  
  x_att <- attributes(x)
  if (is.null(posterior.sample)) {
    sample_list <- x[paste0(x_att$inla.fitted, ".sample")]
  } else {
    sample_list <- list(posterior.sample)
  }
  for (i in seq_along(sample_list)) {
    current_samples <- sample_list[[i]]
    domain_medians <- apply(current_samples, 1, median)
    median_order <- order(domain_medians)
    
    current_samples <- current_samples[median_order, ]
    # get indices for all possible pairs of admin 1 areas
    domain_comb <- combn(length(x_att$domain.names), 2)
    
    plot_dat <- data.frame(
      Domain1 = x_att$domain.names[median_order][domain_comb[1,]],
      Domain2 = x_att$domain.names[median_order][domain_comb[2,]]
    )
    plot_dat$Domain1 <- factor(plot_dat$Domain1,
                               x_att$domain.names[median_order])
    plot_dat$Domain2 <- factor(plot_dat$Domain2,
                               x_att$domain.names[median_order])
    plot_dat$Prob = apply(
      domain_comb, 2, 
      function(x) mean(current_samples[x[1],] > current_samples[x[2],])
    )
    plot_dat$est = NA
    
    median_dat <- data.frame(
      Domain1 = x_att$domain.names[median_order], 
      Domain2 = "Median", 
      Prob = NA,
      est = domain_medians[median_order]
      
    )
    extra_cols <- c("", "Median", "Interval")
    # combine into single tibble for plotting
    plot_dat <- rbind(plot_dat, median_dat)
    
    plot_dat$Domain2 <- 
      factor(plot_dat$Domain2, levels = c(levels(plot_dat$Domain1), extra_cols))
    
    g_heat <- ggplot2::ggplot(data = plot_dat,
                              ggplot2::aes(x = Domain2, y = Domain1, fill = Prob)) + 
      ggplot2::theme_minimal() + 
      ggplot2::theme(legend.position="bottom", 
                     axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1, size = 10), 
                     panel.grid.major = ggplot2::element_blank(), 
                     panel.grid.minor = ggplot2::element_blank(),
                     text = ggplot2::element_text(size = 10)) +
      ggplot2::geom_tile() + 
      ggplot2::geom_text(
        ggplot2::aes(
          label = ifelse(is.na(Prob), "", sprintf("%0.2f", round(Prob, digits = 2))),
          color = Prob > .5
        ), 
        size = 5
      ) + 
      ggplot2::geom_text(
        data = plot_dat,
        ggplot2::aes(
          label = ifelse(is.na(est), "", sprintf("%0.2f", round(est, digits = 2)))
        ),
        size = 5
      ) +
      ggplot2::coord_equal() + 
      ggplot2::scale_fill_viridis_c(name = "Probability of Domain 1 > Domain 2",
                                    na.value = "white", limits = c(0, .5)) + 
      ggplot2::scale_color_manual(values = c("white", "black"),  guide = "none") +
      ggplot2::labs(x = "Domain 2",
                    y = "Domain 1")
    if (!is.null(title)) {
      g_heat <- g_heat + ggplot2::labs(title = title)
    } else if (is.null(posterior.sample)) {
      g_heat <-  g_heat + ggplot2::labs(title = x_att$inla.fitted[i])
    }
    suppressWarnings(print(g_heat))
  }
}


