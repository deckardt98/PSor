#' @title Compute the principal causal effects with an assumed odds ratio sensitivity parameter
#'
#' @description
#' `PSor.fit` is the main function of the package. It estimates the principal causal effects
#' under principal stratification, assuming an odds ratio sensitivity parameter.
#' The function returns point estimates, standard error estimates, and the confidence intervals for
#' both parametric Conditionally Doubly Robust (CDR) estimators and Debiased Machine Learning (DML) estimators.
#'
#' @param out.formula A formula for the outcome model, e.g., `Y ~ X1 + X2`.
#' @param ps.formula A formula for the principal score model (for the intermediate outcome D), e.g., `D ~ X1 + X2`.
#' @param pro.formula A formula for the propensity score model (for the treatment Z), e.g., `Z ~ X1 + X2`.
#' @param df The data frame containing the variables.
#' @param out.name A character string specifying the name of the outcome variable (Y).
#' @param int.name A character string specifying the name of the intermediate variable (D).
#' @param trt.name A character string specifying the name of the treatment variable (Z).
#' @param cov.names A character vector of covariate names.
#' @param or The odds ratio sensitivity parameter. Can be a single numeric value, a vector of values corresponding to each individual, or `Inf` for the monotonicity assumption.
#' @param SLmethods A character vector of learning algorithms to be used by `SuperLearner` for the DML estimators Defaults to `c("SL.glm", "SL.rpart", "SL.nnet")`.
#' @param n.fold The number of folds for cross-fitting in the DML procedure.
#' @param scale The scale for the causal estimand. Can be "RD" (Risk Difference), "RR" (Risk Ratio), or "OR" (Odds Ratio). Defaults to "RD".
#' @param alpha The significance level for constructing confidence intervals. Defaults to 0.05.
#'
#' @return
#' A data frame containing the point estimates, standard errors, and confidence intervals for the causal effects within each principal stratum ("Always-Takers", "Compliers", "Never-Takers", and "Defiers" if applicable). Both CDR and DML results are provided.
#'
#' @importFrom dplyr %>% rename select all_of
#' @importFrom stats as.formula binomial coef family gaussian glm lm plogis predict qnorm update vcov
#' @importFrom geex m_estimate setup_deriv_control grab_design_matrix grab_fixed_formula grab_psiFUN
#' @importFrom SuperLearner SuperLearner
#' @importFrom caret createFolds
#' @importFrom numDeriv jacobian
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # This is a conceptual example, as the function requires specific data.
#' # Assume 'my_data' is a data frame with variables Y, D, Z, X1, X2.
#'
#' results <- PSor.fit(
#'   out.formula = Y ~ X1 + X2,
#'   ps.formula = D ~ X1 + X2,
#'   pro.formula = Z ~ X1 + X2,
#'   df = my_data,
#'   out.name = "Y",
#'   int.name = "D",
#'   trt.name = "Z",
#'   cov.names = c("X1", "X2"),
#'   or = Inf, # Monotonicity assumption
#'   n.fold = 5,
#'   scale = "RD"
#' )
#'
#' print(results)
#' }
PSor.fit <- function(out.formula, ps.formula, pro.formula,
                     df, out.name, int.name, trt.name, cov.names,
                     or, SLmethods = c("SL.glm", "SL.rpart", "SL.nnet"), n.fold, scale = "RD", alpha = 0.05){
  ## Rename the treatment, intermediate outcome, and final outcome variables for simplicity
  df <- df %>%
    rename(
      Z = !!trt.name,
      D = !!int.name,
      Y = !!out.name)
  n <- nrow(df)
  Z <- df$Z
  D <- df$D
  Y <- df$Y
  X <- df %>% select(all_of(cov.names))
  ## add or to the data frame
  if (length(or)>1){
    df$or.fit = or
  } else {
    df$or.fit = rep(or, n)
  }
  ## Rename model formulas for simplicity
  Y_formula <- update(out.formula, as.formula(paste("Y", "~ .")))
  D_formula <- update(ps.formula, as.formula(paste("D", "~ .")))
  Z_formula <- update(pro.formula, as.formula(paste("Z", "~ .")))

  ################ Estimate nuisance functions by parametric approaches ##############
  ## Estimate propensity score by logistic regression
  propen_m <- glm(Z_formula, data = df, family = binomial)
  proscore <- predict(propen_m, newdata=df, type="response")

  ## Estimate principal scores by logistic regressions
  ps_m1 <- glm(D_formula, data = df, subset = (Z==1), family = binomial)
  ps_m0 <- glm(D_formula, data = df, subset = (Z==0), family = binomial)
  prinscore1 <-  predict(ps_m1, newdata = df, type="response")
  prinscore0 <-  predict(ps_m0, newdata = df, type="response")
  prinscore <- list(prinscore1,prinscore0)

  # Estimate outcome mean
  om_m11 <- lm(Y_formula, data = df, subset = (Z==1)&(D==1))
  om_m01 <- lm(Y_formula, data = df, subset = (Z==0)&(D==1))
  om_m10 <- lm(Y_formula, data = df, subset = (Z==1)&(D==0))
  om_m00 <- lm(Y_formula, data = df, subset = (Z==0)&(D==0))
  om11 <-  predict(om_m11, newdata = df, type="response")
  om01 <-  predict(om_m01, newdata = df, type="response")
  om10 <-  predict(om_m10, newdata = df, type="response")
  om00 <-  predict(om_m00, newdata = df, type="response")
  om <- list(om11,om01,om10,om00)

  create_models.11 <- function(data, pi_formula, p0_formula, p1_formula, m01_formula, m11_formula) {
    list(
      pi_model = glm(pi_formula, data = data, family = binomial),
      p0_model = glm(p0_formula, subset = (Z == 0), data = data, family = binomial),
      p1_model = glm(p1_formula, subset = (Z == 1), data = data, family = binomial),
      m01_model = glm(m01_formula, subset = (Z == 0) & (D == 1), data = data, family = gaussian),
      m11_model = glm(m11_formula, subset = (Z == 1) & (D == 1), data = data, family = gaussian)
    )
  }
  create_models.01 <- function(data, pi_formula, p0_formula, p1_formula, m00_formula, m11_formula) {
    list(
      pi_model = glm(pi_formula, data = data, family = binomial),
      p0_model = glm(p0_formula, subset = (Z == 0), data = data, family = binomial),
      p1_model = glm(p1_formula, subset = (Z == 1), data = data, family = binomial),
      m00_model = glm(m00_formula, subset = (Z == 0) & (D == 0), data = data, family = gaussian),
      m11_model = glm(m11_formula, subset = (Z == 1) & (D == 1), data = data, family = gaussian)
    )
  }
  create_models.00 <- function(data, pi_formula, p0_formula, p1_formula, m00_formula, m10_formula) {
    list(
      pi_model = glm(pi_formula, data = data, family = binomial),
      p0_model = glm(p0_formula, subset = (Z == 0), data = data, family = binomial),
      p1_model = glm(p1_formula, subset = (Z == 1), data = data, family = binomial),
      m00_model = glm(m00_formula, subset = (Z == 0) & (D == 0), data = data, family = gaussian),
      m10_model = glm(m10_formula, subset = (Z == 1) & (D == 0), data = data, family = gaussian)
    )
  }
  create_models.10 <- function(data, pi_formula, p0_formula, p1_formula, m01_formula, m10_formula) {
    list(
      pi_model = glm(pi_formula, data = data, family = binomial),
      p0_model = glm(p0_formula, subset = (Z == 0), data = data, family = binomial),
      p1_model = glm(p1_formula, subset = (Z == 1), data = data, family = binomial),
      m01_model = glm(m01_formula, subset = (Z == 0) & (D == 1), data = data, family = gaussian),
      m10_model = glm(m10_formula, subset = (Z == 1) & (D == 0), data = data, family = gaussian)
    )
  }
  formula.all.11 <-  list(pi = Z_formula, p0 = D_formula, p1 = D_formula, m01 = Y_formula, m11 = Y_formula)
  formula.all.01 <- list(pi = Z_formula, p0 = D_formula, p1 = D_formula, m00 = Y_formula, m11 = Y_formula)
  formula.all.00 <- list(pi = Z_formula, p0 = D_formula, p1 = D_formula, m00 = Y_formula, m10 = Y_formula)
  formula.all.10 <- list(pi = Z_formula, p0 = D_formula, p1 = D_formula, m01 = Y_formula, m10 = Y_formula)
  ###################################################
  ############ Parametric CDR estimator #############
  ###################################################
  if (length(or) == 1 & or[1] == Inf) {
    ######### or = inifity, under monotonicity ########
    ### Point estimates
    point.11.all <- ptrc_point11(Z,D,Y,proscore,prinscore,om)
    point.par.11.1 <- point.11.all[1]
    point.par.11.0 <- point.11.all[2]
    point.01.all <- ptrc_point01(Z,D,Y,proscore,prinscore,om)
    point.par.01.1 <- point.01.all[1]
    point.par.01.0 <- point.01.all[2]
    point.00.all <- ptrc_point00(Z,D,Y,proscore,prinscore,om)
    point.par.00.1 <- point.00.all[1]
    point.par.00.0 <- point.00.all[2]
    if (scale == "RD"){
      point.par.11 <- point.11.all[3]
      point.par.01 <- point.01.all[3]
      point.par.00 <- point.00.all[3]
    } else if (scale == "RR") {
      point.par.11 <- log(point.par.11.1)-log(point.par.11.0)
      point.par.01 <- log(point.par.01.1)-log(point.par.01.0)
      point.par.00 <- log(point.par.00.1)-log(point.par.00.0)
    } else if (scale == "OR") {
      point.par.11 <- log(point.par.11.1/(1-point.par.11.1))-log(point.par.11.0/(1-point.par.11.0))
      point.par.01 <- log(point.par.01.1/(1-point.par.01.1))-log(point.par.01.0/(1-point.par.01.0))
      point.par.00 <- log(point.par.00.1/(1-point.par.00.1))-log(point.par.00.0/(1-point.par.00.0))
    }
    ### Roots of the joint estimating equations
    root.est.11 = c(as.vector(coef(propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)),
                    as.vector(coef(om_m01)), as.vector(coef(om_m11)), point.par.11.1, point.par.11.0, point.par.11)
    root.est.01 = c(as.vector(coef(propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)),
                    as.vector(coef(om_m00)), as.vector(coef(om_m11)), point.par.01.1, point.par.01.0, point.par.01)
    root.est.00 = c(as.vector(coef(propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)),
                    as.vector(coef(om_m00)), as.vector(coef(om_m10)), point.par.00.1, point.par.00.0, point.par.00)
    ### Sandwich variance estimator
    estimate.sand.se.11 <- function(data, pi_formula, p0_formula, p1_formula, m01_formula, m11_formula, roots) {
      models <- create_models.11(data, pi_formula, p0_formula, p1_formula, m01_formula, m11_formula)
      m_estimate(
        estFUN = eq_tr11,
        data = data,
        outer_args = list(models = models, scale = scale),
        compute_roots = FALSE,
        roots = roots,
        deriv_control = setup_deriv_control(
          FUN = custom_jacobian,
          method = "simple"
        )
      )
    }
    estimate.sand.se.01 <- function(data, pi_formula, p0_formula, p1_formula, m00_formula, m11_formula, roots) {
      models <- create_models.01(data, pi_formula, p0_formula, p1_formula, m00_formula, m11_formula)
      m_estimate(
        estFUN = eq_tr01,
        data = data,
        outer_args = list(models = models, scale = scale),
        compute_roots = FALSE,
        roots = roots,
        deriv_control = setup_deriv_control(
          FUN = custom_jacobian,
          method = "simple"
        )
      )
    }
    estimate.sand.se.00 <- function(data, pi_formula, p0_formula, p1_formula, m00_formula, m10_formula, roots) {
      models <- create_models.00(data, pi_formula, p0_formula, p1_formula, m00_formula, m10_formula)
      m_estimate(
        estFUN = eq_tr00,
        data = data,
        outer_args = list(models = models, scale = scale),
        compute_roots = FALSE,
        roots = roots,
        deriv_control = setup_deriv_control(
          FUN = custom_jacobian,
          method = "simple"
        )
      )
    }

    ### Compute the standard error
    estimate.11 <- estimate.sand.se.11(df, formula.all.11$pi, formula.all.11$p0, formula.all.11$p1,
                                       formula.all.11$m01, formula.all.11$m11, root.est.11)
    estimate.01 <- estimate.sand.se.01(df, formula.all.01$pi, formula.all.01$p0, formula.all.01$p1,
                                       formula.all.01$m00, formula.all.01$m11, root.est.01)
    estimate.00 <- estimate.sand.se.00(df, formula.all.00$pi, formula.all.00$p0, formula.all.00$p1,
                                       formula.all.00$m00, formula.all.00$m10, root.est.00)

    cov_matrix.11 <- estimate.11@vcov
    cov_matrix.01 <- estimate.01@vcov
    cov_matrix.00 <- estimate.00@vcov

    se.par.11 <- sqrt(cov_matrix.11[ncol(cov_matrix.11), ncol(cov_matrix.11)])
    se.par.01 <- sqrt(cov_matrix.01[ncol(cov_matrix.01), ncol(cov_matrix.01)])
    se.par.00 <- sqrt(cov_matrix.00[ncol(cov_matrix.00), ncol(cov_matrix.00)])

  } else {
    ########### under non-monotonicity ############
    ### Point estimate
    point.11.all <- pqrc_point11(Z, D, Y, proscore, prinscore, om, df$or.fit)
    point.par.11.1 <- point.11.all[1]
    point.par.11.0 <- point.11.all[2]
    point.01.all <- pqrc_point01(Z, D, Y, proscore, prinscore, om, df$or.fit)
    point.par.01.1 <- point.01.all[1]
    point.par.01.0 <- point.01.all[2]
    point.00.all <- pqrc_point00(Z, D, Y, proscore, prinscore, om, df$or.fit)
    point.par.00.1 <- point.00.all[1]
    point.par.00.0 <- point.00.all[2]
    point.10.all <- pqrc_point10(Z, D, Y, proscore, prinscore, om, df$or.fit)
    point.par.10.1 <- point.10.all[1]
    point.par.10.0 <- point.10.all[2]
    if (scale == "RD"){
      point.par.11 <- point.11.all[3]
      point.par.10 <- point.10.all[3]
      point.par.01 <- point.01.all[3]
      point.par.00 <- point.00.all[3]
    } else if (scale == "RR") {
      point.par.11 <- log(point.par.11.1)-log(point.par.11.0)
      point.par.10 <- log(point.par.10.1)-log(point.par.10.0)
      point.par.01 <- log(point.par.01.1)-log(point.par.01.0)
      point.par.00 <- log(point.par.00.1)-log(point.par.00.0)
    } else if (scale == "OR") {
      point.par.11 <- log(point.par.11.1/(1-point.par.11.1))-log(point.par.11.0/(1-point.par.11.0))
      point.par.10 <- log(point.par.10.1/(1-point.par.10.1))-log(point.par.10.0/(1-point.par.10.0))
      point.par.01 <- log(point.par.01.1/(1-point.par.01.1))-log(point.par.01.0/(1-point.par.01.0))
      point.par.00 <- log(point.par.00.1/(1-point.par.00.1))-log(point.par.00.0/(1-point.par.00.0))
    }
    ### Roots of the joint estimating equations
    root.est.11 = c(as.vector(coef(propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)),
                    as.vector(coef(om_m01)), as.vector(coef(om_m11)), point.par.11.1, point.par.11.0, point.par.11)
    root.est.01 = c(as.vector(coef(propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)),
                    as.vector(coef(om_m00)), as.vector(coef(om_m11)), point.par.01.1, point.par.01.0, point.par.01)
    root.est.00 = c(as.vector(coef(propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)),
                    as.vector(coef(om_m00)), as.vector(coef(om_m10)), point.par.00.1, point.par.00.0, point.par.00)
    root.est.10 = c(as.vector(coef(propen_m)), as.vector(coef(ps_m0)), as.vector(coef(ps_m1)),
                    as.vector(coef(om_m01)), as.vector(coef(om_m10)), point.par.10.1, point.par.10.0, point.par.10)
    ### Sandwich variance estimates
    estimate.sand.se.11 <- function(data, pi_formula, p0_formula, p1_formula, m01_formula, m11_formula, roots) {
      models <- create_models.11(data, pi_formula, p0_formula, p1_formula, m01_formula, m11_formula)
      m_estimate(
        estFUN = eq_qr11,
        data = data,
        outer_args = list(models = models, scale = scale),
        compute_roots = FALSE,
        roots = roots,
        deriv_control = setup_deriv_control(
          FUN = custom_jacobian,
          method = "simple"
        )
      )
    }
    estimate.sand.se.01 <- function(data, pi_formula, p0_formula, p1_formula, m00_formula, m11_formula, roots) {
      models <- create_models.01(data, pi_formula, p0_formula, p1_formula, m00_formula, m11_formula)
      m_estimate(
        estFUN = eq_qr01,
        data = data,
        outer_args = list(models = models, scale = scale),
        compute_roots = FALSE,
        roots = roots,
        deriv_control = setup_deriv_control(
          FUN = custom_jacobian,
          method = "simple"
        )
      )
    }
    estimate.sand.se.00 <- function(data, pi_formula, p0_formula, p1_formula, m00_formula, m10_formula, roots) {
      models <- create_models.00(data, pi_formula, p0_formula, p1_formula, m00_formula, m10_formula)
      m_estimate(
        estFUN = eq_qr00,
        data = data,
        outer_args = list(models = models, scale = scale),
        compute_roots = FALSE,
        roots = roots,
        deriv_control = setup_deriv_control(
          FUN = custom_jacobian,
          method = "simple"
        )
      )
    }
    estimate.sand.se.10 <- function(data, pi_formula, p0_formula, p1_formula, m01_formula, m10_formula, roots) {
      models <- create_models.10(data, pi_formula, p0_formula, p1_formula, m01_formula, m10_formula)
      m_estimate(
        estFUN = eq_qr10,
        data = data,
        outer_args = list(models = models, scale = scale),
        compute_roots = FALSE,
        roots = roots,
        deriv_control = setup_deriv_control(
          FUN = custom_jacobian,
          method = "simple"
        )
      )
    }

    ### Compute the standard error
    estimate.11 <- estimate.sand.se.11(df, formula.all.11$pi, formula.all.11$p0, formula.all.11$p1,
                                       formula.all.11$m01, formula.all.11$m11, root.est.11)
    estimate.01 <- estimate.sand.se.01(df, formula.all.01$pi, formula.all.01$p0, formula.all.01$p1,
                                       formula.all.01$m00, formula.all.01$m11, root.est.01)
    estimate.00 <- estimate.sand.se.00(df, formula.all.00$pi, formula.all.00$p0, formula.all.00$p1,
                                       formula.all.00$m00, formula.all.00$m10, root.est.00)
    estimate.10 <- estimate.sand.se.10(df, formula.all.10$pi, formula.all.10$p0, formula.all.10$p1,
                                       formula.all.10$m01, formula.all.10$m10, root.est.10)

    cov_matrix.11 <- estimate.11@vcov
    cov_matrix.01 <- estimate.01@vcov
    cov_matrix.00 <- estimate.00@vcov
    cov_matrix.10 <- estimate.10@vcov

    se.par.11 <- sqrt(cov_matrix.11[ncol(cov_matrix.11), ncol(cov_matrix.11)])
    se.par.01 <- sqrt(cov_matrix.01[ncol(cov_matrix.01), ncol(cov_matrix.01)])
    se.par.00 <- sqrt(cov_matrix.00[ncol(cov_matrix.00), ncol(cov_matrix.00)])
    se.par.10 <- sqrt(cov_matrix.10[ncol(cov_matrix.10), ncol(cov_matrix.10)])
  }

  if (scale == "RD"){
    par.out.11 = point.par.11
    par.out.01 = point.par.01
    par.out.00 = point.par.00
    par.ci.u.11 = par.out.11+qnorm(1-alpha/2)*se.par.11
    par.ci.u.01 = par.out.01+qnorm(1-alpha/2)*se.par.01
    par.ci.u.00 = par.out.00+qnorm(1-alpha/2)*se.par.00
    par.ci.l.11 = par.out.11-qnorm(1-alpha/2)*se.par.11
    par.ci.l.01 = par.out.01-qnorm(1-alpha/2)*se.par.01
    par.ci.l.00 = par.out.00-qnorm(1-alpha/2)*se.par.00
    if(!(length(or) == 1 & or[1] == Inf)){
      par.out.10 = point.par.10
      par.ci.u.10 = par.out.10+qnorm(1-alpha/2)*se.par.10
      par.ci.l.10 = par.out.10-qnorm(1-alpha/2)*se.par.10
    }
  } else if (scale == "RR"|scale == "OR"){
    par.out.11 = exp(point.par.11)
    par.out.01 = exp(point.par.01)
    par.out.00 = exp(point.par.00)
    par.ci.u.11 = exp(point.par.11+qnorm(1-alpha/2)*se.par.11)
    par.ci.u.01 = exp(point.par.01+qnorm(1-alpha/2)*se.par.01)
    par.ci.u.00 = exp(point.par.00+qnorm(1-alpha/2)*se.par.00)
    par.ci.l.11 = exp(point.par.11-qnorm(1-alpha/2)*se.par.11)
    par.ci.l.01 = exp(point.par.01-qnorm(1-alpha/2)*se.par.01)
    par.ci.l.00 = exp(point.par.00-qnorm(1-alpha/2)*se.par.00)
    if(!(length(or) == 1 & or[1] == Inf)){
      par.out.10 = exp(point.par.10)
      par.ci.u.10 = exp(point.par.10+qnorm(1-alpha/2)*se.par.10)
      par.ci.l.10 = exp(point.par.10-qnorm(1-alpha/2)*se.par.10)
    }
  }


  #########################################################
  ########################### DML #########################
  #########################################################
  # Define a function to apply transformations and fit SuperLearner
  fit_superlearner <- function(Y, X, indices, family, SLmethods, i, id) {
    SuperLearner(Y[indices], X = X[indices, ], newX = X[which(id == i),], family = family,
                 SL.library = SLmethods, cvControl = list(V = n.fold))$SL.predict
  }

  # Function to execute for each fold
  process_fold <- function(i, idc) {
    # Define indices for filtering
    idx_F <- which(idc != i)
    idx_Z1 <- which(idc != i & Z == 1)
    idx_Z0 <- which(idc != i & Z == 0)
    idx_Z1D1 <- which(idc != i & Z == 1 & D == 1)
    idx_Z0D1 <- which(idc != i & Z == 0 & D == 1)
    idx_Z1D0 <- which(idc != i & Z == 1 & D == 0)
    idx_Z0D0 <- which(idc != i & Z == 0 & D == 0)

    # Estimate nuisance functions using SuperLearner
    ml_proscore <- fit_superlearner(Z, X, idx_F, family = "binomial", SLmethods, i, idc)
    ml_ps_m1 <- fit_superlearner(D, X, idx_Z1, family = "binomial", SLmethods, i, idc)
    ml_ps_m0 <- fit_superlearner(D, X, idx_Z0, family = "binomial", SLmethods, i, idc)
    ml_om_m11 <- fit_superlearner(Y, X, idx_Z1D1, family = gaussian(), SLmethods, i, idc)
    ml_om_m01 <- fit_superlearner(Y, X, idx_Z0D1, family = gaussian(), SLmethods, i, idc)
    ml_om_m10 <- fit_superlearner(Y, X, idx_Z1D0, family = gaussian(), SLmethods, i, idc)
    ml_om_m00 <- fit_superlearner(Y, X, idx_Z0D0, family = gaussian(), SLmethods, i, idc)

    if (length(or) == 1 & or[1] == Inf) {
      ## SACE
      DML.all.11 <- ptrc_point11(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore,
                                 list(ml_ps_m1, ml_ps_m0),
                                 list(ml_om_m11, ml_om_m01, ml_om_m10, ml_om_m00))
      DML.var.11 = ml_tr_var11(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore,
                               list(ml_ps_m1, ml_ps_m0),
                               list(ml_om_m11, ml_om_m01, ml_om_m10, ml_om_m00))
      ## CACE
      DML.all.01 <- ptrc_point01(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore,
                                 list(ml_ps_m1, ml_ps_m0),
                                 list(ml_om_m11, ml_om_m01, ml_om_m10, ml_om_m00))
      DML.var.01 = ml_tr_var01(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore,
                               list(ml_ps_m1, ml_ps_m0),
                               list(ml_om_m11, ml_om_m01, ml_om_m10, ml_om_m00))
      ## Never-takers
      DML.all.00 <- ptrc_point00(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore,
                                 list(ml_ps_m1, ml_ps_m0),
                                 list(ml_om_m11, ml_om_m01, ml_om_m10, ml_om_m00))
      DML.var.00 = ml_tr_var00(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore,
                               list(ml_ps_m1, ml_ps_m0),
                               list(ml_om_m11, ml_om_m01, ml_om_m10, ml_om_m00))
      DML.num.1.10 = NA
      DML.num.0.10 = NA
      DML.denom.10 = NA
      DML.var.single.10 = NA
    } else {
      DML.all.11 <- pqrc_point11(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore,
                                 list(ml_ps_m1, ml_ps_m0),
                                 list(ml_om_m11, ml_om_m01, ml_om_m10, ml_om_m00), df$or.fit[-idx_F])
      DML.all.01 <- pqrc_point01(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore,
                                 list(ml_ps_m1, ml_ps_m0),
                                 list(ml_om_m11, ml_om_m01, ml_om_m10, ml_om_m00), df$or.fit[-idx_F])
      DML.all.00 <- pqrc_point00(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore,
                                 list(ml_ps_m1, ml_ps_m0),
                                 list(ml_om_m11, ml_om_m01, ml_om_m10, ml_om_m00), df$or.fit[-idx_F])
      DML.var.11 = ml_qr_var11(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore,
                               list(ml_ps_m1, ml_ps_m0),
                               list(ml_om_m11, ml_om_m01, ml_om_m10, ml_om_m00), df$or.fit[-idx_F])
      DML.var.01 = ml_qr_var01(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore,
                               list(ml_ps_m1, ml_ps_m0),
                               list(ml_om_m11, ml_om_m01, ml_om_m10, ml_om_m00), df$or.fit[-idx_F])
      DML.var.00 = ml_qr_var00(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore,
                               list(ml_ps_m1, ml_ps_m0),
                               list(ml_om_m11, ml_om_m01, ml_om_m10, ml_om_m00), df$or.fit[-idx_F])
      DML.all.10 <- pqrc_point10(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore,
                                 list(ml_ps_m1, ml_ps_m0),
                                 list(ml_om_m11, ml_om_m01, ml_om_m10, ml_om_m00), df$or.fit[-idx_F])
      DML.var.10 = ml_qr_var10(Z[-idx_F], D[-idx_F], Y[-idx_F], ml_proscore,
                               list(ml_ps_m1, ml_ps_m0),
                               list(ml_om_m11, ml_om_m01, ml_om_m10, ml_om_m00), df$or.fit[-idx_F])
      DML.num.1.10 <- DML.all.10[4]
      DML.num.0.10 <- DML.all.10[5]
      DML.num.dif.10 = DML.num.1.10-DML.num.0.10
      DML.denom.10 = DML.all.10[6]
      mu.1.temp.DML.10 = DML.num.1.10/DML.denom.10
      mu.0.temp.DML.10 = DML.num.0.10/DML.denom.10
      if (scale == "RD") {
        DML.var.single.10 = DML.var.10[1]
      } else if (scale == "RR"){
        DML.var.single.10 = DML.var.10[3]/mu.1.temp.DML.10^2+DML.var.10[2]/mu.0.temp.DML.10^2-2*DML.var.10[4]/mu.1.temp.DML.10/mu.0.temp.DML.10
      } else if (scale == "OR"){
        DML.var.single.10 = DML.var.10[3]/(mu.1.temp.DML.10*(1-mu.1.temp.DML.10))^2+
          DML.var.10[2]/(mu.0.temp.DML.10*(1-mu.0.temp.DML.10))^2-2*DML.var.10[4]/mu.1.temp.DML.10/(1-mu.1.temp.DML.10)/mu.0.temp.DML.10/(1-mu.0.temp.DML.10)
      }
    }
    ## SACE
    DML.num.1.11 <- DML.all.11[4]
    DML.num.0.11 <- DML.all.11[5]
    DML.num.dif.11 = DML.num.1.11-DML.num.0.11
    DML.denom.11 = DML.all.11[6]
    mu.1.temp.DML.11 = DML.num.1.11/DML.denom.11
    mu.0.temp.DML.11 = DML.num.0.11/DML.denom.11
    ## CACE
    DML.num.1.01 <- DML.all.01[4]
    DML.num.0.01 <- DML.all.01[5]
    DML.num.dif.01 = DML.num.1.01-DML.num.0.01
    DML.denom.01 = DML.all.01[6]
    mu.1.temp.DML.01 = DML.num.1.01/DML.denom.01
    mu.0.temp.DML.01 = DML.num.0.01/DML.denom.01
    ## Never-takers
    DML.num.1.00 <- DML.all.00[4]
    DML.num.0.00 <- DML.all.00[5]
    DML.num.dif.00 = DML.num.1.00-DML.num.0.00
    DML.denom.00 = DML.all.00[6]
    mu.1.temp.DML.00 = DML.num.1.00/DML.denom.00
    mu.0.temp.DML.00 = DML.num.0.00/DML.denom.00
    if (scale == "RD") {
      DML.var.single.11 = DML.var.11[1]
      DML.var.single.01 = DML.var.01[1]
      DML.var.single.00 = DML.var.00[1]
    } else if (scale == "RR"){
      DML.var.single.11 = DML.var.11[3]/mu.1.temp.DML.11^2+DML.var.11[2]/mu.0.temp.DML.11^2-2*DML.var.11[4]/mu.1.temp.DML.11/mu.0.temp.DML.11
      DML.var.single.01 = DML.var.01[3]/mu.1.temp.DML.01^2+DML.var.01[2]/mu.0.temp.DML.01^2-2*DML.var.01[4]/mu.1.temp.DML.01/mu.0.temp.DML.01
      DML.var.single.00 = DML.var.00[3]/mu.1.temp.DML.00^2+DML.var.00[2]/mu.0.temp.DML.00^2-2*DML.var.00[4]/mu.1.temp.DML.00/mu.0.temp.DML.00
    } else if (scale == "OR"){
      DML.var.single.11 = DML.var.11[3]/(mu.1.temp.DML.11*(1-mu.1.temp.DML.11))^2+
        DML.var.11[2]/(mu.0.temp.DML.11*(1-mu.0.temp.DML.11))^2-2*DML.var.11[4]/mu.1.temp.DML.11/(1-mu.1.temp.DML.11)/mu.0.temp.DML.11/(1-mu.0.temp.DML.11)
      DML.var.single.01 = DML.var.01[3]/(mu.1.temp.DML.01*(1-mu.1.temp.DML.01))^2+
        DML.var.01[2]/(mu.0.temp.DML.01*(1-mu.0.temp.DML.01))^2-2*DML.var.01[4]/mu.1.temp.DML.01/(1-mu.1.temp.DML.01)/mu.0.temp.DML.01/(1-mu.0.temp.DML.01)
      DML.var.single.00 = DML.var.00[3]/(mu.1.temp.DML.00*(1-mu.1.temp.DML.00))^2+
        DML.var.00[2]/(mu.0.temp.DML.00*(1-mu.0.temp.DML.00))^2-2*DML.var.00[4]/mu.1.temp.DML.00/(1-mu.1.temp.DML.00)/mu.0.temp.DML.00/(1-mu.0.temp.DML.00)
    }
    ########################### Return results as a vector ######################
    c(DML.num.1.11,DML.num.0.11,DML.denom.11,DML.var.single.11,
      DML.num.1.01,DML.num.0.01,DML.denom.01,DML.var.single.01,
      DML.num.1.00,DML.num.0.00,DML.denom.00,DML.var.single.00,
      DML.num.1.10,DML.num.0.10,DML.denom.10,DML.var.single.10)
  }
  # --- Modifications for retry limit ---
  retry_count <- 0
  max_retries <- 20 # Set a reasonable limit, e.g., 20 times
  results <- NULL # Initialize results to NULL

  ml_error <- TRUE
  n_fold <- n.fold
  # The while loop condition is now updated
  while (ml_error && retry_count < max_retries) {
    idc <- gene_Fi(df,n_fold)
    tryCatch({
      results <- lapply(1:n_fold, process_fold, idc = idc)
      ml_error <- FALSE # Set to FALSE on success to exit the loop

    },
    warning = function(w) {
      message("Warning occurred, retrying. Attempt: ", retry_count + 1)
      retry_count <<- retry_count + 1 # Increment counter on warning
      ml_error <<- TRUE
    },
    error = function(e) {
      message("Error occurred, retrying. Attempt: ", retry_count + 1)
      retry_count <<- retry_count + 1 # Increment counter on error
      ml_error <<- TRUE
    })
  }
  if (is.null(results)) {
    message("Failed to get results after ", max_retries, " attempts. Returning NA.")
    # RETURN NA INSTEAD OF ERRORING
    # This creates a dummy structure to prevent the next lines from failing
    na_results_matrix <- matrix(NA, nrow = ifelse(length(or) == 1 & or[1] == Inf, 3, 4), ncol = 8)
    colnames(na_results_matrix) <- c("CDR.Est", "CDR.SE", "CDR.ci.lower", "CDR.ci.upper",
                                     "DML.Est", "DML.SE", "DML.ci.lower", "DML.ci.upper")
    rownames_vec <- c("Always-Takers (11)", "Compliers (01)", "Never-Takers (00)")
    if (!(length(or) == 1 & or[1] == Inf)) {
      rownames_vec <- c(rownames_vec, "Defiers (10)")
    }
    rownames(na_results_matrix) <- rownames_vec
    return(as.data.frame(na_results_matrix))

  } else {
    DML.num.final.1.11 <- sapply(results, `[`, 1)
    DML.num.final.0.11 <- sapply(results, `[`, 2)
    DML.denom.final.11 = sapply(results, `[`, 3)
    DML.var.final.11 = sapply(results, `[`, 4)

    DML.num.final.1.01 <- sapply(results, `[`, 5)
    DML.num.final.0.01 <- sapply(results, `[`, 6)
    DML.denom.final.01 = sapply(results, `[`, 7)
    DML.var.final.01 = sapply(results, `[`, 8)

    DML.num.final.1.00 <- sapply(results, `[`, 9)
    DML.num.final.0.00 <- sapply(results, `[`, 10)
    DML.denom.final.00 = sapply(results, `[`, 11)
    DML.var.final.00 = sapply(results, `[`, 12)

    DML.num.final.1.10 <- sapply(results, `[`, 13)
    DML.num.final.0.10 <- sapply(results, `[`, 14)
    DML.denom.final.10 = sapply(results, `[`, 15)
    DML.var.final.10 = sapply(results, `[`, 16)

    weig <- as.numeric(table(idc)) / n
    DML.mu.1.out.11 <- (sum(DML.num.final.1.11 * weig, na.rm=TRUE))/(sum(DML.denom.final.11 * weig, na.rm=TRUE))
    DML.mu.0.out.11 <- (sum(DML.num.final.0.11 * weig, na.rm=TRUE))/(sum(DML.denom.final.11 * weig, na.rm=TRUE))
    DML.var.out.11 <- sum(DML.var.final.11 * weig, na.rm=TRUE)/n

    DML.mu.1.out.01 <- (sum(DML.num.final.1.01 * weig, na.rm=TRUE))/(sum(DML.denom.final.01 * weig, na.rm=TRUE))
    DML.mu.0.out.01 <- (sum(DML.num.final.0.01 * weig, na.rm=TRUE))/(sum(DML.denom.final.01 * weig, na.rm=TRUE))
    DML.var.out.01 <- sum(DML.var.final.01 * weig, na.rm=TRUE)/n

    DML.mu.1.out.00 <- (sum(DML.num.final.1.00 * weig, na.rm=TRUE))/(sum(DML.denom.final.00 * weig, na.rm=TRUE))
    DML.mu.0.out.00 <- (sum(DML.num.final.0.00 * weig, na.rm=TRUE))/(sum(DML.denom.final.00 * weig, na.rm=TRUE))
    DML.var.out.00 <- sum(DML.var.final.00 * weig, na.rm=TRUE)/n

    DML.mu.1.out.10 <- (sum(DML.num.final.1.10 * weig, na.rm=TRUE))/(sum(DML.denom.final.10 * weig, na.rm=TRUE))
    DML.mu.0.out.10 <- (sum(DML.num.final.0.10 * weig, na.rm=TRUE))/(sum(DML.denom.final.10 * weig, na.rm=TRUE))
    DML.var.out.10 <- sum(DML.var.final.10 * weig, na.rm=TRUE)/n

    if (scale == "RD"){
      DML.out.11 = DML.mu.1.out.11-DML.mu.0.out.11
      DML.out.01 = DML.mu.1.out.01-DML.mu.0.out.01
      DML.out.00 = DML.mu.1.out.00-DML.mu.0.out.00
      DML.ci.u.11 = DML.out.11+qnorm(1-alpha/2)*sqrt(DML.var.out.11)
      DML.ci.u.01 = DML.out.01+qnorm(1-alpha/2)*sqrt(DML.var.out.01)
      DML.ci.u.00 = DML.out.00+qnorm(1-alpha/2)*sqrt(DML.var.out.00)
      DML.ci.l.11 = DML.out.11-qnorm(1-alpha/2)*sqrt(DML.var.out.11)
      DML.ci.l.01 = DML.out.01-qnorm(1-alpha/2)*sqrt(DML.var.out.01)
      DML.ci.l.00 = DML.out.00-qnorm(1-alpha/2)*sqrt(DML.var.out.00)
      if(!(length(or) == 1 & or[1] == Inf)){
        DML.out.10 = DML.mu.1.out.10-DML.mu.0.out.10
        DML.ci.u.10 = DML.out.10+qnorm(1-alpha/2)*sqrt(DML.var.out.10)
        DML.ci.l.10 = DML.out.10-qnorm(1-alpha/2)*sqrt(DML.var.out.10)
      }
    } else if (scale == "RR"){
      DML.out.11 = DML.mu.1.out.11/DML.mu.0.out.11
      DML.out.01 = DML.mu.1.out.01/DML.mu.0.out.01
      DML.out.00 = DML.mu.1.out.00/DML.mu.0.out.00
      DML.out.log.11 = log(DML.mu.1.out.11)-log(DML.mu.0.out.11)
      DML.out.log.01 = log(DML.mu.1.out.01)-log(DML.mu.0.out.01)
      DML.out.log.00 = log(DML.mu.1.out.00)-log(DML.mu.0.out.00)
      DML.ci.u.11 = exp(DML.out.log.11+qnorm(1-alpha/2)*sqrt(DML.var.out.11))
      DML.ci.u.01 = exp(DML.out.log.01+qnorm(1-alpha/2)*sqrt(DML.var.out.01))
      DML.ci.u.00 = exp(DML.out.log.00+qnorm(1-alpha/2)*sqrt(DML.var.out.00))
      DML.ci.l.11 = exp(DML.out.log.11-qnorm(1-alpha/2)*sqrt(DML.var.out.11))
      DML.ci.l.01 = exp(DML.out.log.01-qnorm(1-alpha/2)*sqrt(DML.var.out.01))
      DML.ci.l.00 = exp(DML.out.log.00-qnorm(1-alpha/2)*sqrt(DML.var.out.00))
      if(!(length(or) == 1 & or[1] == Inf)){
        DML.out.10 = DML.mu.1.out.10/DML.mu.0.out.10
        DML.out.log.10 = log(DML.mu.1.out.10)-log(DML.mu.0.out.10)
        DML.ci.u.10 = exp(DML.out.log.10+qnorm(1-alpha/2)*sqrt(DML.var.out.10))
        DML.ci.l.10 = exp(DML.out.log.10-qnorm(1-alpha/2)*sqrt(DML.var.out.10))
      }
    } else if (scale == "OR"){
      DML.out.11 = DML.mu.1.out.11*(1-DML.mu.0.out.11)/DML.mu.0.out.11/(1-DML.mu.1.out.11)
      DML.out.01 = DML.mu.1.out.01*(1-DML.mu.0.out.01)/DML.mu.0.out.01/(1-DML.mu.1.out.01)
      DML.out.00 = DML.mu.1.out.00*(1-DML.mu.0.out.00)/DML.mu.0.out.00/(1-DML.mu.1.out.00)
      DML.out.log.11 = log(DML.mu.1.out.11/(1-DML.mu.1.out.11))-log(DML.mu.0.out.11/(1-DML.mu.0.out.11))
      DML.out.log.01 = log(DML.mu.1.out.01/(1-DML.mu.1.out.01))-log(DML.mu.0.out.01/(1-DML.mu.0.out.01))
      DML.out.log.00 = log(DML.mu.1.out.00/(1-DML.mu.1.out.00))-log(DML.mu.0.out.00/(1-DML.mu.0.out.00))
      DML.ci.u.11 = exp(DML.out.log.11+qnorm(1-alpha/2)*sqrt(DML.var.out.11))
      DML.ci.u.01 = exp(DML.out.log.01+qnorm(1-alpha/2)*sqrt(DML.var.out.01))
      DML.ci.u.00 = exp(DML.out.log.00+qnorm(1-alpha/2)*sqrt(DML.var.out.00))
      DML.ci.l.11 = exp(DML.out.log.11-qnorm(1-alpha/2)*sqrt(DML.var.out.11))
      DML.ci.l.01 = exp(DML.out.log.01-qnorm(1-alpha/2)*sqrt(DML.var.out.01))
      DML.ci.l.00 = exp(DML.out.log.00-qnorm(1-alpha/2)*sqrt(DML.var.out.00))
      if(!(length(or) == 1 & or[1] == Inf)){
        DML.out.10 = DML.mu.1.out.10*(1-DML.mu.0.out.10)/DML.mu.0.out.10/(1-DML.mu.1.out.10)
        DML.out.log.10 = log(DML.mu.1.out.10/(1-DML.mu.1.out.10))-log(DML.mu.0.out.10/(1-DML.mu.0.out.10))
        DML.ci.u.10 = exp(DML.out.log.10+qnorm(1-alpha/2)*sqrt(DML.var.out.10))
        DML.ci.l.10 = exp(DML.out.log.10-qnorm(1-alpha/2)*sqrt(DML.var.out.10))
      }
    }

    if (length(or) == 1 & or[1] == Inf) {
      out.all.11 <- c(par.out.11, se.par.11, par.ci.l.11, par.ci.u.11,
                      DML.out.11, sqrt(DML.var.out.11), DML.ci.l.11, DML.ci.u.11)
      out.all.01 <- c(par.out.01, se.par.01, par.ci.l.01, par.ci.u.01,
                      DML.out.01, sqrt(DML.var.out.01), DML.ci.l.01, DML.ci.u.01)
      out.all.00 <- c(par.out.00, se.par.00, par.ci.l.00, par.ci.u.00,
                      DML.out.00, sqrt(DML.var.out.00), DML.ci.l.00, DML.ci.u.00)
      results_matrix <- rbind(out.all.11, out.all.01, out.all.00)
      colnames(results_matrix) <- c("CDR.Est", "CDR.SE", "CDR.ci.lower", "CDR.ci.upper",
                                    "DML.Est", "DML.SE", "DML.ci.lower", "DML.ci.upper")
      rownames(results_matrix) <- c("Always-Takers (11)", "Compliers (01)",
                                    "Never-Takers (00)")
      results_table <- as.data.frame(results_matrix)
    } else {
      out.all.11 <- c(par.out.11, se.par.11, par.ci.l.11, par.ci.u.11,
                      DML.out.11, sqrt(DML.var.out.11), DML.ci.l.11, DML.ci.u.11)
      out.all.01 <- c(par.out.01, se.par.01, par.ci.l.01, par.ci.u.01,
                      DML.out.01, sqrt(DML.var.out.01), DML.ci.l.01, DML.ci.u.01)
      out.all.00 <- c(par.out.00, se.par.00, par.ci.l.00, par.ci.u.00,
                      DML.out.00, sqrt(DML.var.out.00), DML.ci.l.00, DML.ci.u.00)
      out.all.10 <- c(par.out.10, se.par.10, par.ci.l.10, par.ci.u.10,
                      DML.out.10, sqrt(DML.var.out.10), DML.ci.l.10, DML.ci.u.10)
      results_matrix <- rbind(out.all.11, out.all.01, out.all.00, out.all.10)
      colnames(results_matrix) <- c("CDR.Est", "CDR.SE", "CDR.ci.lower", "CDR.ci.upper",
                                    "DML.Est", "DML.SE", "DML.ci.lower", "DML.ci.upper")
      rownames(results_matrix) <- c("Always-Takers (11)", "Compliers (01)",
                                    "Never-Takers (00)", "Defiers (10)")
      results_table <- as.data.frame(results_matrix)
    }

    return(round(results_table,3))
  }
}
