#' Internal Helper Functions for PSor Package
#'
#' These functions are not exported for users and are intended for internal use only.
#' They include estimators for different principal strata (always-takers, compliers,
#' never-takers, defiers), variance estimators, and estimating equations for
#' the main `PSor.fit` function.
#'
#' @name internal_helpers
#' @noRd
NULL

######################## Estimating Equations for SACE ###########################

#' Parametric conditionally doubly robust estimator for SACE
#' @param Z Treatment vector.
#' @param D Intermediate outcome vector.
#' @param Y Final outcome vector.
#' @param proscore Propensity scores.
#' @param prinscore Principal scores.
#' @param om Outcome mean predictions.
#' @param theta Odds ratio parameter.
#' @noRd
pqrc_point11 <- function(Z,D,Y,proscore,prinscore,om,theta){
  prinscore1 <- prinscore[[1]]
  prinscore0 <- prinscore[[2]]
  om11 <- om[[1]]
  om01 <- om[[2]]
  om10 <- om[[3]]
  om00 <- om[[4]]

  #Estimate Psi functions
  psi_D1 <- Z/proscore*(D-prinscore1)+prinscore1
  psi_D0 <- (1-Z)/(1-proscore)*(D-prinscore0)+prinscore0
  psi_YD1 <- Z/proscore*(Y*D-om11*prinscore1)+om11*prinscore1
  psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*prinscore0)+om01*prinscore0

  #tau and omega
  ## independence
  ex.ind = prinscore0*prinscore1
  tau.ind <- prinscore0*psi_D1+(psi_D0-prinscore0)*prinscore1
  omega1.ind <- prinscore0*(psi_YD1-om11*psi_D1)+om11*tau.ind
  omega0.ind <- prinscore1*(psi_YD0-om01*psi_D0)+om01*tau.ind
  ## non-independence
  delta <- (1+(theta-1)*(prinscore1+prinscore0))^2-4*theta*(theta-1)*prinscore1*prinscore0
  lambda <- 1+(theta-1)*(prinscore1+prinscore0)-sqrt(delta)
  ex.non.ind = lambda/2/(theta-1)
  tau.non.ind <- ex.non.ind + 1/2/sqrt(delta)*((psi_D0-prinscore0)*(2*theta*prinscore1-lambda)+(psi_D1-prinscore1)*(2*theta*prinscore0-lambda))
  omega1.non.ind <- ex.non.ind/prinscore1*(psi_YD1-om11*psi_D1)+tau.non.ind*om11
  omega0.non.ind <- ex.non.ind/prinscore0*(psi_YD0-om01*psi_D0)+tau.non.ind*om01
  ## combine
  tau = ifelse(theta == 1, tau.ind, tau.non.ind)
  omega1 = ifelse(theta == 1, omega1.ind, omega1.non.ind)
  omega0 = ifelse(theta == 1, omega0.ind, omega0.non.ind)

  #parametric conditionally doubly robust estimator
  pqr1 <- mean(omega1)/mean(tau)
  pqr0 <- mean(omega0)/mean(tau)
  pqrc <- pqr1-pqr0
  denom_est <- mean(tau)
  num_est.1 <- mean(omega1)
  num_est.0 <- mean(omega0)

  return(c(pqr1,pqr0,pqrc,num_est.1,num_est.0,denom_est))
}

#' Variance of the DML estimator for SACE
#' @param Z Treatment vector.
#' @param D Intermediate outcome vector.
#' @param Y Final outcome vector.
#' @param proscore Propensity scores.
#' @param prinscore Principal scores.
#' @param om Outcome mean predictions.
#' @param theta Odds ratio parameter.
#' @noRd
ml_qr_var11 <- function(Z,D,Y,proscore,prinscore,om,theta){
  prinscore1 <- prinscore[[1]]
  prinscore0 <- prinscore[[2]]
  om11 <- om[[1]]
  om01 <- om[[2]]
  om10 <- om[[3]]
  om00 <- om[[4]]

  #Estimate Psi functions
  psi_D1 <- Z/proscore*(D-prinscore1)+prinscore1
  psi_D0 <- (1-Z)/(1-proscore)*(D-prinscore0)+prinscore0
  psi_YD1 <- Z/proscore*(Y*D-om11*prinscore1)+om11*prinscore1
  psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*prinscore0)+om01*prinscore0

  #tau and omega
  ## independence
  ex.ind = prinscore0*prinscore1
  tau.ind <- prinscore0*psi_D1+(psi_D0-prinscore0)*prinscore1
  omega1.ind <- prinscore0*(psi_YD1-om11*psi_D1)+om11*tau.ind
  omega0.ind <- prinscore1*(psi_YD0-om01*psi_D0)+om01*tau.ind
  ## non-independence
  delta <- (1+(theta-1)*(prinscore1+prinscore0))^2-4*theta*(theta-1)*prinscore1*prinscore0
  lambda <- 1+(theta-1)*(prinscore1+prinscore0)-sqrt(delta)
  ex.non.ind = lambda/2/(theta-1)
  tau.non.ind <- ex.non.ind + 1/2/sqrt(delta)*((psi_D0-prinscore0)*(2*theta*prinscore1-lambda)+(psi_D1-prinscore1)*(2*theta*prinscore0-lambda))
  omega1.non.ind <- ex.non.ind/prinscore1*(psi_YD1-om11*psi_D1)+tau.non.ind*om11
  omega0.non.ind <- ex.non.ind/prinscore0*(psi_YD0-om01*psi_D0)+tau.non.ind*om01
  ## combine
  tau = ifelse(theta == 1, tau.ind, tau.non.ind)
  omega1 = ifelse(theta == 1, omega1.ind, omega1.non.ind)
  omega0 = ifelse(theta == 1, omega0.ind, omega0.non.ind)


  #Estimate the denominator
  denom_est <- mean(tau)

  #Compute the variance estimates
  pqr1 <- mean(omega1)/mean(tau)
  pqr0 <- mean(omega0)/mean(tau)
  xi1 <- omega1-pqr1*tau
  xi0 <- omega0-pqr0*tau
  var_single <- mean((xi1-xi0)^2)/denom_est^2
  var.0 = mean(xi0^2)/denom_est^2
  var.1 = mean(xi1^2)/denom_est^2
  cov = mean(xi0*xi1)/denom_est^2

  return(c(var_single,var.0,var.1,cov))
}

#' Parametric triply robust estimator for SACE under monotonicity
#' @param Z Treatment vector.
#' @param D Intermediate outcome vector.
#' @param Y Final outcome vector.
#' @param proscore Propensity scores.
#' @param prinscore Principal scores.
#' @param om Outcome mean predictions.
#' @noRd
ptrc_point11 <- function(Z,D,Y,proscore,prinscore,om){
  prinscore1 <- prinscore[[1]]
  prinscore0 <- prinscore[[2]]
  om11 <- om[[1]]
  om01 <- om[[2]]
  om10 <- om[[3]]
  om00 <- om[[4]]

  #Estimate Psi functions
  psi_D1 <- Z/proscore*(D-prinscore1)+prinscore1
  psi_D0 <- (1-Z)/(1-proscore)*(D-prinscore0)+prinscore0
  psi_YD1 <- Z/proscore*(Y*D-om11*prinscore1)+om11*prinscore1
  psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*prinscore0)+om01*prinscore0

  #parametric triply robust estimator
  ptr1 <- mean(prinscore0/prinscore1*(psi_YD1-om11*psi_D1)+om11*(psi_D0))/mean(psi_D0)
  ptr0 <- mean(psi_YD0)/mean(psi_D0)
  ptrc <- ptr1-ptr0
  denom_est <- mean(psi_D0)
  num_est.1 = mean(prinscore0/prinscore1*(psi_YD1-om11*psi_D1)+om11*(psi_D0))
  num_est.0 = mean(psi_YD0)

  return(c(ptr1,ptr0,ptrc,num_est.1,num_est.0,denom_est))
}

#' Variance of the DML estimator for SACE under monotonicity
#' @param Z Treatment vector.
#' @param D Intermediate outcome vector.
#' @param Y Final outcome vector.
#' @param proscore Propensity scores.
#' @param prinscore Principal scores.
#' @param om Outcome mean predictions.
#' @noRd
ml_tr_var11 <- function(Z,D,Y,proscore,prinscore,om){
  prinscore1 <- prinscore[[1]]
  prinscore0 <- prinscore[[2]]
  om11 <- om[[1]]
  om01 <- om[[2]]
  om10 <- om[[3]]
  om00 <- om[[4]]

  #Estimate Psi functions
  psi_D1 <- Z/proscore*(D-prinscore1)+prinscore1
  psi_D0 <- (1-Z)/(1-proscore)*(D-prinscore0)+prinscore0
  psi_YD1 <- Z/proscore*(Y*D-om11*prinscore1)+om11*prinscore1
  psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*prinscore0)+om01*prinscore0

  #Estimate the denominator
  denom_est <- mean(psi_D0)

  #Compute variance estimates
  ptr1 <- mean(prinscore0/prinscore1*(psi_YD1-om11*psi_D1)+om11*(psi_D0))/mean(psi_D0)
  ptr0 <- mean(psi_YD0)/mean(psi_D0)
  var_single <- mean(((prinscore0/prinscore1*(psi_YD1-om11*psi_D1)+om11*(psi_D0)-ptr1*psi_D0)-(psi_YD0-ptr0*psi_D0))^2)/denom_est^2
  var.0 = mean((psi_YD0-ptr0*psi_D0)^2)/denom_est^2
  var.1 = mean((prinscore0/prinscore1*(psi_YD1-om11*psi_D1)+om11*(psi_D0)-ptr1*psi_D0)^2)/denom_est^2
  cov = mean((prinscore0/prinscore1*(psi_YD1-om11*psi_D1)+om11*(psi_D0)-ptr1*psi_D0)*(psi_YD0-ptr0*psi_D0))/denom_est^2

  return(c(var_single,var.0,var.1,cov))
}

#' Estimating equations for parametric multiply robust estimator for SACE
#' @param data Data frame.
#' @param models List of fitted models.
#' @param scale Causal estimand scale.
#' @noRd
eq_qr11 <- function(data, models, scale){
  Z <- data$Z
  D <- data$D
  Y <- data$Y
  oddsratio = data$or.fit

  Xpi <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$pi_model))
  Xp0 <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$p0_model))
  Xp1 <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$p1_model))
  Xm01 <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$m01_model))
  Xm11 <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$m11_model))

  pi_pos <- 1:ncol(Xpi)
  p0_pos <- (max(pi_pos) + 1):(max(pi_pos) + ncol(Xp0))
  p1_pos <- (max(p0_pos) + 1):(max(p0_pos) + ncol(Xp1))
  m01_pos <- (max(p1_pos) + 1):(max(p1_pos) + ncol(Xm01))
  m11_pos <- (max(m01_pos) + 1):(max(m01_pos) + ncol(Xm11))

  pi_scores <- grab_psiFUN(models$pi_model, data)
  p0_scores <- grab_psiFUN(models$p0_model, data)
  p1_scores <- grab_psiFUN(models$p1_model, data)
  m01_scores <- grab_psiFUN(models$m01_model, data)
  m11_scores <- grab_psiFUN(models$m11_model, data)

  function(theta){
    proscore <- plogis(Xpi %*% theta[pi_pos])
    p0 <- plogis(Xp0 %*% theta[p0_pos])
    p1 <- plogis(Xp1 %*% theta[p1_pos])
    om01 <- Xm01 %*% theta[m01_pos]
    om11 <- Xm11 %*% theta[m11_pos]
    psi_D0 <- (1-Z)/(1-proscore)*(D-p0)+p0
    psi_D1 <- Z/proscore*(D-p1)+p1
    psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*p0)+om01*p0
    psi_YD1 <- Z/proscore*(Y*D-om11*p1)+om11*p1
    ## independence
    ex.ind = p0*p1
    tau.ind <- p0*psi_D1+(psi_D0-p0)*p1
    omega1.ind <- p0*(psi_YD1-om11*psi_D1)+om11*tau.ind
    omega0.ind <- p1*(psi_YD0-om01*psi_D0)+om01*tau.ind
    ## non-independence
    delta <- (1+(oddsratio-1)*(p1+p0))^2-4*oddsratio*(oddsratio-1)*p1*p0
    lambda <- 1+(oddsratio-1)*(p1+p0)-sqrt(delta)
    ex.non.ind = lambda/2/(oddsratio-1)
    tau.non.ind <- ex.non.ind + 1/2/sqrt(delta)*((psi_D0-p0)*(2*oddsratio*p1-lambda)+(psi_D1-p1)*(2*oddsratio*p0-lambda))
    omega1.non.ind <- ex.non.ind/p1*(psi_YD1-om11*psi_D1)+tau.non.ind*om11
    omega0.non.ind <- ex.non.ind/p0*(psi_YD0-om01*psi_D0)+tau.non.ind*om01
    ## combine
    tau = ifelse(oddsratio == 1, tau.ind, tau.non.ind)
    omega1 = ifelse(oddsratio == 1, omega1.ind, omega1.non.ind)
    omega0 = ifelse(oddsratio == 1, omega0.ind, omega0.non.ind)
    if (scale == "RD") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m01_scores(theta[m01_pos])*I(Z==0)*I(D==1),
        m11_scores(theta[m11_pos])*I(Z==1)*I(D==1),
        omega1-theta[max(m11_pos)+1]*tau,
        omega0-theta[max(m11_pos)+2]*tau,
        theta[max(m11_pos)+1]-theta[max(m11_pos)+2]-theta[max(m11_pos)+3])
    } else if (scale == "RR") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m01_scores(theta[m01_pos])*I(Z==0)*I(D==1),
        m11_scores(theta[m11_pos])*I(Z==1)*I(D==1),
        omega1-theta[max(m11_pos)+1]*tau,
        omega0-theta[max(m11_pos)+2]*tau,
        log(theta[max(m11_pos)+1])-log(theta[max(m11_pos)+2])-theta[max(m11_pos)+3])
    } else if (scale == "OR") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m01_scores(theta[m01_pos])*I(Z==0)*I(D==1),
        m11_scores(theta[m11_pos])*I(Z==1)*I(D==1),
        omega1-theta[max(m11_pos)+1]*tau,
        omega0-theta[max(m11_pos)+2]*tau,
        log(theta[max(m11_pos)+1]/(1-theta[max(m11_pos)+1]))-log(theta[max(m11_pos)+2]/(1-theta[max(m11_pos)+2]))-theta[max(m11_pos)+3])
    }
  }
}

#' Estimating equations for parametric triply robust estimator for SACE under monotonicity
#' @param data Data frame.
#' @param models List of fitted models.
#' @param scale Causal estimand scale.
#' @noRd
eq_tr11 <- function(data, models, scale){
  Z <- data$Z
  D <- data$D
  Y <- data$Y

  Xpi <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$pi_model))
  Xp0 <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$p0_model))
  Xp1 <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$p1_model))
  Xm01 <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$m01_model))
  Xm11 <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$m11_model))

  pi_pos <- 1:ncol(Xpi)
  p0_pos <- (max(pi_pos) + 1):(max(pi_pos) + ncol(Xp0))
  p1_pos <- (max(p0_pos) + 1):(max(p0_pos) + ncol(Xp1))
  m01_pos <- (max(p1_pos) + 1):(max(p1_pos) + ncol(Xm01))
  m11_pos <- (max(m01_pos) + 1):(max(m01_pos) + ncol(Xm11))

  pi_scores <- grab_psiFUN(models$pi_model, data)
  p0_scores <- grab_psiFUN(models$p0_model, data)
  p1_scores <- grab_psiFUN(models$p1_model, data)
  m01_scores <- grab_psiFUN(models$m01_model, data)
  m11_scores <- grab_psiFUN(models$m11_model, data)

  function(theta){
    proscore <- plogis(Xpi %*% theta[pi_pos])
    p0 <- plogis(Xp0 %*% theta[p0_pos])
    p1 <- plogis(Xp1 %*% theta[p1_pos])
    om01 <- Xm01 %*% theta[m01_pos]
    om11 <- Xm11 %*% theta[m11_pos]
    psi_D0 <- (1-Z)/(1-proscore)*(D-p0)+p0
    psi_D1 <- Z/proscore*(D-p1)+p1
    psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*p0)+om01*p0
    psi_YD1 <- Z/proscore*(Y*D-om11*p1)+om11*p1
    if (scale == "RD") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m01_scores(theta[m01_pos])*I(Z==0)*I(D==1),
        m11_scores(theta[m11_pos])*I(Z==1)*I(D==1),
        p0/p1*psi_YD1+om11*(psi_D0-p0/p1*psi_D1)-theta[max(m11_pos)+1]*psi_D0,
        psi_YD0-theta[max(m11_pos)+2]*psi_D0,
        theta[max(m11_pos)+1]-theta[max(m11_pos)+2]-theta[max(m11_pos)+3])
    } else if (scale == "RR") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m01_scores(theta[m01_pos])*I(Z==0)*I(D==1),
        m11_scores(theta[m11_pos])*I(Z==1)*I(D==1),
        p0/p1*psi_YD1+om11*(psi_D0-p0/p1*psi_D1)-theta[max(m11_pos)+1]*psi_D0,
        psi_YD0-theta[max(m11_pos)+2]*psi_D0,
        log(theta[max(m11_pos)+1])-log(theta[max(m11_pos)+2])-theta[max(m11_pos)+3])
    } else if (scale == "OR") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m01_scores(theta[m01_pos])*I(Z==0)*I(D==1),
        m11_scores(theta[m11_pos])*I(Z==1)*I(D==1),
        p0/p1*psi_YD1+om11*(psi_D0-p0/p1*psi_D1)-theta[max(m11_pos)+1]*psi_D0,
        psi_YD0-theta[max(m11_pos)+2]*psi_D0,
        log(theta[max(m11_pos)+1]/(1-theta[max(m11_pos)+1]))-log(theta[max(m11_pos)+2]/(1-theta[max(m11_pos)+2]))-theta[max(m11_pos)+3])
    }
  }
}


######################## Estimating Equations for CACE ###########################

#' Parametric conditionally doubly robust estimator for CACE
#' @param Z Treatment vector.
#' @param D Intermediate outcome vector.
#' @param Y Final outcome vector.
#' @param proscore Propensity scores.
#' @param prinscore Principal scores.
#' @param om Outcome mean predictions.
#' @param theta Odds ratio parameter.
#' @noRd
pqrc_point01 <- function(Z,D,Y,proscore,prinscore,om,theta){
  prinscore1 <- prinscore[[1]]
  prinscore0 <- prinscore[[2]]
  om11 <- om[[1]]
  om01 <- om[[2]]
  om10 <- om[[3]]
  om00 <- om[[4]]

  #Estimate Psi functions
  psi_D1 <- Z/proscore*(D-prinscore1)+prinscore1
  psi_D0 <- (1-Z)/(1-proscore)*(D-prinscore0)+prinscore0
  psi_YD1 <- Z/proscore*(Y*D-om11*prinscore1)+om11*prinscore1
  psi_Y1_D0 <- (1-Z)/(1-proscore)*(Y*(1-D)-om00*(1-prinscore0))+om00*(1-prinscore0)

  #tau and omega
  ## independence
  ex.ind <- prinscore1*(1-prinscore0)
  tau.ind <- psi_D1*(1-prinscore0)+(prinscore0-psi_D0)*prinscore1
  omega1.ind <- ex.ind/prinscore1*(psi_YD1-om11*psi_D1)+tau.ind*om11
  omega0.ind <- ex.ind/(1-prinscore0)*(psi_Y1_D0-om00*(1-psi_D0))+tau.ind*om00
  ## non-independence
  delta <- (1+(theta-1)*(prinscore1+prinscore0))^2-4*theta*(theta-1)*prinscore1*prinscore0
  lambda <- 1+(theta-1)*(prinscore0-prinscore1)-sqrt(delta)
  ex.non.ind = lambda/2/(1-theta)
  tau.non.ind <- ex.non.ind - 1/2/sqrt(delta)*((psi_D0-prinscore0)*(2*prinscore1-lambda)+(psi_D1-prinscore1)*(2*((theta-1)*(prinscore0-prinscore1)+prinscore0-sqrt(delta))-lambda))
  omega1.non.ind <- ex.non.ind/prinscore1*(psi_YD1-om11*psi_D1)+tau.non.ind*om11
  omega0.non.ind <- ex.non.ind/(1-prinscore0)*(psi_Y1_D0-om00*(1-psi_D0))+tau.non.ind*om00
  ## combine
  tau = ifelse(theta == 1, tau.ind, tau.non.ind)
  omega1 = ifelse(theta == 1, omega1.ind, omega1.non.ind)
  omega0 = ifelse(theta == 1, omega0.ind, omega0.non.ind)

  #parametric conditionally doubly robust estimator
  pqr1 <- mean(omega1)/mean(tau)
  pqr0 <- mean(omega0)/mean(tau)
  pqrc <- pqr1-pqr0
  denom_est <- mean(tau)
  num_est.1 <- mean(omega1)
  num_est.0 <- mean(omega0)

  return(c(pqr1,pqr0,pqrc,num_est.1,num_est.0,denom_est))
}

#' Variance of the DML estimator for CACE
#' @param Z Treatment vector.
#' @param D Intermediate outcome vector.
#' @param Y Final outcome vector.
#' @param proscore Propensity scores.
#' @param prinscore Principal scores.
#' @param om Outcome mean predictions.
#' @param theta Odds ratio parameter.
#' @noRd
ml_qr_var01 <- function(Z,D,Y,proscore,prinscore,om,theta){
  prinscore1 <- prinscore[[1]]
  prinscore0 <- prinscore[[2]]
  om11 <- om[[1]]
  om01 <- om[[2]]
  om10 <- om[[3]]
  om00 <- om[[4]]

  #Estimate Psi functions
  psi_D1 <- Z/proscore*(D-prinscore1)+prinscore1
  psi_D0 <- (1-Z)/(1-proscore)*(D-prinscore0)+prinscore0
  psi_YD1 <- Z/proscore*(Y*D-om11*prinscore1)+om11*prinscore1
  psi_Y1_D0 <- (1-Z)/(1-proscore)*(Y*(1-D)-om00*(1-prinscore0))+om00*(1-prinscore0)

  #tau and omega
  ## independence
  ex.ind <- prinscore1*(1-prinscore0)
  tau.ind <- psi_D1*(1-prinscore0)+(prinscore0-psi_D0)*prinscore1
  omega1.ind <- ex.ind/prinscore1*(psi_YD1-om11*psi_D1)+tau.ind*om11
  omega0.ind <- ex.ind/(1-prinscore0)*(psi_Y1_D0-om00*(1-psi_D0))+tau.ind*om00
  ## non-independence
  delta <- (1+(theta-1)*(prinscore1+prinscore0))^2-4*theta*(theta-1)*prinscore1*prinscore0
  lambda <- 1+(theta-1)*(prinscore0-prinscore1)-sqrt(delta)
  ex.non.ind = lambda/2/(1-theta)
  tau.non.ind <- ex.non.ind - 1/2/sqrt(delta)*((psi_D0-prinscore0)*(2*prinscore1-lambda)+(psi_D1-prinscore1)*(2*((theta-1)*(prinscore0-prinscore1)+prinscore0-sqrt(delta))-lambda))
  omega1.non.ind <- ex.non.ind/prinscore1*(psi_YD1-om11*psi_D1)+tau.non.ind*om11
  omega0.non.ind <- ex.non.ind/(1-prinscore0)*(psi_Y1_D0-om00*(1-psi_D0))+tau.non.ind*om00
  ## combine
  tau = ifelse(theta == 1, tau.ind, tau.non.ind)
  omega1 = ifelse(theta == 1, omega1.ind, omega1.non.ind)
  omega0 = ifelse(theta == 1, omega0.ind, omega0.non.ind)


  #Estimate the denominator
  denom_est <- mean(tau)

  #Compute the variance estimates
  pqr1 <- mean(omega1)/mean(tau)
  pqr0 <- mean(omega0)/mean(tau)
  xi1 <- omega1-pqr1*tau
  xi0 <- omega0-pqr0*tau
  var_single <- mean((xi1-xi0)^2)/denom_est^2
  var.0 = mean(xi0^2)/denom_est^2
  var.1 = mean(xi1^2)/denom_est^2
  cov = mean(xi0*xi1)/denom_est^2

  return(c(var_single,var.0,var.1,cov))
}

#' Parametric triply robust estimator for CACE under monotonicity
#' @param Z Treatment vector.
#' @param D Intermediate outcome vector.
#' @param Y Final outcome vector.
#' @param proscore Propensity scores.
#' @param prinscore Principal scores.
#' @param om Outcome mean predictions.
#' @noRd
ptrc_point01 <- function(Z,D,Y,proscore,prinscore,om){
  prinscore1 <- prinscore[[1]]
  prinscore0 <- prinscore[[2]]
  om11 <- om[[1]]
  om01 <- om[[2]]
  om10 <- om[[3]]
  om00 <- om[[4]]

  #Estimate Psi functions
  psi_D1 <- Z/proscore*(D-prinscore1)+prinscore1
  psi_D0 <- (1-Z)/(1-proscore)*(D-prinscore0)+prinscore0
  psi_YD1 <- Z/proscore*(Y*D-om11*prinscore1)+om11*prinscore1
  psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*prinscore0)+om01*prinscore0
  psi_Y1_D0 <- (1-Z)/(1-proscore)*(Y*(1-D)-om00*(1-prinscore0))+om00*(1-prinscore0)

  #parametric triply robust estimator
  ptr1 <- mean((prinscore1-prinscore0)/prinscore1*psi_YD1-om11*(psi_D0-prinscore0/prinscore1*psi_D1))/mean(psi_D1-psi_D0)
  ptr0 <- mean((prinscore1-prinscore0)/(1-prinscore0)*psi_Y1_D0-om00*(1-psi_D1-(1-prinscore1)/(1-prinscore0)*(1-psi_D0)))/mean(psi_D1-psi_D0)
  ptrc <- ptr1-ptr0
  denom_est <- mean(psi_D1-psi_D0)
  num_est.1 = mean((prinscore1-prinscore0)/prinscore1*psi_YD1-om11*(psi_D0-prinscore0/prinscore1*psi_D1))
  num_est.0 = mean((prinscore1-prinscore0)/(1-prinscore0)*psi_Y1_D0-om00*(1-psi_D1-(1-prinscore1)/(1-prinscore0)*(1-psi_D0)))

  return(c(ptr1,ptr0,ptrc,num_est.1,num_est.0,denom_est))
}

#' Variance of the DML estimator for CACE under monotonicity
#' @param Z Treatment vector.
#' @param D Intermediate outcome vector.
#' @param Y Final outcome vector.
#' @param proscore Propensity scores.
#' @param prinscore Principal scores.
#' @param om Outcome mean predictions.
#' @noRd
ml_tr_var01 <- function(Z,D,Y,proscore,prinscore,om){
  prinscore1 <- prinscore[[1]]
  prinscore0 <- prinscore[[2]]
  om11 <- om[[1]]
  om01 <- om[[2]]
  om10 <- om[[3]]
  om00 <- om[[4]]

  #Estimate Psi functions
  psi_D1 <- Z/proscore*(D-prinscore1)+prinscore1
  psi_D0 <- (1-Z)/(1-proscore)*(D-prinscore0)+prinscore0
  psi_YD1 <- Z/proscore*(Y*D-om11*prinscore1)+om11*prinscore1
  psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*prinscore0)+om01*prinscore0
  psi_Y1_D0 <- (1-Z)/(1-proscore)*(Y*(1-D)-om00*(1-prinscore0))+om00*(1-prinscore0)

  #Estimate the denominator
  denom_est <- mean(psi_D1-psi_D0)

  #Compute variance estimates
  ptr1 <- mean((prinscore1-prinscore0)/prinscore1*psi_YD1-om11*(psi_D0-prinscore0/prinscore1*psi_D1))/mean(psi_D1-psi_D0)
  ptr0 <- mean((prinscore1-prinscore0)/(1-prinscore0)*psi_Y1_D0-om00*(1-psi_D1-(1-prinscore1)/(1-prinscore0)*(1-psi_D0)))/mean(psi_D1-psi_D0)
  var_single <- mean((((prinscore1-prinscore0)/prinscore1*psi_YD1-om11*(psi_D0-prinscore0/prinscore1*psi_D1)-ptr1*(psi_D1-psi_D0))-((prinscore1-prinscore0)/(1-prinscore0)*psi_Y1_D0-om00*(1-psi_D1-(1-prinscore1)/(1-prinscore0)*(1-psi_D0))-ptr0*(psi_D1-psi_D0)))^2)/denom_est^2
  var.0 = mean(((prinscore1-prinscore0)/(1-prinscore0)*psi_Y1_D0-om00*(1-psi_D1-(1-prinscore1)/(1-prinscore0)*(1-psi_D0))-ptr0*(psi_D1-psi_D0))^2)/denom_est^2
  var.1 = mean(((prinscore1-prinscore0)/prinscore1*psi_YD1-om11*(psi_D0-prinscore0/prinscore1*psi_D1)-ptr1*(psi_D1-psi_D0))^2)/denom_est^2
  cov = mean(((prinscore1-prinscore0)/prinscore1*psi_YD1-om11*(psi_D0-prinscore0/prinscore1*psi_D1)-ptr1*(psi_D1-psi_D0))*((prinscore1-prinscore0)/(1-prinscore0)*psi_Y1_D0-om00*(1-psi_D1-(1-prinscore1)/(1-prinscore0)*(1-psi_D0))-ptr0*(psi_D1-psi_D0)))/denom_est^2

  return(c(var_single,var.0,var.1,cov))
}

#' Estimating equations for parametric multiply robust estimator for CACE
#' @param data Data frame.
#' @param models List of fitted models.
#' @param scale Causal estimand scale.
#' @noRd
eq_qr01 <- function(data, models, scale){
  Z <- data$Z
  D <- data$D
  Y <- data$Y
  oddsratio = data$or.fit

  Xpi <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$pi_model))
  Xp0 <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$p0_model))
  Xp1 <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$p1_model))
  Xm00 <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$m00_model))
  Xm11 <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$m11_model))

  pi_pos <- 1:ncol(Xpi)
  p0_pos <- (max(pi_pos) + 1):(max(pi_pos) + ncol(Xp0))
  p1_pos <- (max(p0_pos) + 1):(max(p0_pos) + ncol(Xp1))
  m00_pos <- (max(p1_pos) + 1):(max(p1_pos) + ncol(Xm00))
  m11_pos <- (max(m00_pos) + 1):(max(m00_pos) + ncol(Xm11))

  pi_scores <- grab_psiFUN(models$pi_model, data)
  p0_scores <- grab_psiFUN(models$p0_model, data)
  p1_scores <- grab_psiFUN(models$p1_model, data)
  m00_scores <- grab_psiFUN(models$m00_model, data)
  m11_scores <- grab_psiFUN(models$m11_model, data)

  function(theta){
    proscore <- plogis(Xpi %*% theta[pi_pos])
    p0 <- plogis(Xp0 %*% theta[p0_pos])
    p1 <- plogis(Xp1 %*% theta[p1_pos])
    om00 <- Xm00 %*% theta[m00_pos]
    om11 <- Xm11 %*% theta[m11_pos]
    psi_D0 <- (1-Z)/(1-proscore)*(D-p0)+p0
    psi_D1 <- Z/proscore*(D-p1)+p1
    psi_YD1 <- Z/proscore*(Y*D-om11*p1)+om11*p1
    psi_Y1_D0 <- (1-Z)/(1-proscore)*(Y*(1-D)-om00*(1-p0))+om00*(1-p0)
    ## independence
    ex.ind = p1*(1-p0)
    tau.ind <- psi_D1*(1-p0)+(p0-psi_D0)*p1
    omega1.ind <- ex.ind/p1*(psi_YD1-om11*psi_D1)+tau.ind*om11
    omega0.ind <- ex.ind/(1-p0)*(psi_Y1_D0-om00*(1-psi_D0))+tau.ind*om00
    ## non-independence
    delta <- (1+(oddsratio-1)*(p1+p0))^2-4*oddsratio*(oddsratio-1)*p1*p0
    lambda <- 1+(oddsratio-1)*(p0-p1)-sqrt(delta)
    ex.non.ind = -lambda/2/(oddsratio-1)
    tau.non.ind <- ex.non.ind - 1/2/sqrt(delta)*((psi_D0-p0)*(2*p1-lambda)+(psi_D1-p1)*(2*((oddsratio-1)*(p0-p1)+p0-sqrt(delta))-lambda))
    omega1.non.ind <- ex.non.ind/p1*(psi_YD1-om11*psi_D1)+tau.non.ind*om11
    omega0.non.ind <- ex.non.ind/(1-p0)*(psi_Y1_D0-om00*(1-psi_D0))+tau.non.ind*om00
    ## combine
    tau = ifelse(oddsratio == 1, tau.ind, tau.non.ind)
    omega1 = ifelse(oddsratio == 1, omega1.ind, omega1.non.ind)
    omega0 = ifelse(oddsratio == 1, omega0.ind, omega0.non.ind)
    if (scale == "RD") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m00_scores(theta[m00_pos])*I(Z==0)*I(D==0),
        m11_scores(theta[m11_pos])*I(Z==1)*I(D==1),
        omega1-theta[max(m11_pos)+1]*tau,
        omega0-theta[max(m11_pos)+2]*tau,
        theta[max(m11_pos)+1]-theta[max(m11_pos)+2]-theta[max(m11_pos)+3])
    } else if (scale == "RR") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m00_scores(theta[m00_pos])*I(Z==0)*I(D==0),
        m11_scores(theta[m11_pos])*I(Z==1)*I(D==1),
        omega1-theta[max(m11_pos)+1]*tau,
        omega0-theta[max(m11_pos)+2]*tau,
        log(theta[max(m11_pos)+1])-log(theta[max(m11_pos)+2])-theta[max(m11_pos)+3])
    } else if (scale == "OR") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m00_scores(theta[m00_pos])*I(Z==0)*I(D==0),
        m11_scores(theta[m11_pos])*I(Z==1)*I(D==1),
        omega1-theta[max(m11_pos)+1]*tau,
        omega0-theta[max(m11_pos)+2]*tau,
        log.odds(theta[max(m11_pos)+1])-log.odds(theta[max(m11_pos)+2])-theta[max(m11_pos)+3])
    }
  }
}

#' Estimating equations for parametric triply robust estimator for CACE under monotonicity
#' @param data Data frame.
#' @param models List of fitted models.
#' @param scale Causal estimand scale.
#' @noRd
eq_tr01 <- function(data, models, scale){
  Z <- data$Z
  D <- data$D
  Y <- data$Y

  Xpi <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$pi_model))
  Xp0 <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$p0_model))
  Xp1 <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$p1_model))
  Xm00 <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$m00_model))
  Xm11 <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$m11_model))

  pi_pos <- 1:ncol(Xpi)
  p0_pos <- (max(pi_pos) + 1):(max(pi_pos) + ncol(Xp0))
  p1_pos <- (max(p0_pos) + 1):(max(p0_pos) + ncol(Xp1))
  m00_pos <- (max(p1_pos) + 1):(max(p1_pos) + ncol(Xm00))
  m11_pos <- (max(m00_pos) + 1):(max(m00_pos) + ncol(Xm11))

  pi_scores <- grab_psiFUN(models$pi_model, data)
  p0_scores <- grab_psiFUN(models$p0_model, data)
  p1_scores <- grab_psiFUN(models$p1_model, data)
  m00_scores <- grab_psiFUN(models$m00_model, data)
  m11_scores <- grab_psiFUN(models$m11_model, data)

  function(theta){
    proscore <- plogis(Xpi %*% theta[pi_pos])
    p0 <- plogis(Xp0 %*% theta[p0_pos])
    p1 <- plogis(Xp1 %*% theta[p1_pos])
    om00 <- Xm00 %*% theta[m00_pos]
    om11 <- Xm11 %*% theta[m11_pos]
    psi_D0 <- (1-Z)/(1-proscore)*(D-p0)+p0
    psi_D1 <- Z/proscore*(D-p1)+p1
    psi_YD1 <- Z/proscore*(Y*D-om11*p1)+om11*p1
    psi_Y1_D0 <- (1-Z)/(1-proscore)*(Y*(1-D)-om00*(1-p0))+om00*(1-p0)
    if (scale == "RD") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m00_scores(theta[m00_pos])*I(Z==0)*I(D==0),
        m11_scores(theta[m11_pos])*I(Z==1)*I(D==1),
        (p1-p0)/p1*psi_YD1-om11*(psi_D0-p0/p1*psi_D1)-theta[max(m11_pos)+1]*(psi_D1-psi_D0),
        (p1-p0)/(1-p0)*psi_Y1_D0-om00*(1-psi_D1-(1-p1)/(1-p0)*(1-psi_D0))-theta[max(m11_pos)+2]*(psi_D1-psi_D0),
        theta[max(m11_pos)+1]-theta[max(m11_pos)+2]-theta[max(m11_pos)+3])
    } else if (scale == "RR") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m00_scores(theta[m00_pos])*I(Z==0)*I(D==0),
        m11_scores(theta[m11_pos])*I(Z==1)*I(D==1),
        (p1-p0)/p1*psi_YD1-om11*(psi_D0-p0/p1*psi_D1)-theta[max(m11_pos)+1]*(psi_D1-psi_D0),
        (p1-p0)/(1-p0)*psi_Y1_D0-om00*(1-psi_D1-(1-p1)/(1-p0)*(1-psi_D0))-theta[max(m11_pos)+2]*(psi_D1-psi_D0),
        log(theta[max(m11_pos)+1])-log(theta[max(m11_pos)+2])-theta[max(m11_pos)+3])
    } else if (scale == "OR") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m00_scores(theta[m00_pos])*I(Z==0)*I(D==0),
        m11_scores(theta[m11_pos])*I(Z==1)*I(D==1),
        (p1-p0)/p1*psi_YD1-om11*(psi_D0-p0/p1*psi_D1)-theta[max(m11_pos)+1]*(psi_D1-psi_D0),
        (p1-p0)/(1-p0)*psi_Y1_D0-om00*(1-psi_D1-(1-p1)/(1-p0)*(1-psi_D0))-theta[max(m11_pos)+2]*(psi_D1-psi_D0),
        log.odds(theta[max(m11_pos)+1])-log.odds(theta[max(m11_pos)+2])-theta[max(m11_pos)+3])
    }
  }
}

########################Estimating Equations for Never-takers###########################

#' Parametric conditionally doubly robust estimator for Never-Takers
#' @param Z Treatment vector.
#' @param D Intermediate outcome vector.
#' @param Y Final outcome vector.
#' @param proscore Propensity scores.
#' @param prinscore Principal scores.
#' @param om Outcome mean predictions.
#' @param theta Odds ratio parameter.
#' @noRd
pqrc_point00 <- function(Z,D,Y,proscore,prinscore,om,theta){
  prinscore1 <- prinscore[[1]]
  prinscore0 <- prinscore[[2]]
  om11 <- om[[1]]
  om01 <- om[[2]]
  om10 <- om[[3]]
  om00 <- om[[4]]

  #Estimate Psi functions
  psi_D1 <- Z/proscore*(D-prinscore1)+prinscore1
  psi_D0 <- (1-Z)/(1-proscore)*(D-prinscore0)+prinscore0
  psi_YD1 <- Z/proscore*(Y*D-om11*prinscore1)+om11*prinscore1
  psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*prinscore0)+om01*prinscore0
  psi_Y1_D1 <- Z/proscore*(Y*(1-D)-om10*(1-prinscore1))+om10*(1-prinscore1)
  psi_Y1_D0 <- (1-Z)/(1-proscore)*(Y*(1-D)-om00*(1-prinscore0))+om00*(1-prinscore0)

  #tau and omega
  ## independence
  ex.ind = (1-prinscore0)*(1-prinscore1)
  tau.ind <- (1-psi_D1)*(1-prinscore0)+(prinscore0-psi_D0)*(1-prinscore1)
  omega1.ind <- ex.ind/(1-prinscore1)*(psi_Y1_D1-om10*(1-psi_D1))+tau.ind*om10
  omega0.ind <- ex.ind/(1-prinscore0)*(psi_Y1_D0-om00*(1-psi_D0))+tau.ind*om00
  ## non-independence
  delta <- (1+(theta-1)*(prinscore1+prinscore0))^2-4*theta*(theta-1)*prinscore1*prinscore0
  lambda <- 1+(theta-1)*(2-prinscore0-prinscore1)-sqrt(delta)
  ex.non.ind = lambda/2/(theta-1)
  tau.non.ind <- ex.non.ind + 1/2/sqrt(delta)*((psi_D0-prinscore0)*(2*(prinscore1+(theta-1)*(1-prinscore0)-sqrt(delta))-lambda)+(psi_D1-prinscore1)*(2*(prinscore0+(theta-1)*(1-prinscore1)-sqrt(delta))-lambda))
  omega1.non.ind <- ex.non.ind/(1-prinscore1)*(psi_Y1_D1-om10*(1-psi_D1))+tau.non.ind*om10
  omega0.non.ind <- ex.non.ind/(1-prinscore0)*(psi_Y1_D0-om00*(1-psi_D0))+tau.non.ind*om00
  ## combine
  tau = ifelse(theta == 1, tau.ind, tau.non.ind)
  omega1 = ifelse(theta == 1, omega1.ind, omega1.non.ind)
  omega0 = ifelse(theta == 1, omega0.ind, omega0.non.ind)

  #parametric conditionally doubly robust estimator
  pqr1 <- mean(omega1)/mean(tau)
  pqr0 <- mean(omega0)/mean(tau)
  pqrc <- pqr1-pqr0
  denom_est <- mean(tau)
  num_est.1 <- mean(omega1)
  num_est.0 <- mean(omega0)

  return(c(pqr1,pqr0,pqrc,num_est.1,num_est.0,denom_est))
}

#' Variance of the DML estimator for Never-Takers
#' @param Z Treatment vector.
#' @param D Intermediate outcome vector.
#' @param Y Final outcome vector.
#' @param proscore Propensity scores.
#' @param prinscore Principal scores.
#' @param om Outcome mean predictions.
#' @param theta Odds ratio parameter.
#' @noRd
ml_qr_var00 <- function(Z,D,Y,proscore,prinscore,om,theta){
  prinscore1 <- prinscore[[1]]
  prinscore0 <- prinscore[[2]]
  om11 <- om[[1]]
  om01 <- om[[2]]
  om10 <- om[[3]]
  om00 <- om[[4]]

  #Estimate Psi functions
  psi_D1 <- Z/proscore*(D-prinscore1)+prinscore1
  psi_D0 <- (1-Z)/(1-proscore)*(D-prinscore0)+prinscore0
  psi_YD1 <- Z/proscore*(Y*D-om11*prinscore1)+om11*prinscore1
  psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*prinscore0)+om01*prinscore0
  psi_Y1_D1 <- Z/proscore*(Y*(1-D)-om10*(1-prinscore1))+om10*(1-prinscore1)
  psi_Y1_D0 <- (1-Z)/(1-proscore)*(Y*(1-D)-om00*(1-prinscore0))+om00*(1-prinscore0)

  #tau and omega
  ## independence
  ex.ind = (1-prinscore0)*(1-prinscore1)
  tau.ind <- (1-psi_D1)*(1-prinscore0)+(prinscore0-psi_D0)*(1-prinscore1)
  omega1.ind <- ex.ind/(1-prinscore1)*(psi_Y1_D1-om10*(1-psi_D1))+tau.ind*om10
  omega0.ind <- ex.ind/(1-prinscore0)*(psi_Y1_D0-om00*(1-psi_D0))+tau.ind*om00
  ## non-independence
  delta <- (1+(theta-1)*(prinscore1+prinscore0))^2-4*theta*(theta-1)*prinscore1*prinscore0
  lambda <- 1+(theta-1)*(2-prinscore0-prinscore1)-sqrt(delta)
  ex.non.ind = lambda/2/(theta-1)
  tau.non.ind <- ex.non.ind + 1/2/sqrt(delta)*((psi_D0-prinscore0)*(2*(prinscore1+(theta-1)*(1-prinscore0)-sqrt(delta))-lambda)+(psi_D1-prinscore1)*(2*(prinscore0+(theta-1)*(1-prinscore1)-sqrt(delta))-lambda))
  omega1.non.ind <- ex.non.ind/(1-prinscore1)*(psi_Y1_D1-om10*(1-psi_D1))+tau.non.ind*om10
  omega0.non.ind <- ex.non.ind/(1-prinscore0)*(psi_Y1_D0-om00*(1-psi_D0))+tau.non.ind*om00
  ## combine
  tau = ifelse(theta == 1, tau.ind, tau.non.ind)
  omega1 = ifelse(theta == 1, omega1.ind, omega1.non.ind)
  omega0 = ifelse(theta == 1, omega0.ind, omega0.non.ind)


  #Estimate the denominator
  denom_est <- mean(tau)

  #Compute the variance estimates
  pqr1 <- mean(omega1)/mean(tau)
  pqr0 <- mean(omega0)/mean(tau)
  xi1 <- omega1-pqr1*tau
  xi0 <- omega0-pqr0*tau
  var_single <- mean((xi1-xi0)^2)/denom_est^2
  var.0 = mean(xi0^2)/denom_est^2
  var.1 = mean(xi1^2)/denom_est^2
  cov = mean(xi0*xi1)/denom_est^2

  return(c(var_single,var.0,var.1,cov))
}

#' Parametric triply robust estimator for Never-Takers under monotonicity
#' @param Z Treatment vector.
#' @param D Intermediate outcome vector.
#' @param Y Final outcome vector.
#' @param proscore Propensity scores.
#' @param prinscore Principal scores.
#' @param om Outcome mean predictions.
#' @noRd
ptrc_point00 <- function(Z,D,Y,proscore,prinscore,om){
  prinscore1 <- prinscore[[1]]
  prinscore0 <- prinscore[[2]]
  om11 <- om[[1]]
  om01 <- om[[2]]
  om10 <- om[[3]]
  om00 <- om[[4]]

  #Estimate Psi functions
  psi_D1 <- Z/proscore*(D-prinscore1)+prinscore1
  psi_D0 <- (1-Z)/(1-proscore)*(D-prinscore0)+prinscore0
  psi_YD1 <- Z/proscore*(Y*D-om11*prinscore1)+om11*prinscore1
  psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*prinscore0)+om01*prinscore0
  psi_Y1_D1 <- Z/proscore*(Y*(1-D)-om10*(1-prinscore1))+om10*(1-prinscore1)
  psi_Y1_D0 <- (1-Z)/(1-proscore)*(Y*(1-D)-om00*(1-prinscore0))+om00*(1-prinscore0)

  #parametric quadruply robust estimator
  ptr1 <- mean(psi_Y1_D1)/mean(1-psi_D1)
  ptr0 <- mean((1-prinscore1)/(1-prinscore0)*psi_Y1_D0+om00*(1-psi_D1-(1-prinscore1)/(1-prinscore0)*(1-psi_D0)))/mean(1-psi_D1)
  ptrc <- ptr1-ptr0
  denom_est <- mean(1-psi_D1)
  num_est.1 = mean(psi_Y1_D1)
  num_est.0 = mean((1-prinscore1)/(1-prinscore0)*psi_Y1_D0+om00*(1-psi_D1-(1-prinscore1)/(1-prinscore0)*(1-psi_D0)))

  return(c(ptr1,ptr0,ptrc,num_est.1,num_est.0,denom_est))
}

#' Variance of the DML estimator for Never-Takers under monotonicity
#' @param Z Treatment vector.
#' @param D Intermediate outcome vector.
#' @param Y Final outcome vector.
#' @param proscore Propensity scores.
#' @param prinscore Principal scores.
#' @param om Outcome mean predictions.
#' @noRd
ml_tr_var00 <- function(Z,D,Y,proscore,prinscore,om){
  prinscore1 <- prinscore[[1]]
  prinscore0 <- prinscore[[2]]
  om11 <- om[[1]]
  om01 <- om[[2]]
  om10 <- om[[3]]
  om00 <- om[[4]]

  #Estimate Psi functions
  psi_D1 <- Z/proscore*(D-prinscore1)+prinscore1
  psi_D0 <- (1-Z)/(1-proscore)*(D-prinscore0)+prinscore0
  psi_YD1 <- Z/proscore*(Y*D-om11*prinscore1)+om11*prinscore1
  psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*prinscore0)+om01*prinscore0
  psi_Y1_D1 <- Z/proscore*(Y*(1-D)-om10*(1-prinscore1))+om10*(1-prinscore1)
  psi_Y1_D0 <- (1-Z)/(1-proscore)*(Y*(1-D)-om00*(1-prinscore0))+om00*(1-prinscore0)

  #Estimate the denominator
  denom_est <- mean(1-psi_D1)

  #Compute variance estimates
  ptr1 <- mean(psi_Y1_D1)/mean(1-psi_D1)
  ptr0 <- mean((1-prinscore1)/(1-prinscore0)*psi_Y1_D0+om00*(1-psi_D1-(1-prinscore1)/(1-prinscore0)*(1-psi_D0)))/mean(1-psi_D1)
  var_single <- mean(((psi_Y1_D1-ptr1*(1-psi_D1))-((1-prinscore1)/(1-prinscore0)*psi_Y1_D0+om00*(1-psi_D1-(1-prinscore1)/(1-prinscore0)*(1-psi_D0))-ptr0*(1-psi_D1)))^2)/denom_est^2
  var.0 = mean(((1-prinscore1)/(1-prinscore0)*psi_Y1_D0+om00*(1-psi_D1-(1-prinscore1)/(1-prinscore0)*(1-psi_D0))-ptr0*(1-psi_D1))^2)/denom_est^2
  var.1 = mean((psi_Y1_D1-ptr1*(1-psi_D1))^2)/denom_est^2
  cov = mean((psi_Y1_D1-ptr1*(1-psi_D1))*((1-prinscore1)/(1-prinscore0)*psi_Y1_D0+om00*(1-psi_D1-(1-prinscore1)/(1-prinscore0)*(1-psi_D0))-ptr0*(1-psi_D1)))/denom_est^2

  return(c(var_single,var.0,var.1,cov))
}

#' Estimating equations for parametric multiply robust estimator for Never-Takers
#' @param data Data frame.
#' @param models List of fitted models.
#' @param scale Causal estimand scale.
#' @noRd
eq_qr00 <- function(data, models, scale){
  Z <- data$Z
  D <- data$D
  Y <- data$Y
  oddsratio = data$or.fit

  Xpi <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$pi_model))
  Xp0 <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$p0_model))
  Xp1 <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$p1_model))
  Xm00 <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$m00_model))
  Xm10 <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$m10_model))

  pi_pos <- 1:ncol(Xpi)
  p0_pos <- (max(pi_pos) + 1):(max(pi_pos) + ncol(Xp0))
  p1_pos <- (max(p0_pos) + 1):(max(p0_pos) + ncol(Xp1))
  m00_pos <- (max(p1_pos) + 1):(max(p1_pos) + ncol(Xm00))
  m10_pos <- (max(m00_pos) + 1):(max(m00_pos) + ncol(Xm10))

  pi_scores <- grab_psiFUN(models$pi_model, data)
  p0_scores <- grab_psiFUN(models$p0_model, data)
  p1_scores <- grab_psiFUN(models$p1_model, data)
  m00_scores <- grab_psiFUN(models$m00_model, data)
  m10_scores <- grab_psiFUN(models$m10_model, data)

  function(theta){
    proscore <- plogis(Xpi %*% theta[pi_pos])
    p0 <- plogis(Xp0 %*% theta[p0_pos])
    p1 <- plogis(Xp1 %*% theta[p1_pos])
    om00 <- Xm00 %*% theta[m00_pos]
    om10 <- Xm10 %*% theta[m10_pos]
    psi_D0 <- (1-Z)/(1-proscore)*(D-p0)+p0
    psi_D1 <- Z/proscore*(D-p1)+p1
    psi_Y1_D1 <- Z/proscore*(Y*(1-D)-om10*(1-p1))+om10*(1-p1)
    psi_Y1_D0 <- (1-Z)/(1-proscore)*(Y*(1-D)-om00*(1-p0))+om00*(1-p0)
    ## independence
    ex.ind = (1-p0)*(1-p1)
    tau.ind <- (1-psi_D1)*(1-p0)+(p0-psi_D0)*(1-p1)
    omega1.ind <- ex.ind/(1-p1)*(psi_Y1_D1-om10*(1-psi_D1))+tau.ind*om10
    omega0.ind <- ex.ind/(1-p0)*(psi_Y1_D0-om00*(1-psi_D0))+tau.ind*om00
    ## non-independence
    delta <- (1+(oddsratio-1)*(p1+p0))^2-4*oddsratio*(oddsratio-1)*p1*p0
    lambda <- 1+(oddsratio-1)*(2-p0-p1)-sqrt(delta)
    ex.non.ind = lambda/2/(oddsratio-1)
    tau.non.ind <- ex.non.ind + 1/2/sqrt(delta)*((psi_D0-p0)*(2*(p1+(oddsratio-1)*(1-p0)-sqrt(delta))-lambda)+(psi_D1-p1)*(2*(p0+(oddsratio-1)*(1-p1)-sqrt(delta))-lambda))
    omega1.non.ind <- ex.non.ind/(1-p1)*(psi_Y1_D1-om10*(1-psi_D1))+tau.non.ind*om10
    omega0.non.ind <- ex.non.ind/(1-p0)*(psi_Y1_D0-om00*(1-psi_D0))+tau.non.ind*om00
    ## combine
    tau = ifelse(oddsratio == 1, tau.ind, tau.non.ind)
    omega1 = ifelse(oddsratio == 1, omega1.ind, omega1.non.ind)
    omega0 = ifelse(oddsratio == 1, omega0.ind, omega0.non.ind)
    if (scale == "RD") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m00_scores(theta[m00_pos])*I(Z==0)*I(D==0),
        m10_scores(theta[m10_pos])*I(Z==1)*I(D==0),
        omega1-theta[max(m10_pos)+1]*tau,
        omega0-theta[max(m10_pos)+2]*tau,
        theta[max(m10_pos)+1]-theta[max(m10_pos)+2]-theta[max(m10_pos)+3])
    } else if (scale == "RR") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m00_scores(theta[m00_pos])*I(Z==0)*I(D==0),
        m10_scores(theta[m10_pos])*I(Z==1)*I(D==0),
        omega1-theta[max(m10_pos)+1]*tau,
        omega0-theta[max(m10_pos)+2]*tau,
        log(theta[max(m10_pos)+1])-log(theta[max(m10_pos)+2])-theta[max(m10_pos)+3])
    } else if (scale == "OR") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m00_scores(theta[m00_pos])*I(Z==0)*I(D==0),
        m10_scores(theta[m10_pos])*I(Z==1)*I(D==0),
        omega1-theta[max(m10_pos)+1]*tau,
        omega0-theta[max(m10_pos)+2]*tau,
        log.odds(theta[max(m10_pos)+1])-log.odds(theta[max(m10_pos)+2])-theta[max(m10_pos)+3])
    }
  }
}

#' Estimating equations for parametric triply robust estimator for Never-Takers under monotonicity
#' @param data Data frame.
#' @param models List of fitted models.
#' @param scale Causal estimand scale.
#' @noRd
eq_tr00 <- function(data, models, scale){
  Z <- data$Z
  D <- data$D
  Y <- data$Y

  Xpi <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$pi_model))
  Xp0 <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$p0_model))
  Xp1 <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$p1_model))
  Xm00 <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$m00_model))
  Xm10 <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$m10_model))

  pi_pos <- 1:ncol(Xpi)
  p0_pos <- (max(pi_pos) + 1):(max(pi_pos) + ncol(Xp0))
  p1_pos <- (max(p0_pos) + 1):(max(p0_pos) + ncol(Xp1))
  m00_pos <- (max(p1_pos) + 1):(max(p1_pos) + ncol(Xm00))
  m10_pos <- (max(m00_pos) + 1):(max(m00_pos) + ncol(Xm10))

  pi_scores <- grab_psiFUN(models$pi_model, data)
  p0_scores <- grab_psiFUN(models$p0_model, data)
  p1_scores <- grab_psiFUN(models$p1_model, data)
  m00_scores <- grab_psiFUN(models$m00_model, data)
  m10_scores <- grab_psiFUN(models$m10_model, data)

  function(theta){
    proscore <- plogis(Xpi %*% theta[pi_pos])
    p0 <- plogis(Xp0 %*% theta[p0_pos])
    p1 <- plogis(Xp1 %*% theta[p1_pos])
    om00 <- Xm00 %*% theta[m00_pos]
    om10 <- Xm10 %*% theta[m10_pos]
    psi_D0 <- (1-Z)/(1-proscore)*(D-p0)+p0
    psi_D1 <- Z/proscore*(D-p1)+p1
    psi_Y1_D1 <- Z/proscore*(Y*(1-D)-om10*(1-p1))+om10*(1-p1)
    psi_Y1_D0 <- (1-Z)/(1-proscore)*(Y*(1-D)-om00*(1-p0))+om00*(1-p0)
    if (scale == "RD") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m00_scores(theta[m00_pos])*I(Z==0)*I(D==0),
        m10_scores(theta[m10_pos])*I(Z==1)*I(D==0),
        psi_Y1_D1-theta[max(m10_pos)+1]*(1-psi_D1),
        (1-p1)/(1-p0)*psi_Y1_D0+om00*(1-psi_D1-(1-p1)/(1-p0)*(1-psi_D0))-theta[max(m10_pos)+2]*(1-psi_D1),
        theta[max(m10_pos)+1]-theta[max(m10_pos)+2]-theta[max(m10_pos)+3])
    } else if (scale == "RR") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m00_scores(theta[m00_pos])*I(Z==0)*I(D==0),
        m10_scores(theta[m10_pos])*I(Z==1)*I(D==0),
        psi_Y1_D1-theta[max(m10_pos)+1]*(1-psi_D1),
        (1-p1)/(1-p0)*psi_Y1_D0+om00*(1-psi_D1-(1-p1)/(1-p0)*(1-psi_D0))-theta[max(m10_pos)+2]*(1-psi_D1),
        log(theta[max(m10_pos)+1])-log(theta[max(m10_pos)+2])-theta[max(m10_pos)+3])
    } else if (scale == "OR") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m00_scores(theta[m00_pos])*I(Z==0)*I(D==0),
        m10_scores(theta[m10_pos])*I(Z==1)*I(D==0),
        psi_Y1_D1-theta[max(m10_pos)+1]*(1-psi_D1),
        (1-p1)/(1-p0)*psi_Y1_D0+om00*(1-psi_D1-(1-p1)/(1-p0)*(1-psi_D0))-theta[max(m10_pos)+2]*(1-psi_D1),
        log.odds(theta[max(m10_pos)+1])-log.odds(theta[max(m10_pos)+2])-theta[max(m10_pos)+3])
    }
  }
}

########################Estimating Equations for defiers###########################

#' Parametric conditionally doubly robust estimator for Defiers
#' @param Z Treatment vector.
#' @param D Intermediate outcome vector.
#' @param Y Final outcome vector.
#' @param proscore Propensity scores.
#' @param prinscore Principal scores.
#' @param om Outcome mean predictions.
#' @param theta Odds ratio parameter.
#' @noRd
pqrc_point10 <- function(Z,D,Y,proscore,prinscore,om,theta){
  prinscore1 <- prinscore[[1]]
  prinscore0 <- prinscore[[2]]
  om11 <- om[[1]]
  om01 <- om[[2]]
  om10 <- om[[3]]
  om00 <- om[[4]]

  #Estimate Psi functions
  psi_D1 <- Z/proscore*(D-prinscore1)+prinscore1
  psi_D0 <- (1-Z)/(1-proscore)*(D-prinscore0)+prinscore0
  psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*prinscore0)+om01*prinscore0
  psi_Y1_D1 <- Z/proscore*(Y*(1-D)-om10*(1-prinscore1))+om10*(1-prinscore1)

  #tau and omega
  ## independence
  ex.ind = prinscore0*(1-prinscore1)
  tau.ind <- prinscore0*(1-psi_D1)+(psi_D0-prinscore0)*(1-prinscore1)
  omega1.ind <- ex.ind/(1-prinscore1)*(psi_Y1_D1-om10*(1-psi_D1))+tau.ind*om10
  omega0.ind <- ex.ind/prinscore0*(psi_YD0-om01*psi_D0)+tau.ind*om01
  ## non-independence
  delta <- (1+(theta-1)*(prinscore1+prinscore0))^2-4*theta*(theta-1)*prinscore1*prinscore0
  lambda <- 1+(theta-1)*(prinscore1-prinscore0)-sqrt(delta)
  ex.non.ind = -lambda/2/(theta-1)
  tau.non.ind <- ex.non.ind - 1/2/sqrt(delta)*((psi_D0-prinscore0)*(2*((theta-1)*(prinscore1-prinscore0)+prinscore1-sqrt(delta))-lambda)+(psi_D1-prinscore1)*(2*prinscore0-lambda))
  omega1.non.ind <- ex.non.ind/(1-prinscore1)*(psi_Y1_D1-om10*(1-psi_D1))+tau.non.ind*om10
  omega0.non.ind <- ex.non.ind/prinscore0*(psi_YD0-om01*psi_D0)+tau.non.ind*om01
  ## combine
  tau = ifelse(theta == 1, tau.ind, tau.non.ind)
  omega1 = ifelse(theta == 1, omega1.ind, omega1.non.ind)
  omega0 = ifelse(theta == 1, omega0.ind, omega0.non.ind)

  #parametric conditionally doubly robust estimator
  pqr1 <- mean(omega1)/mean(tau)
  pqr0 <- mean(omega0)/mean(tau)
  pqrc <- pqr1-pqr0
  denom_est <- mean(tau)
  num_est.1 <- mean(omega1)
  num_est.0 <- mean(omega0)

  return(c(pqr1,pqr0,pqrc,num_est.1,num_est.0,denom_est))
}

#' Variance of the DML estimator for Defiers
#' @param Z Treatment vector.
#' @param D Intermediate outcome vector.
#' @param Y Final outcome vector.
#' @param proscore Propensity scores.
#' @param prinscore Principal scores.
#' @param om Outcome mean predictions.
#' @param theta Odds ratio parameter.
#' @noRd
ml_qr_var10 <- function(Z,D,Y,proscore,prinscore,om,theta){
  prinscore1 <- prinscore[[1]]
  prinscore0 <- prinscore[[2]]
  om11 <- om[[1]]
  om01 <- om[[2]]
  om10 <- om[[3]]
  om00 <- om[[4]]

  #Estimate Psi functions
  psi_D1 <- Z/proscore*(D-prinscore1)+prinscore1
  psi_D0 <- (1-Z)/(1-proscore)*(D-prinscore0)+prinscore0
  psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*prinscore0)+om01*prinscore0
  psi_Y1_D1 <- Z/proscore*(Y*(1-D)-om10*(1-prinscore1))+om10*(1-prinscore1)

  #tau and omega
  ## independence
  ex.ind = prinscore0*(1-prinscore1)
  tau.ind <- prinscore0*(1-psi_D1)+(psi_D0-prinscore0)*(1-prinscore1)
  omega1.ind <- ex.ind/(1-prinscore1)*(psi_Y1_D1-om10*(1-psi_D1))+tau.ind*om10
  omega0.ind <- ex.ind/prinscore0*(psi_YD0-om01*psi_D0)+tau.ind*om01
  ## non-independence
  delta <- (1+(theta-1)*(prinscore1+prinscore0))^2-4*theta*(theta-1)*prinscore1*prinscore0
  lambda <- 1+(theta-1)*(prinscore1-prinscore0)-sqrt(delta)
  ex.non.ind = -lambda/2/(theta-1)
  tau.non.ind <- ex.non.ind - 1/2/sqrt(delta)*((psi_D0-prinscore0)*(2*((theta-1)*(prinscore1-prinscore0)+prinscore1-sqrt(delta))-lambda)+(psi_D1-prinscore1)*(2*prinscore0-lambda))
  omega1.non.ind <- ex.non.ind/(1-prinscore1)*(psi_Y1_D1-om10*(1-psi_D1))+tau.non.ind*om10
  omega0.non.ind <- ex.non.ind/prinscore0*(psi_YD0-om01*psi_D0)+tau.non.ind*om01
  ## combine
  tau = ifelse(theta == 1, tau.ind, tau.non.ind)
  omega1 = ifelse(theta == 1, omega1.ind, omega1.non.ind)
  omega0 = ifelse(theta == 1, omega0.ind, omega0.non.ind)


  #Estimate the denominator
  denom_est <- mean(tau)

  #Compute the variance estimates
  pqr1 <- mean(omega1)/mean(tau)
  pqr0 <- mean(omega0)/mean(tau)
  xi1 <- omega1-pqr1*tau
  xi0 <- omega0-pqr0*tau
  var_single <- mean((xi1-xi0)^2)/denom_est^2
  var.0 = mean(xi0^2)/denom_est^2
  var.1 = mean(xi1^2)/denom_est^2
  cov = mean(xi0*xi1)/denom_est^2

  return(c(var_single,var.0,var.1,cov))
}

#' Estimating equations for parametric multiply robust estimator for Defiers
#' @param data Data frame.
#' @param models List of fitted models.
#' @param scale Causal estimand scale.
#' @noRd
eq_qr10 <- function(data, models, scale){
  Z <- data$Z
  D <- data$D
  Y <- data$Y
  oddsratio = data$or.fit

  Xpi <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$pi_model))
  Xp0 <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$p0_model))
  Xp1 <- grab_design_matrix(data = data,
                            rhs_formula = grab_fixed_formula(models$p1_model))
  Xm01 <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$m01_model))
  Xm10 <- grab_design_matrix(data = data,
                             rhs_formula = grab_fixed_formula(models$m10_model))

  pi_pos <- 1:ncol(Xpi)
  p0_pos <- (max(pi_pos) + 1):(max(pi_pos) + ncol(Xp0))
  p1_pos <- (max(p0_pos) + 1):(max(p0_pos) + ncol(Xp1))
  m01_pos <- (max(p1_pos) + 1):(max(p1_pos) + ncol(Xm01))
  m10_pos <- (max(m01_pos) + 1):(max(m01_pos) + ncol(Xm10))

  pi_scores <- grab_psiFUN(models$pi_model, data)
  p0_scores <- grab_psiFUN(models$p0_model, data)
  p1_scores <- grab_psiFUN(models$p1_model, data)
  m01_scores <- grab_psiFUN(models$m01_model, data)
  m10_scores <- grab_psiFUN(models$m10_model, data)

  function(theta){
    proscore <- plogis(Xpi %*% theta[pi_pos])
    p0 <- plogis(Xp0 %*% theta[p0_pos])
    p1 <- plogis(Xp1 %*% theta[p1_pos])
    om01 <- Xm01 %*% theta[m01_pos]
    om10 <- Xm10 %*% theta[m10_pos]
    psi_D0 <- (1-Z)/(1-proscore)*(D-p0)+p0
    psi_D1 <- Z/proscore*(D-p1)+p1
    psi_YD0 <- (1-Z)/(1-proscore)*(Y*D-om01*p0)+om01*p0
    psi_Y1_D1 <- Z/proscore*(Y*(1-D)-om10*(1-p1))+om10*(1-p1)
    ## independence
    ex.ind = p0*(1-p1)
    tau.ind <- p0*(1-psi_D1)+(psi_D0-p0)*(1-p1)
    omega1.ind <- ex.ind/(1-p1)*(psi_Y1_D1-om10*(1-psi_D1))+tau.ind*om10
    omega0.ind <- ex.ind/p0*(psi_YD0-om01*psi_D0)+tau.ind*om01
    ## non-independence
    delta <- (1+(oddsratio-1)*(p1+p0))^2-4*oddsratio*(oddsratio-1)*p1*p0
    lambda <- 1+(oddsratio-1)*(p1-p0)-sqrt(delta)
    ex.non.ind = -lambda/2/(oddsratio-1)
    tau.non.ind <- ex.non.ind - 1/2/sqrt(delta)*((psi_D0-p0)*(2*((oddsratio-1)*(p1-p0)+p1-sqrt(delta))-lambda)+(psi_D1-p1)*(2*p0-lambda))
    omega1.non.ind <- ex.non.ind/(1-p1)*(psi_Y1_D1-om10*(1-psi_D1))+tau.non.ind*om10
    omega0.non.ind <- ex.non.ind/p0*(psi_YD0-om01*psi_D0)+tau.non.ind*om01
    ## combine
    tau = ifelse(oddsratio == 1, tau.ind, tau.non.ind)
    omega1 = ifelse(oddsratio == 1, omega1.ind, omega1.non.ind)
    omega0 = ifelse(oddsratio == 1, omega0.ind, omega0.non.ind)
    if (scale == "RD") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m01_scores(theta[m01_pos])*I(Z==0)*I(D==1),
        m10_scores(theta[m10_pos])*I(Z==1)*I(D==0),
        omega1-theta[max(m10_pos)+1]*tau,
        omega0-theta[max(m10_pos)+2]*tau,
        theta[max(m10_pos)+1]-theta[max(m10_pos)+2]-theta[max(m10_pos)+3])
    } else if (scale == "RR") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m01_scores(theta[m01_pos])*I(Z==0)*I(D==1),
        m10_scores(theta[m10_pos])*I(Z==1)*I(D==0),
        omega1-theta[max(m10_pos)+1]*tau,
        omega0-theta[max(m10_pos)+2]*tau,
        log(theta[max(m10_pos)+1])-log(theta[max(m10_pos)+2])-theta[max(m10_pos)+3])
    } else if (scale == "OR") {
      c(pi_scores(theta[pi_pos]),
        p0_scores(theta[p0_pos])*I(Z==0),
        p1_scores(theta[p1_pos])*I(Z==1),
        m01_scores(theta[m01_pos])*I(Z==0)*I(D==1),
        m10_scores(theta[m10_pos])*I(Z==1)*I(D==0),
        omega1-theta[max(m10_pos)+1]*tau,
        omega0-theta[max(m10_pos)+2]*tau,
        log.odds(theta[max(m10_pos)+1])-log.odds(theta[max(m10_pos)+2])-theta[max(m10_pos)+3])
    }
  }
}


######################## Miscellaneous utility functions ###########################

#' Custom jacobian function for geex
#' @param func A function to differentiate.
#' @param x A vector of values at which to differentiate.
#' @param ... Additional arguments passed to numDeriv::jacobian.
#' @noRd
custom_jacobian <- function(func, x, ...) {
  numDeriv::jacobian(func, x, ...)
}

#' Generate stratified folds for cross-fitting
#' @param df Data frame.
#' @param n_group Number of folds.
#' @noRd
gene_Fi <- function(df,n_group){
  # Note: The original code had a bug returning df$G which doesn't exist.
  # Corrected to return the group vector.
  df$strata <- interaction(df$Z, df$D)
  folds <- createFolds(df$strata, k = n_group, list = TRUE, returnTrain = FALSE)
  df$Group <- NA
  for (i in seq_along(folds)) {
    df$Group[folds[[i]]] <- i
  }
  return(df$Group)
}

#' Inverse logit (expit) function
#' @param x A numeric vector.
#' @noRd
expit <- function(x){return(exp(x)/(1+exp(x)))}

#' Log-odds function
#' @param x A numeric vector of probabilities.
#' @noRd
log.odds <- function(x){return(log(x/(1-x)))}
