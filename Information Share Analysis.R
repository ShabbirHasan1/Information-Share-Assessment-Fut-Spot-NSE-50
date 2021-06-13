library(readxl)
library(readr)
setwd("/Users/parthsarthy/Desktop/Research/Stock Data")
getwd()

library(vars)
library(urca)

MVARselect <- function (y, lag.max = 10, type = c("const", "trend", "both", 
                                                  "none"), season = NULL, exogen = NULL) 
{
  y <- as.matrix(y)
  if (any(is.na(y))) 
    stop("\nNAs in y.\n")
  colnames(y) <- make.names(colnames(y))
  K <- ncol(y)
  lag.max <- abs(as.integer(lag.max))
  type <- match.arg(type)
  lag <- abs(as.integer(lag.max + 1))
  ylagged <- embed(y, lag)[, -c(1:K)]
  yendog <- y[-c(1:lag.max), ]
  sample <- nrow(ylagged)
  rhs <- switch("const", const = rep(1, sample), trend = seq(lag.max + 
                                                               1, length = sample), both = cbind(rep(1, sample), seq(lag.max + 
                                                                                                                       1, length = sample)), none = NULL)
  if (!(is.null(season))) {
    season <- abs(as.integer(season))
    dum <- (diag(season) - 1/season)[, -season]
    dums <- dum
    while (nrow(dums) < sample) {
      dums <- rbind(dums, dum)
    }
    dums <- dums[1:sample, ]
    rhs <- cbind(rhs, dums)
  }
  if (!(is.null(exogen))) {
    exogen <- as.matrix(exogen)
    if (!identical(nrow(exogen), nrow(y))) {
      stop("\nDifferent row size of y and exogen.\n")
    }
    if (is.null(colnames(exogen))) {
      colnames(exogen) <- paste("exo", 1:ncol(exogen), 
                                sep = "")
      warning(paste("No column names supplied in exogen, using:", 
                    paste(colnames(exogen), collapse = ", "), ", instead.\n"))
    }
    colnames(exogen) <- make.names(colnames(exogen))
    rhs <- cbind(rhs, exogen[-c(1:lag.max), ])
  }
  idx <- seq(K, K * lag.max, K)
  if (!is.null(rhs)) {
    detint <- ncol(as.matrix(rhs))
  }
  else {
    detint <- 0
  }
  criteria <- matrix(NA, nrow = 4, ncol = lag.max)
  rownames(criteria) <- c("AIC(n)", "HQ(n)", "SC(n)", "FPE(n)")
  colnames(criteria) <- paste(seq(1:lag.max))
  for (i in 1:lag.max) {
    ys.lagged <- cbind(ylagged[, c(1:idx[i])], rhs)
    sampletot <- nrow(y)
    nstar <- ncol(ys.lagged)
    resids <- lm.fit(x=ys.lagged, y=yendog)$residuals
    sigma.det <- det(crossprod(resids)/sample)
    criteria[1, i] <- log(sigma.det) + (2/sample) * (i * 
                                                       K^2 + K * detint)
    criteria[2, i] <- log(sigma.det) + (2 * log(log(sample))/sample) * 
      (i * K^2 + K * detint)
    criteria[3, i] <- log(sigma.det) + (log(sample)/sample) * 
      (i * K^2 + K * detint)
    criteria[4, i] <- ((sample + nstar)/(sample - nstar))^K * 
      sigma.det
  }
  order <- apply(criteria, 1, which.min)
  return(list(selection = order, criteria = criteria))
}


pdshare <- function(x, override.lags = NULL, lag.max = 10) {
  stopifnot(ncol(x)==2)
  stopifnot(is.numeric(x[,1]))
  stopifnot(is.numeric(x[,2]))
  if (is.null(override.lags)){
    if(lag.max<2) stop("Minimum lags should be 2")
  } else {
    if(override.lags<2) stop("Minimum lags should be 2")
  }
  cnames <- colnames(x)
  pdshare.computation <- function(x, nlag) {
    cointest <- ca.jo(x, K = nlag, type = "eigen", ecdet = "const",
                      spec = "transitory")  
    k <- cointest@lag
    vecm <- cajorls(cointest)
    varm <- vec2var(cointest)
    vma <- Psi(varm)
    ## converts level VARS to VMA model and gives orthogonalised psi
    ## matrix. 
    
    ## Head towards IS
    ## We need Psi(1), the matrix that captures the long run impact 
    ## of a disturbance on each of the prices.
    ## Psi(1) = beta.orthogonal*
    ##        [inverse(transpose(alpha.orthogonal) *gamma*
    ## (beta.orthogonal))] * transpose(alpha.orthogonal)
    
    ## the beta_orthogonal and alpha_orthogonal vectors :
    beta.ort <- as.vector(c(-cointest@V[2,1], cointest@V[1,1]))
    alpha.ort <- as.vector(c(-cointest@W[2,1], cointest@W[1,1]))
    
    ## initializing the parameters of gamma matrix
    aa <- bb <- cc <- dd <- 0
    for (i in 1:(k-1)) {
      aa <- aa + vecm$rlm$coefficients[2*i,1]
      bb <- bb + vecm$rlm$coefficients[2*i+1,1]
      cc <- cc + vecm$rlm$coefficients[2*i,2]
      dd <- dd + vecm$rlm$coefficients[2*i+1,2]
    }
    gamma.1 <- matrix(c(1-aa, -bb, -cc, 1-dd), nrow = 2, ncol = 2, byrow
                      = TRUE) 
    
    b <- as.numeric(t(alpha.ort) %*% gamma.1 %*% beta.ort)
    psi <- (beta.ort %*% t(alpha.ort))/b
    
    ## Information share is: (psi[1,]* f)_j^2) /
    ##                       (psi[1,]*omega*transpose(psi[1,])) 
    ## where f is the cholesky factorization of the omega matrix. 
    
    f <- vma[,,1]
    omega <- f %*% t(f)
    psi <- t(psi[1,])
    n <- psi %*% f
    d <- psi %*% omega %*% t(psi)
    
    list(ishares = c((n[, 1]^2)/d, (n[, 2]^2)/d), alpha.ort = alpha.ort, 
         omega = omega, lags = varm$p)
  }
  
  # Choosing the number of lags
  if (is.null(override.lags)) {
    nlag <- MVARselect(x, lag.max=lag.max)$selection[1] 
  } else {
    nlag <- override.lags
  }
  
  # First do the supplied ordering
  tmp <- pdshare.computation(x, nlag)
  is.original.ordering <- as.data.frame(tmp$ishares)
  component.share <- as.data.frame(abs(tmp$alpha.ort)/sum(abs(tmp$alpha.ort)))
  var.covar.matrix <- tmp$omega
  lags.used <- tmp$lags
  
  # Do the reverse ordering
  tmp <- pdshare.computation(x[,c(2,1)], nlag)
  is.reversed.ordering <- as.data.frame(tmp$ishares)
  
  
  rownames(var.covar.matrix) <- colnames(var.covar.matrix) <-
    rownames(component.share) <- rownames(is.original.ordering) <- cnames
  rownames(is.reversed.ordering) <- c(cnames[2], cnames[1])
  colnames(is.original.ordering) <- colnames(is.reversed.ordering) <-
    "IS" 
  colnames(component.share) <- "CS"
  
  
  list(is.original.ordering = is.original.ordering,
       is.reversed.ordering = is.reversed.ordering,
       component.share = component.share,
       var.covar.matrix = var.covar.matrix,
       lags.used = lags.used)
}

library(dplyr)

all_stocks <- list.files(pattern = '.csv')
all_stocks

RIL <- read.csv(file = 'RIl.csv', header = TRUE)
RIL <- RIL[-1,]

ADSEZ <- read.csv(file = 'ADSEZ.csv', header = TRUE)
ADSEZ <- ADSEZ[-1,]

APNT <- read.csv(file = 'APNT.csv', header = TRUE)
APNT <- APNT[-1,]

AXSB <- read.csv(file = 'AXSB.csv', header = TRUE)
AXSB <- AXSB[-1,]

BAF <- read.csv(file = 'BAF.csv', header = TRUE)
BAF <- BAF[-1,]

BHARTI <- read.csv(file = 'BHARTI.csv', header = TRUE)
BHARTI <- BHARTI[-1,]

BHIN <- read.csv(file = 'BHIN.csv', header = TRUE)
BHIN <- BHIN[-1,]

BJAUT <- read.csv(file = 'BJAUT.csv', header = TRUE)
BJAUT <- BJAUT[-1,]

BJFIN <- read.csv(file = 'BJFIN.csv', header = TRUE)
BJFIN <- BJFIN[-1,]

BPCL <- read.csv(file = 'BPCL.csv', header = TRUE)
BPCL <- BPCL[-1,]

BRIT <- read.csv(file = 'BRIT.csv', header = TRUE)
BRIT <- BRIT[-1,]

CIPLA <- read.csv(file = 'CIPLA.csv', header = TRUE)
CIPLA <- CIPLA[-1,]

COAL <- read.csv(file = 'COAL.csv', header = TRUE)
COAL <- COAL[-1,]

DRRD <- read.csv(file = 'DRRD.csv', header = TRUE)
DRRD <- DRRD[-1,]

EIM <- read.csv(file = 'EIM.csv', header = TRUE)
EIM <- EIM[-1,]

GAIL <- read.csv(file = 'GAIL.csv', header = TRUE)
GAIL <- GAIL[-1,]

GRASIM <- read.csv(file = 'GRASIM.csv', header = TRUE)
GRASIM <- GRASIM[-1,]

HCLT <- read.csv(file = 'HCLT.csv', header = TRUE)
HCLT <- HCLT[-1,]

HDFCB <- read.csv(file = 'HDFCB.csv', header = TRUE)
HDFCB <- HDFCB[-1,]

HDFCLIFE <- read.csv(file = 'HDFCLIFE.csv', header = TRUE)
HDFCLIFE <- HDFCLIFE[-1,]

HMCL <- read.csv(file = 'HMCL.csv', header = TRUE)
HMCL <- HMCL[-1,]

HNDL <- read.csv(file = 'HNDL.csv', header = TRUE)
HNDL <- HNDL[-1,]

HUVR <- read.csv(file = 'HUVR.csv', header = TRUE)
HUVR <- HUVR[-1,]

ICICIBC <- read.csv(file = 'ICICIBC.csv', header = TRUE)
ICICIBC <- ICICIBC[-1,]

HNDL <- read.csv(file = 'HNDL.csv', header = TRUE)
HNDL <- HNDL[-1,]

IIB <- read.csv(file = 'IIB.csv', header = TRUE)
IIB <- IIB[-1,]

INFO <- read.csv(file = 'INFO.csv', header = TRUE)
INFO <- INFO[-1,]

IOCL <- read.csv(file = 'IOCL.csv', header = TRUE)
IOCL <- IOCL[-1,]

ITC <- read.csv(file = 'ITC.csv', header = TRUE)
ITC <- ITC[-1,]

JSTL <- read.csv(file = 'JSTL.csv', header = TRUE)
JSTL <- JSTL[-1,]

KMB <- read.csv(file = 'KMB.csv', header = TRUE)
KMB <- KMB[-1,]

LT <- read.csv(file = 'LT.csv', header = TRUE)
LT <- LT[-1,]

MM <- read.csv(file = 'MM.csv', header = TRUE)
MM <- MM[-1,]

MSIL <- read.csv(file = 'MSIL.csv', header = TRUE)
MSIL <- MSIL[-1,]

NEST <- read.csv(file = 'NEST.csv', header = TRUE)
NEST <- NEST[-1,]

NTPC <- read.csv(file = 'NTPC.csv', header = TRUE)
NTPC <- NTPC[-1,]

ONGC <- read.csv(file = 'ONGC.csv', header = TRUE)
ONGC <- ONGC[-1,]

PWGR <- read.csv(file = 'PWGR.csv', header = TRUE)
PWGR <- PWGR[-1,]

LT <- read.csv(file = 'LT.csv', header = TRUE)
LT <- LT[-1,]

SBIN <- read.csv(file = 'SBIN.csv', header = TRUE)
SBIN <- SBIN[-1,]

SRCM <- read.csv(file = 'SRCM.csv', header = TRUE)
SRCM <- SRCM[-1,]

SUNP <- read.csv(file = 'SUNP.csv', header = TRUE)
SUNP <- SUNP[-1,]

TATA <- read.csv(file = 'TATA.csv', header = TRUE)
TATA <- TATA[-1,]

TCS <- read.csv(file = 'TCS.csv', header = TRUE)
TCS <- TCS[-1,]

TECHM <- read.csv(file = 'TECHM.csv', header = TRUE)
TECHM <- TECHM[-1,]

TTAN <- read.csv(file = 'TTAN.csv', header = TRUE)
TTAN <- TTAN[-1,]

TTMT <- read.csv(file = 'TTMT.csv', header = TRUE)
TTMT <- TTMT[-1,]

UPLL <- read.csv(file = 'UPLL.csv', header = TRUE)
UPLL <- UPLL[-1,]

UTCEM <- read.csv(file = 'UTCEM.csv', header = TRUE)
UTCEM <- UTCEM[-1,]

WPRO <- read.csv(file = 'WPRO.csv', header = TRUE)
WPRO <- WPRO[-1,]

Z <- read.csv(file = 'Z.csv', header = TRUE)
Z <- Z[-1,]

x1 = Z %>%
  select(Stock_log,Futures_log)

Z_IS<- pdshare(x1,lag.max = 120)

x2 = RIL %>%
  select(Stock_log,Futures_log)

RIL_IS <- pdshare(x2, lag.max = 120)

pdshare(x2, lag.max = 120)

x3 = ADSEZ %>%
  select(Stock_log,Futures_log)

ADSEZ_IS <- pdshare(x3, lag.max = 120)

x4 = APNT %>%
  select(Stock_log,Futures_log)

APNT_IS <- pdshare(x4, lag.max = 120)

x5 = AXSB %>%
  select(Stock_log,Futures_log)

AXSB_IS <- pdshare(x5, lag.max = 120)

x6 = BAF %>%
  select(Stock_log,Futures_log)

BAF_IS <- pdshare(x6, lag.max = 120)

x7 = BHARTI %>%
  select(Stock_log,Futures_log)

BHARTI_IS <- pdshare(x7, lag.max = 120)

x8 = BHIN %>%
  select(Stock_log,Futures_log)

BHIN_IS <- pdshare(x8, lag.max = 120)

x9 = BJAUT %>%
  select(Stock_log,Futures_log)

BJAUT_IS <- pdshare(x9, lag.max = 120)

x10 = BJFIN %>%
  select(Stock_log,Futures_log)

BJFIN_IS <- pdshare(x10, lag.max = 120)

x11 = BPCL %>%
  select(Stock_log,Futures_log)

BPCL_IS <- pdshare(x11, lag.max = 120)

x12 = BRIT %>%
  select(Stock_log,Futures_log)

BRIT_IS <- pdshare(x12, lag.max = 120)

x13 = CIPLA %>%
  select(Stock_log,Futures_log)

CIPLA_IS <- pdshare(x13, lag.max = 120)

x14 = COAL %>%
  select(Stock_log,Futures_log)

COAL_IS <- pdshare(x14, lag.max = 120)

x15 = DRRD %>%
  select(Stock_log,Futures_log)

DRRD_IS <- pdshare(x15, lag.max = 120)

x16 = EIM %>%
  select(Stock_log,Futures_log)

EIM_IS <- pdshare(x16, lag.max = 120)

x17 = GAIL %>%
  select(Stock_log,Futures_log)

GAIL_IS <- pdshare(x17, lag.max = 120)

x18 = GRASIM %>%
  select(Stock_log,Futures_log)

GRASIM_IS <- pdshare(x18, lag.max = 120)

x19 = HCLT %>%
  select(Stock_log,Futures_log)

HCLT_IS <- pdshare(x19, lag.max = 120)

x20 = HDFCB %>%
  select(Stock_log,Futures_log)

HDFCB_IS <- pdshare(x20, lag.max = 120)

x21 = HDFCLIFE %>%
  select(Stock_log,Futures_log)

HDFCLIFE_IS <- pdshare(x21, lag.max = 120)

x22 = HMCL %>%
  select(Stock_log,Futures_log)

HMCL_IS <- pdshare(x22, lag.max = 120)

x23 = HNDL %>%
  select(Stock_log,Futures_log)

HNDL_IS <- pdshare(x23, lag.max = 120)

x24 = HUVR %>%
  select(Stock_log,Futures_log)

HUVR_IS <- pdshare(x24, lag.max = 120)

x25 = ICICIBC %>%
  select(Stock_log,Futures_log)

ICICIBC_IS <- pdshare(x25, lag.max = 120)

x26 = IIB %>%
  select(Stock_log,Futures_log)

IIB_IS <- pdshare(x26, lag.max = 120)

x27 = INFO %>%
  select(Stock_log,Futures_log)

INFO_IS <- pdshare(x27, lag.max = 120)

x28 = IOCL %>%
  select(Stock_log,Futures_log)

IOCL_IS <- pdshare(x28, lag.max = 120)

x29 = ITC %>%
  select(Stock_log,Futures_log)

ITC_IS <- pdshare(x29, lag.max = 120)

x30 = JSTL %>%
  select(Stock_log,Futures_log)

JSTL_IS <- pdshare(x30, lag.max = 120)

x31 = KMB %>%
  select(Stock_log,Futures_log)

KMB_IS <- pdshare(x31, lag.max = 120)

x32 = LT %>%
  select(Stock_log,Futures_log)

LT_IS <- pdshare(x32, lag.max = 120)

x33 = MM %>%
  select(Stock_log,Futures_log)

MM_IS <- pdshare(x33, lag.max = 120)

x34 = MSIL %>%
  select(Stock_log,Futures_log)

MSIL_IS <- pdshare(x34, lag.max = 120)

x35 = NEST %>%
  select(Stock_log,Futures_log)

NEST_IS <- pdshare(x35, lag.max = 120)

x36 = NTPC %>%
  select(Stock_log,Futures_log)

NTPC_IS <- pdshare(x36, lag.max = 120)

  
x37 = ONGC %>%
  select(Stock_log,Futures_log)

ONGC_IS <- pdshare(x37, lag.max = 120)

x38 = PWGR %>%
  select(Stock_log,Futures_log)

PWGR_IS <- pdshare(x38, lag.max = 120)

x39 = SBIN %>%
  select(Stock_log,Futures_log)

SBIN_IS <- pdshare(x39, lag.max = 120)

x40 = SRCM %>%
  select(Stock_log,Futures_log)

SRCM_IS <- pdshare(x40, lag.max = 120)

x41 = SUNP %>%
  select(Stock_log,Futures_log)

SUNP_IS <- pdshare(x41, lag.max = 120)

x42 = TATA %>%
  select(Stock_log,Futures_log)

TATA_IS <- pdshare(x42, lag.max = 120)

x43 = TCS %>%
  select(Stock_log,Futures_log)

TCS_IS <- pdshare(x43, lag.max = 120)

x44 = TECHM %>%
  select(Stock_log,Futures_log)

TECHM_IS <- pdshare(x44, lag.max = 120)

x45 = TTAN %>%
  select(Stock_log,Futures_log)

TTAN_IS <- pdshare(x45, lag.max = 120)

x46 = TTMT %>%
  select(Stock_log,Futures_log)

TTMT_IS <- pdshare(x46, lag.max = 120)

x47 = UPLL %>%
  select(Stock_log,Futures_log)

UPLL_IS <- pdshare(x47, lag.max = 120)

x48 = UTCEM %>%
  select(Stock_log,Futures_log)

UTCEM_IS <- pdshare(x48, lag.max = 120)

x49 = WPRO %>%
  select(Stock_log,Futures_log)

WPRO_IS <- pdshare(x49, lag.max = 120)

IS <- c(ADSEZ_IS, APNT_IS, AXSB_IS, BAF_IS, BHARTI_IS, BHIN_IS, BJAUT_IS, BJFIN_IS, BPCL_IS, BRIT_IS, CIPLA_IS, COAL_IS, DRRD_IS, EIM_IS, GAIL_IS, GRASIM_IS, HCLT_IS, HDFCB_IS, HDFCLIFE_IS, HMCL_IS, HNDL_IS, HUVR_IS, ICICIBC_IS, IIB_IS, INFO_IS, IOCL_IS, ITC_IS, JSTL_IS, KMB_IS, LT_IS, MM_IS, MSIL_IS, NEST_IS, NTPC_IS, ONGC_IS, PWGR_IS, RIL_IS, SBIN_IS, SRCM_IS, SUNP_IS, TATA_IS, TCS_IS, TECHM_IS, TTAN_IS, TTMT_IS, UPLL_IS, UTCEM_IS, WPRO_IS, Z_IS)

library(xlsx)
lapply(IS, function(x) write.table(data.frame(x), 'IS Results.csv'  , append= T, sep=',' ))


