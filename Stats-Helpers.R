# function for is nan
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
  
lmer.std.beta <- function(mod) {
   b <- fixef(mod)[-1]
   sd.x <- apply(getME(mod,"X")[,-1],2,sd)
   sd.y <- sd(getME(mod,"y"))
   b*sd.x/sd.y
}

confint.rlmerMod <- function(object, parm, level = 0.95) {
  beta <- fixef(object)
  if (missing(parm))
    parm <- names(beta)
  se <- sqrt(diag(vcov(object)))
  z <- qnorm((1 + level) / 2)
  ctab <- cbind(beta - z * se, beta + z * se)
  colnames(ctab) <- stats:::format.perc(c((1 - level) / 2, (1 + level) /
                                            2),
                                        digits = 3)
  return(ctab[parm, ])
}

rlmerMod.to.glm <- function(mod) {
  family <- gaussian()
  link <- family$link
  family <- family$family
  cl <- mod@call
  cl$control <- glm.control(epsilon = 1)
  .m <- match(
    c(
      "formula",
      "family",
      "data",
      "weights",
      "subset",
      "na.action",
      "offset",
      "model",
      "contrasts"
    ),
    names(cl),
    0L
  )
  cl <- cl[c(1L, .m)]
  cl[[1L]] <- as.name("glm")
  cl$formula <- effects:::fixFormula(as.formula(cl$formula))
  mod2 <- eval(cl)
  mod2$coefficients <- lme4::fixef(mod)
  mod2$vcov <- as.matrix(vcov(mod))
  mod2$linear.predictors <- model.matrix(mod2) %*% mod2$coefficients
  mod2$fitted.values <- mod2$family$linkinv(mod2$linear.predictors)
  mod2$weights <-
    as.vector(with(mod2, prior.weights * (
      family$mu.eta(linear.predictors) ^ 2 / family$variance(fitted.values)
    )))
  mod2$residuals <-
    with(mod2, prior.weights * (y - fitted.values) / weights)
  class(mod2) <- c("fakeglm", class(mod2))
  mod2
}


overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp ^ 2)
  prat <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(
    chisq = Pearson.chisq,
    ratio = prat,
    rdf = rdf,
    p = pval
  )
}

diag.plot.lmer <- function(model) {
  fitted.vals <- fitted(model)
  resid.vals <- resid(model)
  data.f <- data.frame(fitted.vals, resid.vals)
  g1 <-
    ggplot(data.f, aes(x = fitted.vals, y = resid.vals)) + geom_point() + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) + xlab('Fitted values') + ylab('Residuals')
  g2 <-
    ggqqplot(resid(model)) + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) + xlab('Theoretical') + ylab('Residuals')
  g3 <-
    gghistogram(resid(model), bins = 10) + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) + xlab('Residuals') + ylab('Count')
  figure <-
    ggarrange(
      g1,
      ggarrange(
        g2,
        g3,
        labels = c("B", "C"),
        ncol = 2,
        nrow = 1
      ),
      labels = c("A"),
      ncol = 1,
      nrow = 2
    )
}

diag.plot.lme <- function(model) {
  fitted.vals <- fitted(model)
  resid.vals <- resid(model, type = "pearson")
  data.f <- data.frame(fitted.vals, resid.vals)
  g1 <-
    ggplot(data.f, aes(x = fitted.vals, y = resid.vals)) + geom_point() + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) + xlab('Fitted values') + ylab('Residuals')
  g2 <-
    ggqqplot(resid(model)) + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) + xlab('Theoretical') + ylab('Residuals')
  g3 <-
    gghistogram(resid(model), bins = 10) + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    ) + xlab('Residuals') + ylab('Count')
  figure <-
    ggarrange(
      g1,
      ggarrange(
        g2,
        g3,
        labels = c("B", "C"),
        ncol = 2,
        nrow = 1
      ),
      labels = c("A"),
      ncol = 1,
      nrow = 2
    )
}

get.table.emmeans <- function(emmeans.object) {
  d <- as.data.frame(emmeans.object)
  d$t <- (d$emmean - d$Baseline) / d$SE
  d$p <- 2 * pt(abs(d$t), d$df, lower.tail = FALSE)
  if (sum(d$df == Inf) > 0) {
    d$df <- NULL
    names(d) <-
      c(
        "Group",
        "Time",
        "Baseline",
        "Mean",
        "SE",
        "95% CI lower",
        "95% CI upper",
        "Z-value",
        "P-value"
      )
    d[, seq(3, 8)] <- round(d[, seq(3, 8)], 1)
    d[, 9] <- round(d[, 9], 5)
  } else {
    d <- d[, c(seq(1, 5), 7, 8, 6, 9, 10)]
    names(d) <-
      c(
        "Group",
        "Time",
        "Baseline",
        "Mean",
        "SE",
        "95% CI lower",
        "95% CI upper",
        "df",
        "T-value",
        "P-value"
      )
    d[, seq(3, 9)] <- round(d[, seq(3, 9)], 1)
    d[, 10] <- round(d[, 10], 5)
  }
  
  
  c <-
    as.data.frame(summary(
      pairs(emmeans.object),
      infer = c(T, T),
      adjust = "none"
    ))
  
  names(c) <-
    c(
      "Contrast",
      "Time",
      "Baseline",
      "Mean",
      "SE",
      "df",
      "95% CI lower",
      "95% CI upper",
      "T-value",
      "P-value"
    )
  
  c <- c[, c(1, 2, 3, 4, 5, 7, 8, 6, 9, 10)]
  
  c[, seq(3, 9)] <- round(c[, seq(3, 9)], 1)
  
  c[, 10] <- round(c[, 10], 5)
  
  c$Baseline <- NULL
  
  t <- list(d, c)
}

var.table.model <- function(model) {
  var <- as.data.frame(unlist(insight::get_variance(model)))
  var$Variance <-var$`unlist(insight::get_variance(model))`
  var$`unlist(insight::get_variance(model))` <- NULL
  var$Component <- row.names(var)
  row.names(var) <- NULL
  var <- var[, c(2, 1)]
  
  var.total <- var[var$Component == "var.fixed", 2] + var[var$Component == "var.random", 2] + var[var$Component == "var.residual", 2]
  var$Percentage <- var$Variance / var.total * 100
  
  var$Variance <- round(var$Variance, 4)
  var$Percentage <- round(var$Percentage, 1)
  
  return (var)
}

make.factor <- function(Data.vector, levels, labels) {
  Result.vector <- as.factor(Data.vector)
  Result.vector <- factor(Result.vector, levels = levels, labels = labels)
  Result.vector <- fct_explicit_na(Result.vector, na_level = "Missing")
}

present.factor <- function(Datasource, var.name) {
  gg <- ggplot(Datasource, aes(x = Datasource[, var.name])) + geom_bar() + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none") + xlab(var.name) + ylab('N')
  
  total <- length(Datasource[, var.name])
  stats <- as.data.frame(Datasource %>% group_by(.dots = var.name) %>% summarise(N = n(), Percentage = round(n()/total*100, 1)))
  
  return(list(graph = gg, table = stats))
}

present.continuous <- function(Datasource, var.name, na.rm = FALSE) {
  gHist <- gghistogram(Datasource[, var.name], bins = 7) + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none") + xlab(var.name) + ylab('N')
  gQQ <- ggqqplot(Datasource[, var.name]) + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none")
  
  gBP <- ggboxplot(Datasource[, var.name]) + coord_flip() + scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none") + xlab('') + ylab(var.name)
  
  gg <- ggarrange(ggarrange(gHist, gQQ, nrow = 1, ncol = 2), gBP, nrow = 2, ncol = 1)
  
  stats <- data.frame(list(Mean = mean(Datasource[, var.name], na.rm = na.rm), Median = median(Datasource[, var.name], na.rm = na.rm), Min = min(Datasource[, var.name], na.rm = na.rm), Max = max(Datasource[, var.name], na.rm = na.rm), SD = sd(Datasource[, var.name], na.rm = na.rm), IQR = IQR(Datasource[, var.name], na.rm = na.rm)))
  return(list(graph = gg, table = stats))
}

Univariate.summary.table <- function(Datasource, exclude.vars, na.rm = FALSE, round.to = 1) {
  stats.table <- data.frame(list(Variable = NA, Category = NA, N = NA, Percentage = NA, Mean = NA, Median = NA, Min = NA, Max = NA, SD = NA, IQR = NA))
  for (var.name in names(Datasource)) {
    
    if(sum(exclude.vars == var.name) > 0) {
      next
    }
    
    if(is.factor(Datasource[, var.name])) {
      total <- length(Datasource[, var.name])
      stats <- as.data.frame(Datasource %>% group_by(.dots = var.name) %>% summarise(N = n(), Percentage = round(n()/total*100, 1)))
      stats$Mean <- NA
      stats$Median <- NA
      stats$Min <- NA
      stats$Max <- NA
      stats$SD <- NA
      stats$IQR <- NA
      stats$Variable <- var.name
      stats <- stats[, c(10, seq(1, 9))]
      names(stats)[2] <- "Category"
      stats[seq(2, length(stats[, 1])), 1] <- NA
      stats.table <- rbind(stats.table, stats)
    } else {
      stats <- data.frame(list(N = NA, Percentage = NA, Mean = mean(Datasource[, var.name], na.rm = na.rm), Median = median(Datasource[, var.name], na.rm = na.rm), Min = min(Datasource[, var.name], na.rm = na.rm), Max = max(Datasource[, var.name], na.rm = na.rm), SD = sd(Datasource[, var.name], na.rm = na.rm), IQR = IQR(Datasource[, var.name], na.rm = na.rm)))
      stats$Variable <- var.name
      stats$Category <- NA
      stats <- stats[, c(9, 10, seq(1, 8))]
      stats.table <- rbind(stats.table, stats)
    }
  }
  
  stats.table[, seq(3, 10)] <- round(stats.table[, seq(3, 10)], round.to)
  stats.table <- stats.table[-c(1), ]
  rownames(stats.table) <- NULL
  
  return(stats.table)
}

Plot.factor.vs.factor <- function(Datasource, factor.name.x, factor.name.y) {
  gg<- ggplot(Datasource, aes(x = Datasource[, factor.name.x], y = Datasource[, factor.name.y])) + geom_jitter(width = 0.2, height = 0.2) +
    scale_color_nejm() + scale_fill_nejm() + theme_minimal() + theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none") + xlab(factor.name.x) + ylab(factor.name.y)
  
  return(gg)
}

standard.p.value <- function(value, sigD = 3) {
  text <- if(value < 10^(-sigD)) {paste0("<", 10^(-sigD))} else {round(value, sigD)}
  text
}


