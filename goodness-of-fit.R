# Goodness of fit and model accuracy stats

## Association measures
gofCalc = function(targ, pred, thres = 0.5, ...) {
    if ((n = length(targ)) != length(pred)) stop("targ and pred length vary.")

    pos = pred[targ == 1]
    neg = pred[targ == 0]
    pair = length(pos) * length(neg)

    conc = do.call(sum, lapply(pos, function(x) sum(x > neg)))
    disc = do.call(sum, lapply(pos, function(x) sum(x < neg)))

    data.frame(check.names = FALSE,
        "Concordance" = conc / pair,
        "Discordance" = disc / pair,
        "Tied" = 1 - (conc + disc) / pair,
        "Number of Pairs (mil)" = pair / 1000000,
        "Somers D" = (conc - disc) / pair,
        "Goodman - Kruskal's Gamma" = (conc - disc) / (conc + disc),
        "Kendall's Tau A" = 2 * (conc - disc) / (n * (n - 1)),
        "Kendall's Tau B" = cor(targ, pred, method = "kendall"),
        "Spearman's rho" = cor(targ, pred, method = "spearman")
    )
}

## Pseudo-r-squared values
rSquared = function(.model) {
    vF = var(model.matrix(.model) %*% fixef(.model))
    vT = vF + sum(as.numeric(VarCorr(.model)))
    data.frame(check.names = FALSE,
        "R-squared (marginal)" = vF / (vT + pi ^ 2 / 3),
        "R-squared (conditional)" = vT / (vT + pi ^ 2 / 3)    )
}

## Main wrapper function
gofTest = function(.model, ...) {
    if (! "glmerMod" %in% class(.model)) stop("Model should be class glmerMod.")
    table = data.frame()
    gof = gofCalc(targ = .model@resp$y, pred = fitted(.model), ...)
    rsq = rSquared(.model)
    table = cbind(gof, rsq, row.names = "Final model")
    format(as.data.frame(t(round(table, 4))), scientific = F, digits = 4)
}
