# Set up environment -----------------------------------------------------------

## Libraries
library(ggplot2)
library(gridExtra)
library(lme4)
library(extrafont)
library(pbapply)

#Import system font segoe UI, only for first time
if (! "Open Sans" %in% fonts()) font_import(pattern="opensans", prompt=F)


# Cross validation models ------------------------------------------------------

## Final model with coraviates
final = function(.model, ...) {
    targ = .model@resp$y
    pred = fitted(.model)
    accu = modelAccuracy(targ, pred, ...)
    p = plotRoc(targ, pred, title = paste("Final model"), ...)
    list(accu, p)
}

## Leave-p-out internal cross validation
leavePOut = function(.model, p = 10, ...) {
    data = randomize(model.frame(.model), ...)
    ncv = split(rownames(data), seq_len(nrow(data)) %% ceiling(nrow(data) / p))
    df = trainModel(.model, data, ncv[1:10], "leave-p-out")
    accu = modelAccuracy(df$targ, df$pred)
    p = plotRoc(df$targ, df$pred, paste0("Leave-", p, "-out CV"))
    list(accu, p)
}

# k-fold internal cross validation
kFold = function(.model, k = 10, ...) {
    data = randomize(model.frame(.model), ...)
    ncv = split(rownames(data), seq_len(nrow(data)) %% k)
    df = trainModel(.model, data, ncv, "k-fold")
    accu = modelAccuracy(df$targ, df$pred)
    p = plotRoc(df$targ, df$pred, title = paste0(k, "-fold CV"))
    list(accu, p)
}

# Accuracy and ROC functions ---------------------------------------------------

## Goodness of fit and model accuracy stats
modelAccuracy = function(targ, pred, thres = 0.5, ...) {
    if ((n = length(targ)) != length(pred)) stop("targ and pred length vary.")
    pos = targ == 1
    neg = targ == 0
    a = sum(pred >= thres & pos)
    c = sum(pos) - a
    d = sum(pred < thres & neg)
    b = sum(neg) - d
    pExp = ((a + c) * (a + b) + (b + d) * (c + d)) / n^2
    mccDen = sqrt(as.numeric((a + b)) * (a + c) * (b + d) * (c + d))
    data.frame(
        "Sensitivity" = a / sum(pos),
        "Specificity" = d / sum(neg),
        "Correct class" = (a + d) / n,
        "Cohen's kappa" = ((a + d) / n - pExp) / (1 - pExp),
        "Matthew's Correlation Coefficient" = (a * d - b * c) / mccDen,
        "True Skill Statistic" = (a * d - b * c) / ((a + c) * (b + d)),
        check.names = FALSE
    )
}

## Receiver Operator Characteristics Curve
plotRoc = function(targ, pred, title = "", bin = 20, ...) {
    pos = targ == 1; neg = targ == 0
    roc = data.frame(th = seq(0, 1, length.out = bin + 1))
    roc$tpr = sapply(th, function(x) sum(pred >= x & pos) / sum(pos)),
    roc$fpr = sapply(th, function(x) sum(pred >= x & neg) / sum(neg)),
}

## Actual plotting function for ROC
plotROCFunc = function(roc, thres = 0.5, ...) {
    roc$hline = roc$tpr[roc$th == thres],
    roc$vline = roc$fpr[roc$th == thres],
    roc$auc = sum(
            sapply(2:nrow(roc),
                function(y) (roc$tpr[y] + roc$tpr[y - 1]) / 2
            ) * -1 * diff(roc$fpr)
        )
    )
    
    .plot = ggplot(roc, aes(x = fpr, y = tpr)) +
        geom_line(col = "#f55252", size = 1.5) +
        geom_point(size = 2, alpha = 0.75) +
        coord_fixed() +
        geom_line(aes(th, th),  col = "blue",  size = 1) +
        geom_hline(aes(yintercept = hline), linetype = 2) +
        geom_vline(aes(xintercept = vline), linetype = 2) +
        geom_text(aes(label = paste("AUC =", round(auc, 2))),
            size = 10, x = .65, y = .2) +
        labs(title = title, xlab = "FPR", ylab = "TPR") +
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 0),
            text = element_text(family = "Open Sans"),
            axis.title = element_text(size = 20),
            axis.text = element_text(size = 16),
            plot.margin = unit(c(1, 2, 1, 1), "lines")
        )
}



# Helper functions -------------------------------------------------------------

## Row randomization for cross validation
randomize = function(.data, seed = 10, ...) {
    set.seed(seed)
    .data = .data[sample(nrow(.data)), ]
    row.names(.data) = NULL
    .data
}

## Training for internal validation
trainModel = function(.model, .data, splitList, message = "training") {
    message(paste("Processing", message, ":", length(splitList), "iterations."))
    out = do.call(rbind,
        pblapply(splitList, function(x) {
            fit = glmer(formula(.model), .data[-as.numeric(x), ], family(.model))
            pred = predict(.., .data[x, ], type = "resp", allow.new.levels = T)
            data.frame(targ = .data[x, 1], pred = pred)
        })
    )
    out[complete.cases(out), ]
}



# Cross validation parent function ---------------------------------------------
xVal = function(.model, p = NULL, k = NULL, extdata = NULL, noCov = F, ...) {
    if (! "glmerMod" %in% class(.model)) stop("Model should be class glmerMod.")

    roc = list()
    table = data.frame()

    if (exists(".model")) {
      out = final(.model, ...)
      table = rbind(table, "Final model" = out[[1]])
      roc[["Final model"]] = out[[2]]
    }

    for (i in unique(p)) {
      out = leavePOut(.model, i, ...)
      table = rbind(table, "Leave-p-out CV" = out[[1]])
      roc[[paste0("Leave-", i, "-out CV")]] = out[[2]]
    }

    for (i in unique(k)) {
      out = kFold(.model, i, ...)
      table = rbind(table, "k-fold CV" = out[[1]])
      roc[[paste0(i, "-fold CV")]] = out[[2]]
    }

    table = as.data.frame(t(round(table, 4)))
    print(format(table, scientific = F, digits = 4))
    
    gridExtra::grid.arrange(grobs = roc, ncol = 3)
}
