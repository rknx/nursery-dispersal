# Set up environment -----------------------------------------------------------

## Libraries
library(lme4)


# Model preparation ------------------------------------------------------------

## For null model
modelNull = function(.model, nullVar = NULL) {
    if (is.null(nullVar)) {
         formulaEls = sapply(
             findbars(formula(.model)),
             function(x) paste0("(", deparse(x), ")")
         )
    } else {
         formulaEls = paste0("(", deparse(str2lang(paste("1 |", nullVar))), ")")
    }
    newFormula = reformulate(formulaEls, response = formula(.model)[[2]])
    update(.model, newFormula)
}

## For main effect only model
modelMain = function(.model) {
    formulaEls = c(
        all.vars(nobars(formula(.model)[[3]])),
        sapply(
            findbars(formula(.model)[[3]]),
            function(x) paste0("(", deparse(x), ")")
        )
    )
    newFormula = reformulate(formulaEls, response = formula(.model)[[2]])
    update(.model, newFormula)
}

## For final model with new dataset
reModel = function(.model) {
    glmer(
        formula(.model),
        na.omit(model.frame(.model)),
        family(.model),
        glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    )
}


# Main wrapper function --------------------------------------------------------

modelSelection = function(.model, nullVar = NULL) {
    final = reModel(.model) # !important: anova doesn't play nice with NAs.
    null = modelNull(final, nullVar)
    main = modelMain(final)
    aov = anova(null, main, final)
    stats = setNames(aov,
        c(
            "Number of parameters", "AIC", "BIC", "Log likelihood",
            "Deviance", "LR Chi-sq", "DF (LR Chi-sq)", "Pr(> LR Chi-sq)"
        )
    )
    out = format(as.data.frame(t(round(stats, 4))), scientific = F, digits = 6)
    setNames(out, c("Null Model", "Main Effects Model", "Final Model"))
}
