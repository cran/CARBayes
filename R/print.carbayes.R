print.carbayes <- function(x,...)
{
#### Print out the model fitted
cat("\n#################\n")
cat("#### Model fitted\n")
cat("#################\n")
cat(x$model)
cat("Regression equation - ")
print(x$formula)

#### Print out the results
cat("\n############\n")
cat("#### Results\n")
cat("############\n")
cat("Posterior quantiles and DIC\n\n")
print(x$summary.results)
cat("\nDIC = ", x$modelfit[1], "     ", "p.d = ", x$modelfit[2], "\n")
return(invisible(x))
}