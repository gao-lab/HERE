report.expr <- function(expr){
    expr.substituted <- substitute(expr)
    cat(gsub(pattern="[ \t]{2,}", replacement=" ", x=paste0(collapse="", capture.output(print(expr.substituted)))), "\n") ## remove newlines and replace all successive whitespaces with a single whitespace
    eval(expr.substituted)
}
