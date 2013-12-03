arc <-
function (x, y, theta1, theta2, rad, border, ...) 
{
    args <- lapply(as.list(match.call())[-1L], eval, parent.frame())
    names <- if (is.null(names(args))) 
        character(length(args))
    else names(args)
    dovec <- names %in% vectorize.args
    do.call("mapply", c(FUN = FUN, args[dovec], MoreArgs = list(args[!dovec]), 
        SIMPLIFY = SIMPLIFY, USE.NAMES = USE.NAMES))
}
