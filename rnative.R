
dyn.load("rnative.so")

basicFunction <- function(x, y){  .Call("basicFunction", x, y) }

getVersion <- function(){ .Call("GetVersion")}

metrics <- function(xs){ .Call("metrics", xs) }

# Compute: a * xs + b * ys + c
linearcomb <- function(xs, ys, a, b, c){
    .Call("linearcomb", xs, ys, a, b, c)
}


testPackage <- function(){
    print(sprintf("basicFunction(4, 5) = %.3f", basicFunction(4, 5)))
    print(sprintf("basicFunction(10, 6) = %.3f", basicFunction(10, 6)))

    metrics(8)
    metrics(c(3.4, 9.0, 10.0, -15.75, 20.65))
}

testPackage()

# x = .Call("computeStatistics", c(3.0, 6.0, 8.0, -5.0, 10.0))
# x

