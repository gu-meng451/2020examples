#' ---
#' title : Let's Play with Scoping in Julia
#' author : T Fitzgerald
#' date : 15 April 2020
#' options:
#'  doctype : md2pdf
#'  fig_ext : .pdf
#'  template : tex.tpl
#' ---

#' # Introduction
#'
#' Scoping in Julia can seem complicated, and while it follows Python that may
#' not help if don't know Python.  Certain structures like
#' while/for/if blocks all have locally scoped variables.  We can pass things
#' through this in several ways.  See the [documentation](https://docs.julialang.org/en/v1/manual/variables-and-scoping/)
#' for more info.

#' ## Using a function
#' The first way I'll use is a function.  Inside
#' a function, everything is in local scope

function func1( x, a::Real)
    iter = 0
    flag = 0
    while flag == 0
        iter += 1

        if iter >= a
            flag = 1
        end
    end

    # now we have counted (very sloppily) up to an integer around a
    # I'm going to return that integer + x
    return x + iter
end;

#+ term=true
func1(0.5, pi)

#' ## Using a let block
#' Now we'll use a let block and inside this block things are local. I could
#' place everything inside the let block, but I'll demonstate passing things in
#' and out of one as well as the dreaded `global` keyword

x = 0.5 # it's ok to define this out here since we are not changing it
iter_outside = 0 # this we will need to pass in since we are changing it
a = pi # this one is also not changing
let iter = iter_outside
    flag = 0
    while flag == 0
        iter += 1

        if iter >= a
            flag = 1
        end
    end

    # now we have counted (very sloppily) up to an integer around a
    # I'm going to return that integer + x in a new variable.  I need to use the
    # global keyword so that y exists after the let block.
    global y = x + iter
end

#+ term=true
y

#' I personally don't like this very much, but this is world we live in.
