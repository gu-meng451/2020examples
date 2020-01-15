# Julia

Download it at [julialang.org](git@github.com:gu-ensc481/2020examples.git)

This is a really great open-source project that aims to make technical computing fast for development (like Python or Matlab), and fast for execution (like C).  In this folder there are (or will be) examples on getting started in Julia 1.x.  Many of them will be in the form of a Jupyter Notebook.

## Documentation
All the docs of the core of Julia are post [here](https://docs.julialang.org/).
There are also many tutorials [here](https://julialang.org/learning/).
If you're just getting started, I suggest sitting down with some popcorn and working along side this [video](https://www.youtube.com/watch?v=8h8rQyEpiZA&t).
Another great resource is the stuff posted by [MIT](https://github.com/mitmath/julia-mit).

## The Package manager
Julia has a built in way to load packages (libraries, functions, etc) from a wide variety of places.  There is an official registry that users can submit their packages to, and we can access them through `Pkg`.  The project's documentation is [here](https://julialang.github.io/Pkg.jl/v1/).  Of particular importance is how to add a package.  The first time is done with this
```julia
julia> using Pkg
julia> Pkg.add("LinearAlgebra")
```
and then to use the package, all that's needed is
```julia
julia> using LinearAlgebra
```

### Environments
It's a good idea to not install too many things in the root.  Instead you can make an Environment for each project or directory you are working in and keep the loaded packages local to there.  The documentation on that is [here](https://julialang.github.io/Pkg.jl/v1/environments/).

### Packages to think about using
There are many useful packages, and some have functions that (if coming from Matlab) we may be surprised aren't in the core.  The following two are examples of this
- `Printf` [docs](https://docs.julialang.org/en/v1.1/stdlib/Printf/), has the Matlab equivalent of writing formatted text like `fprintf`.
- `LinearAlgebra` [docs](https://docs.julialang.org/en/v1.1/stdlib/LinearAlgebra/), has functions like trace, determinant, and inverse.

Other packages are useful for specific tasks
- For symbolic stuff: `SymPy`
- For plotting start with `Plots`, I tend to use `PyPlot`
- For working with polynomials there's `Polynomials`
- For simple root finding like Matlab's `fzero` there's the package `Roots`


## Jupyter and Julia
**(This is not the only way to do this, but I recommend it as the first way.)** Launch the Julia REPL and type
```julia
julia> using Pkg
julia> Pkg.add("IJulia")
```
This will go out on the Internet and download the package and its requirements.  As well as precompile (setup for use) the package.  After these steps, which you only need to do one, Jupyter is installed inside your Julia install and will show Julia as one of the kernel options.

To launch Jupyter from inside Julia type the followingThen to launch Jupyter we need to type
```julia
julia> using IJulia
julia> notebook()
```
