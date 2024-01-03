# Overview

Here we provide functions required to do a computer search for minimal examples of graphs, of varying regularity, which contain cospectral yet non-similar pairs of vertices. The results of this search are reported alongside other results in our paper:

Title:
Cospectral vertices, walk-regular planar graphs and the echolocation problem

Authors:
Shi-Lei Kong
Emmett L. Wyman
Yakun Xi

Link:
...

# Dependencies

This code is written in the Julia language (https://julialang.org) and is intended to be executed in the Julia REPL.

This code uses the Graphs package, which can be found at (https://github.com/JuliaGraphs/Graphs.jl). Follow the link for installation instructions.

# How to use

Download the minimal_cospectral.jl file and open it in an IDE of your choice. The authors use Virtual Studio Code. Start the Julia REPL and run the file. You should now be able to use the functions.

# The basic functions

Here we describe the basic functions used below. First, we have two functions which produce exhaustive lists of non-isomorphic graphs.

* `get_graphs(n::Int, filter::Function)` produces an exhaustive list of non-isomorphic graphs on `n` vertices. The `filter` argument is a Boolean-valued function on the Graph type (e.g. `is_connected`), and is used to narrow our search.
* `get_reg_graphs(n::Int, d::Int, filter::Function)` works similarly to `get_graphs`, but returns only `d`-regular graphs on `n` vertices.

Next, we need a function which allows us to determine if two vertices are similar. To this end, we have:

* `get_orbits(g::Graph)` returns a list `orbits` of sets such that `orbit[i]` is the orbit of vertex `i` under the isomorphism group of the graph `g`.
* `is_vertex_transitive(g::Graph)` returns `true` if the graph `g` is vertex-transitive and `false` otherwise. This function calls `get_orbits`.

Next, we need functions which determine when pairs of vertices are cospectral, and of these which are also not similar.

* `get_loops(g::Graph)` returns an n by n-1 matrix `L` of integers, where n is the number of vertices in `g`. The entry `L[i,k]` is equal to `(A^k)_{i,i}`, i.e. the number of walks of length `k` starting and ending at vertex `i`. By Theorem 1.1 in our paper, vertices `i` and `j` are cospectral if and only if `L[i,:] == L[j,:]`.
* `is_walk_regular(g::Graphs)` returns `true` if `g` is walk-regular and `false` otherwise. This function calls `get_loops` and checks if all vertics are cospectral.
* `ncvs(g::Graph)` returns a list of all pairs `(i,j)` of vertices `i` and `j` which are cospectral but not similar. Calls both `get_loops` and `get_orbits`.
* `has_ncvs(g::Graph)` is a function which returns `true` if the graph `g` contains non-similar cospectral vertices. This is fed as a filter to the search functions above. Calls `nonsimilar_cospectral_vertices` to check if the result is empty.

# How to search for minimal examples

The theorems in the paper which we verify with this code are as follows. (All graphs are finite, simple, undirected, and loopless.)

> Theorem 1.2. The smallest graph which contains a pair of vertices which are cospectral but not similar has eight vertices.

> Theorem 1.3. The smallest regular graph containing a pair of non-similar cospectral vertices has ten vertices.

> Theorem 1.4. The smallest walk-regular graph containing a pair of non-similar vertices has twelve vertices.

## Proof of Theorem 1.2

First, we generate a list `graphlist` of lists of graphs, such that `graphlist[n]` for `n` from 1 to 7 is the list of all graphs on `n` vertices having non-similar cospectral vertices. We expect these lists to be empty.

```
graphlists = [get_graphs(n, has_ncvs) for n in 1:7]
```

This will likely take a second or two to execute, and the output reveals all lists in `graphlists` are empty. To verify there exist such graphs on 8 vertices, we can enter

```
length(get_graphs(8, has_ncvs))
```

to which the REPL should respond after a few seconds by printing `126`. That is, there are `126` such minimal examples.

## Proof of Theorem 1.3

First, we verify there are no d-regular graphs on n vertices which contain nonsimilar cospectral pairs, with d ranging from 0 to n-1 and n between 8 and 9. We have already checked the more general case for smaller n in Theorem 1.2. Run:

```
for n in 8:9, d in 0:n-1
    print("n, d = " * string(n) * ", " * string(d) * ": ")
    println(isempty(get_reg_graphs(n,d,has_ncvs)))
end
```

Next, we verify there exist regular graphs on 10 vertices with nonsimilar cospectral pairs by running:

```
reggraphlists = [get_reg_graphs(10,d,has_ncvs) for d in 1:8]
```

This will take a few seconds to run and will produce a list of graph lists, where `reggraphlist[d]` lists the d-regular graphs on 10 vertices which contain nonsimilar cospectral vertices. Note excluding the d = 0 case (empty graph) and the n-1 case (complete graph) lose us nothing. We write

```
map(length, reggraphlists)
```

which reveals there are 3 examples of degree 3, 22 examples of degree 4, 22 examples of degree 5, and 3 examples of degree 6, and none others. Furthermore, for the sake of applying Proposition 3.1, we will need to verify that each of these graphs are connected. For this we may simply write

```
connectedlists = [filter(is_connected, gs) for gs in reggraphlists]
map(length, connectedlists)
```

and note the output is the same.

## Proof of Theorem 1.4

We begin by defining a new filter in the REPL:

```
is_wrnvt(g) = is_walk_regular(g) && !is_vertex_transitive(g)
```

which returns `true` if the graph `g` is walk-regular and non-vertex-transitive and `false` otherwise. This will be our filter. We start by checking that there are no such graphs on 10 or 11 vertices by entering:

```
for n in 10:11, d in 0:n-1
    print("n, d = " * string(n) * ", " * string(d) * ": ")
    println(isempty(get_reg_graphs(n,d,is_wrnvt)))
end
```

(Recall that all walk-regular graphs are also regular.) Finally, we must find our minimal examples amongst graphs on 12 vertices. Similar to before, we write

```
wrgraphlists = [get_reg_graphs(12,d,is_wrnvt) for d in 1:10]
```

This takes a few minutes to execute. Running

```
map(length, wrgraphlists)
```

shows there is one walk-regular, non-vertex-transitive graph on 12 vertices each of degree 4, 5, 6, and 7, and no others.

# Checks for validity

Many of the functions are easy to check for validity. The hardest to check are likely the search functions `get_graphs` and `get_reg_graphs`, which use a recursive search with some optimizations. Looking into the code, we see indeed that we use the `has_isomorph` function from the Graphs package to ensure that no two graphs in the resulting list are isomorphic.

Thankfully, we don't need to prove the search is exhaustive for all cases, we only need to show it's exhaustive for the cases we consider.

To verify `get_graphs(n)` is exhaustive, we check that the length of the result is equal to the number of non-isomorphic graphs on `n` unlabeled nodes from the Online Encyclopedia of Integer Sequences (https://oeis.org/A000088) for `n` from 1 to 8. The function `test()` does just this. This function takes a minute or two to execute.

To verify `get_reg_graphs(n)` is exhaustive, we do very much the same thing, checking against the number of non-isomorphic regular graphs on `n` unlabeled nodes (https://oeis.org/A005176). The case `n == 12` takes at least multiple days to run on my personal computer, and I have not let it finish running. Thankfully, show this case is exhaustive isn't as important to the result as the cases `n < 12`.
