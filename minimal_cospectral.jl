import Graphs: Graph, adjacency_matrix, nv, is_connected, is_tree
import Graphs.Experimental: has_isomorph, all_isomorph
using LinearAlgebra
import Base: open, write, close, read, readdir

"""
    get_graphs(d_seq::Vector{Int}, filter::Function=((G::Graph)->true)) -> list::Vector{Graph}

Returns an exhaustive list of simple, undirected, loop-less graphs with unlabeled
vertices of the given degree sequence `d_seq`. The list contains no more than one
graph from each isomorphism class. This function does not use `d_seq` directly,
but rather a sorted copy.

The optional `filter` function is used to exclude graphs from the list. This is
useful for speeding up searches for special kinds of graphs. `filter` is a
function taking an object of type `Graph` and returning a `Bool`.

# Example
We can produce an exhaustive list of `d`-regular graphs on `n` vertices by typing
`get_graphs(fill(d,n))`.

```julia-repl
julia> n, d = 6, 2
(6, 2)

julia> get_graphs(fill(d,n))
2-element Vector{Any}:
 {6, 6} undirected simple Int64 graph
 {6, 6} undirected simple Int64 graph
```

We can give the filter `is_connected` if we only want to consider connected graphs.

```julia-repl
julia> get_graphs(fill(d,n), is_connected)
1-element Vector{Any}:
 {6, 6} undirected simple Int64 graph
```
"""
function get_graphs(d_seq::Vector{Int}, filter::Function=((G::Graph)->true))

    d_seq = sort(d_seq) # Creates a sorted copy of the degree sequence.
    rel(i,j) = d_seq[i] == d_seq[j] # Optimization for non-regular searches.
    n0 = length(d_seq)
    A = falses(n0,n0)
    list = []

    # Convenience function.
    degree(i::Int) = @views count(A[i,i+1:end]) + count(A[1:i-1, i])

    # Define function for recursive search.
    function recur(n::Int, i::Int)
        degn = degree(n)

        # We found a graph! See if we add it to the list.
        if n == 1 && degn == d_seq[1]
            G = Graph(A + transpose(A))

            if !filter(G)
                return
            end

            if any(has_isomorph(G,H; vertex_relation=rel) for H in list)
                return
            end

            push!(list, G)

        # We finished vertex n. Move to the next and reset after.
        elseif degn == d_seq[n]
            recur(n-1, n-2)
            A[1:n-1,n-1] .= 0

        # Continue assigning edges to vertex n, if there's room.
        elseif i > 0 && i >= d_seq[n] - degn

            roomavailable = degree(i) < d_seq[i]
            @views identicaltoiplus1 = (i + 1 < n) && (A[i+1,n] == 0) && (d_seq[i+1] == d_seq[i]) && (A[i, n+1:end] == A[i+1, n+1:end])
            # ^^^ This line is an optimization.

            if roomavailable && !identicaltoiplus1
                A[i,n] = 1
                recur(n, i-1)
            end

            A[i,n] = 0
            recur(n, i-1)
        end
    end

    # Run search.
    recur(n0, n0-1)

    return list
end

"""
    get_reg_graphs(n::Int, d::Int, filter)

This is a convenience function equal to `get_graphs(fill(d,n), filter)`. It
returns an exhaustive list of non-isomorphic copies of `d`-regular graphs on
`n` vertices.
"""
get_reg_graphs(n::Int, d::Int, filter::Function=((G::Graph)->true)) = get_graphs(fill(d,n), filter)

"""
    get_graphs(n::Int, filter)

This is a convenience function which calls `get_graphs(d_seq, filter)` for each
nondecreasing sequence `d_seq` of `n` elements in the set {0,...,n-1}, and
returns the concatenation of the results.
"""
function get_graphs(n::Int, filter::Function=((G::Graph)->true))

    d_seq = zeros(Int, n)
    graph_list = []

    append_graphs() = append!(graph_list, get_graphs(d_seq, filter))

    append_graphs()
    while d_seq[1] < n-1
        i = n
        while d_seq[i] == n-1
            i -= 1
        end

        d_seq[i] += 1
        d_seq[i+1:n] .= d_seq[i]

        if isodd(sum(d_seq))
            continue
        end

        append_graphs()
    end

    return graph_list
end

"""
    get_orbits(G::Graph) -> orbits::Vector{Set{Int}}

Returns a list `orbits` of sets for which `orbits[i]` is the orbit of `i` under the action of the
automorphism group of G.
"""
function get_orbits(G::Graph)

    n = nv(G)
    orb = [Set(i) for i in 1:n]
    done = falses(n)

    for i in 1:n

        if done[i]
            continue
        end

        for j in i+1:n
            if (j in orb[i]) || done[j]
                continue
            end

            if has_isomorph(G, G, vertex_relation=((a,b)->(a==i)==(b==j)))
                union!(orb[i], orb[j])
                orb[j] = orb[i]
            end
        end

        for j in orb[i]
            done[j] = true
        end
    end

    return orb
end

"""
    get_loops(G::Graph)

Returns an n by n-1 matrix `L` of integers, where n is the number of vertices
in `G`, and where `L[i,k]` is the count of walks in `G` of length `k` which start and
end at vertex `i`.

The returned matrix `L` is used to determine if two vertices `i` and `j` are cospectral.
In particular, `i` and `j` are cospectral if and only if `L[i,:] == L[j,:]`.
"""
function get_loops(G::Graph)
    n = nv(G)
    A = Matrix{Int}(adjacency_matrix(G))
    B = LinearAlgebra.I
    L = zeros(Int, n, n-1)

    for i in 1:n-1
        B *= A
        for j in 1:n
            L[j,i] = B[j,j]
        end
    end

    return L
end

"""
    is_walk_regular(G::Graph)

Returns `true` if and only if the graph `G` is walk-regular.
"""
function is_walk_regular(G::Graph)
    L = get_loops(G)
    return all(L[1,:] == L[i,:] for i in 2:nv(G))
end

"""
    is_vertex_transitive(G::Graph)

Returns `true` if and only if the graph `G` is vertex-transitive.
"""
is_vertex_transitive(G::Graph) = length(get_orbits(G)[1]) == nv(G)

"""
    ncvs(G::Graph)

Returns an exhaustive list of pairs `(i,j)` of vertices in `G`, with `i < j`, which are
cospectral but not similar.
"""
function ncvs(G::Graph)
    n = nv(G)
    L = get_loops(G)
    orbits = get_orbits(G)
    @views return [(i,j) for i in 1:n for j in (i+1):n if ( L[i,:] == L[j,:] && !(j in orbits[i]) )]
end

"""
    has_ncvs(G::Graph)

Returns `true` if there are any pairs of vertices in `G` which are cospectral but not similar.
"""
has_ncvs(G::Graph) = !isempty(ncvs(G))

### READING AND WRITING FILES ###

# Filepaths
ncvs_examples_dir = "examples/ncvs/"
reg_examples_dir = "examples/regular/"
walkreg_examples_dir = "examples/walkregular/"

function save_csv(g::Graph, filename::String)
    # Convert adjacency matrix to string.
    A = adjacency_matrix(g)
    contents = ""
    rows, columns = size(A)
    for i in 1:rows
        for j in 1:columns
            contents *= string(A[i,j])
            if j < columns
                contents *= ","
            end
        end
        if i < rows
            contents *= "\n"
        end
    end

    # Save string as file.
    io = open(filename, "w")
    write(io, contents)
    close(io)
end

function save_ncvs_examples()
    i = 1
    for g in get_graphs(8, has_ncvs)
        save_csv(g, ncvs_examples_dir * "ncvs" * string(i) * ".csv")
        i += 1
    end
end

function save_reg_examples()
    i = 1
    for d in 0:9, g in get_reg_graphs(10, d, has_ncvs)
        save_csv(g, reg_examples_dir * "reg" * string(i) * ".csv")
        i += 1
    end
end

function save_walkreg_examples()
    i = 1
    criterion(g) = is_walk_regular(g) && !is_vertex_transitive(g)
    for d in 0:11, g in get_reg_graphs(12, d, criterion)
        save_csv(g, walkreg_examples_dir * "walkreg" * string(i) * ".csv")
        i += 1
    end
end

# To do: load graphs.

### CHECKS FOR VALIDITY ###

function test()
    graph_count = [1, 2, 4, 11, 34, 156, 1044, 12346]

    for n in 1:8
        message = "n = " * string(n) * " test: "
        message *= (length(get_graphs(n)) == graph_count[n] ? "SUCCESS" : "FAILED")
        println(message)
    end
end

function reg_test()
    reg_graph_count = [1, 2, 2, 4, 3, 8, 6, 22, 26, 176, 546, 19002]

    for n in 1:12

        tally = 0
        for d in 0:n-1
            tally += length(get_reg_graphs(n,d))
        end

        message = "n = " * string(n) * " test: "
        message *= (tally == reg_graph_count[n] ? "SUCCESS" : "FAILED")
        println(message)
    end
end
