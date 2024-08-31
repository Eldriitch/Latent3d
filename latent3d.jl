module Latent3D
using LinearAlgebra, Graphs, Polyhedra, GLPK, PlotlyJS

const eps = 0.05
const floattol = 0.01

  ### Basic definitions ###
# Struct for all scatterers, so things are appropriately extensible should we 
# wish to relax assumptions of the dilute model. 
struct Scatterer 
    position::Array{Float64}
    cap::Float64
end
function Scatterer(v::Vector{Float64})
    Scatterer(v, 4*pi)
end
function Scatterer(x::Real, y::Real, z::Real)
    Scatterer([Float64(x),Float64(y),Float64(z)], 4*pi)
end
function Scatterer(x::Real, y::Real, z::Real, R::Real)
    Scatterer([Float64(x),Float64(y), Float64(z)], 4*pi*R)
end
Base.:(==)(a::Scatterer, b::Scatterer) = a.position == b.position && a.cap == b.cap
Base.hash(a::Scatterer, h::UInt) = hash(a.cap, hash(a.position, h))

# Function to generate a configuration of N distinct identical scatterers.
# For convenience the configuration is centered on the centroid.
function generate_configuration(N::Integer, Rmax::Integer)
    scattererarray = [Scatterer(0,0,0) for i in 1:N]
    while length(unique([i.position for i in scattererarray])) != length(scattererarray)
        for i in 1:N
            scattererarray[i] = Scatterer(rand(0:5), rand(0:5), rand(0:5), rand(1:Rmax))
        end
    end
    c = compute_centroid(scattererarray)
    for i in 1:N
        scattererarray[i] = Scatterer(scattererarray[i].position - c, scattererarray[i].cap)
    end
    scattererarray
end
# Function to compute the capacitance matrix of a configuration
function compute_capacitance_matrix(scattererarray::Array{Scatterer})
    N = length(scattererarray)
    capacitancematrix = Array{Float64}(undef, N, N)
    for i in 1:N
        for j in 1:N
            if i == j
                capacitancematrix[i,j] = eps*scattererarray[i].cap
            else 
                r = scattererarray[i].position - scattererarray[j].position
                d = sqrt(r[1]^2+r[2]^2+r[3]^2)
                capacitancematrix[i,j] = -1*eps^2*scattererarray[i].cap^2/(4*pi*d)
            end
        end
    end
    capacitancematrix
end
# Function to check for latent symmetry
function checkforlatentsymmetry(scattererarray::Array{Scatterer})
    N = length(scattererarray)
    C = compute_capacitance_matrix(scattererarray)
    L = N > 10 ? 10 : N-1
    M = Array{Float64}(undef, N, L)
    for j in 1:L
        powerC = C^j
        for i in 1:N
            M[i,j] = powerC[i,i]
        end
    end
    if M != unique(M, dims=1)
        return [M, find_duplicate_rows(M)]
    else
        return nothing
    end
end
  ### Symmetry Detection ###
# Partitions the array on distances from the origin.
function distance_partition(scattererarray::Array{Scatterer})
    Y = Dict(D => (D.position[1]^2+D.position[2]^2+D.position[3]^2) for D in scattererarray)
    V = values(Y)
    X = Dict{Float64, Array{Scatterer}}()
    for v in V
        X[v] = findall(x -> x == v, Y)
    end
    X
end
# Extracts the graph of the convex hull of a set of points
function extract_graph(positions)
    hull = vrep(positions)
    halfspacerep = doubledescription(hull)
    graph = Graph()
    p = collect(points(hull))
    n = length(p)
    for v in p
        add_vertex!(graph)
    end
    planes = collect(hyperplanes(halfspacerep))
    for I in halfspaces(halfspacerep)
        iplane = HyperPlane(I.a, I.β)
        push!(planes, iplane)
    end
    for I in planes
        for J in planes
            if I != J
                cap = I ∩ J
                for v in 1:n
                    for w in 1:n
                        if v != w && (p[v] in cap && p[w] in cap)
                            add_edge!(graph, v, w)
                        end
                    end
                end
            end
        end
    end
    graph
end
# converts a Graphs.jl graph into the format used by nauty
const WORDSIZE = 64
function graph_to_nautygraph(g::Graph)
    num_vertices = nv(g)
    num_setwords = div(num_vertices - 1, WORDSIZE) + 1
    graph_arr = BitArray(undef, num_setwords * WORDSIZE, num_vertices)
    fill!(graph_arr, false)
    for (j, i) in enumerate(g.fadjlist)
        for value in i
            graph_arr[end-value+1,j] = true
        end
    end
    graph_arr.chunks
end
# Calls to nauty to compute the automorphism group of a graph, if enough space
# has been allocated. 
function find_automorphisms(g::Graph)
    g_nauty = graph_to_nautygraph(g)
    n = nv(g)
    target_array_size = 100*n
    target_array = Vector{Cint}(undef, target_array_size)
    target_array = - ones(Cint, target_array_size)
    @ccall "./nautyinterface.so".find_automorphism_group(
        g_nauty::Ref{Culong},
        n::Cint,
        target_array::Ref{Cint},
        target_array_size::Culonglong
        )::Cvoid
    unique(reshape(target_array, (n, :)), dims=2)[:, 1:end-1]
end
# Finds all of the symmetries of af an array of scatterers
# For arrays with full rank, we have two cases, where there are
# only n=1,2 distance sets, and when there is at least one n>=3 distance set.
# in both cases, we take a convex hull of a good subset, which will have finitely many
# automorphisms. Then we find the associated matrices, and check them. 
function find_all_symmetries(scattererarray::Array{Scatterer})
    all_positions = scatterer_positions(scattererarray)
    Y = distance_partition(scattererarray)
    symmetries = Array{Matrix{Float64}}(undef, 0)
    if !any([length(D)>=3 for D in Y])
        positions = all_positions
        v = vrep(positions)
        v = removevredundancy(v, GLPK.Optimizer)
        positions = collect(points(v))
        g = extract_graph(positions)
    else
        F = filter(x -> length(x)>=3, D)
        i = argmin(length, keys(F))
        positions = scatterer_positions(F[i])
        g = extract_graph(positions)
    end
    perms = find_automorphisms(g)
    for i in 1:size(perms)[2]
        p = perms[:, i]
        M = matrix_of_permutation(positions, p)
        # print(M); print("\n")
        # print(is_orthogonal(M)); print("\n")
        # print(approx_equal_perm(map(v->M*v, all_positions), all_positions)); print("\n")
        if is_orthogonal(M) && approx_equal_perm(map(v ->M*v, all_positions), all_positions)
            push!(symmetries, M)
        end
    end
    symmetries
end
  ### Scripts Functions ###
function generate_latent_config()
    while true
        S = generate_configuration(5, 1)
        if !has_full_rank(S) || checkforlatentsymmetry(S) isa Nothing || length(find_all_symmetries(S)) > 1 
            continue 
        end
        return S
    end
end
function long_run()
    open("commutator.txt", "w") do file
        i = 1
        while true
            M = generate_configuration(5, 1)
            Q = compute_capacitance_matrix(M) 
            a = rand(1:5)
            b = rand(1:5)
            while a != b 
                b = rand(1:5)
            end
            R = compute_isospectral_reduction(0.5, [rand(1:5), rand(1:5)], Q)
            X = commute_space(R)
            if length(X) > 2
                println(i)
                println(file, M)
                println(file, Q)
                println(file, X)
                i += 1
            end
        end
    end
end
function fixed_config_test()
    A = []
    open("fixedconfighexagon.txt", "w") do file
        for i in -5:5, j in -5:5, k in -5:5, u in -5:5, v in -5:5, w in -5:5
            if [i,j,k] == [u,v,w]
                continue
            end
            X = [Scatterer(1.0,0.0,-1.0),Scatterer(1.0,-1.0,0.0),Scatterer(0.0,-1.0,1.0), Scatterer(-1.0,0.0,1.0), Scatterer(-1.0,1.0,0.0), Scatterer(0.0,1.0,-1.0)]
            # X = [Scatterer(1.2, -0.4, -0.2), Scatterer(0.2, 0.6, -0.2), Scatterer(0.2, -0.4, 0.8), Scatterer(-1.3, -0.9, 1.3), Scatterer(-0.3, 1.1, -1.7)]
            # X = [Scatterer(2.5,2.5,2.5),Scatterer(3.5,2.5,2.5), Scatterer(1.5,2.5,2.5), Scatterer(2.5,3.5,2.5), Scatterer(2.5, 1.5,2.5)]
            # X = [Scatterer(2.5 + 0.5, 2.5 - sqrt((5+2*sqrt(5))/20), 2.5),
                 # Scatterer(2.5 - 0.5, 2.5 - sqrt((5+2*sqrt(5))/20), 2.5), 
                 # Scatterer(2.5 + (1+sqrt(5))/4, 2.5 + sqrt((5 - sqrt(5))/40), 2.5),
                 # Scatterer(2.5 - (1+sqrt(5))/4, 2.5 + sqrt((5 - sqrt(5))/40), 2.5),
                 # Scatterer(2.5, 2.5 + sqrt((5 + sqrt(5))/10), 2.5)]
            push!(X, Scatterer(i,j,k))
            push!(X, Scatterer(u,v,w))
            C = compute_centroid(X)
            X = [Scatterer(D.position - C, D.cap) for D in X]
            if !has_full_rank(X) || checkforlatentsymmetry(X) isa Nothing || length(find_all_symmetries(X)) > 1
                continue
            end
            println(file, "A: ", [i,j,k])
            println(file, "B: ", [u,v,w])
            println(file, checkforlatentsymmetry(X)[2][1])
            push!(A, X)
        end
    end
    A
end
function generalised_symmetry_run()
    open("generalised_symmetry_run.txt", "w") do file
        for k in 1:100
            println(k)
            while true
                n = rand(5:15)
                M = generate_configuration(n, 1)
                C = compute_capacitance_matrix(M)
                U = Array{Float64}(undef, (n, n-1))
                for j in 1:n-1
                    powerC = C^j
                    for i in 1:n
                        U[i,j] = powerC[i,i]
                    end
                end
                if isdisjoint(eachrow(U), eachrow(2*U))
                    continue
                end
                println(file, M)
                println(file, C)
                break
            end
        end
    end
end
  ### Helper Functions ### 
# Computes the centroid
function compute_centroid(scattererarray::Array{Scatterer})
    C = zeros(Float64, 3)
    for D in scattererarray
        C = C + D.position
    end
    C ./ float(length(scattererarray))
end
function scatterer_positions(scattererarray::Array{Scatterer})
    N = length(scattererarray)
    positions = Array{Vector{Float64}}(undef, N) 
    for i in 1:N
        positions[i] = scattererarray[i].position
    end
    positions
end
function is_orthogonal(M::Matrix{Float64})
     isapprox(M*M',I; atol=floattol)
end
function matrix_of_permutation(position_array::Array{Vector{Float64}}, perm::Vector{Int32})
    N = length(position_array)
    ivals = rand(1:N, 3)
    u1 = position_array[ivals[1]]
    u2 = position_array[ivals[2]]
    u3 = position_array[ivals[3]]

    while isapprox(det([u1 u2 u3]), 0.0; atol=floattol)
        ivals = rand(1:N, 3)
        u1 = position_array[ivals[1]]
        u2 = position_array[ivals[2]]
        u3 = position_array[ivals[3]]
    end
    v1 = position_array[perm[ivals[1]]+1]
    v2 = position_array[perm[ivals[2]]+1]
    v3 = position_array[perm[ivals[3]]+1]
    # print([u1 u2 u3]); print("\n")
    # print(det([u1 u2 u3])); print("\n")
    M = [v1 v2 v3]*inv([u1 u2 u3])
    M
end
function has_full_rank(scattererarray::Array{Scatterer})
    p = scatterer_positions(scattererarray)
    v = vrep(p)
    h = doubledescription(v)
    H = polyhedron(hrep(allhalfspaces(h)))
    detecthlinearity!(H, optimizer_with_attributes(GLPK.Optimizer, "presolve"=>GLPK.GLP_ON))
    dim(H, current=true) == 3
end
function approx_equal_perm(u::Vector{Vector{Float64}}, w::Vector{Vector{Float64}})
    N = length(u)
    # print([[norm(u[i]-w[j]) for j in 1:N] for i in 1:N]); print("\n")
    all([norm([norm(u[i]-w[j]) for j in 1:N],-Inf)<floattol for i in 1:N])   
end
function find_duplicate_rows(matrix::Array{Float64, 2})
    n, m = size(matrix)
    row_map = Dict{Vector{Float64}, Vector{Int}}()

    for i in 1:n
        row = matrix[i, :]
        if haskey(row_map, row)
            push!(row_map[row], i)
        else
            row_map[row] = [i]
        end
    end
    duplicate_indices = []
    for (row, indices) in row_map
        if length(indices) > 1
            push!(duplicate_indices, indices)
        end
    end
    duplicate_indices
end
function plot_scatterer(scattererarray::Array{Scatterer}, latentindexes)
    position_array = scatterer_positions(scattererarray)
    N = length(position_array)
    xs = [position_array[i][1] for i in 1:N]
    ys = [position_array[i][2] for i in 1:N]
    zs = [position_array[i][3] for i in 1:N]
    latentxs = [position_array[i][1] for i in latentindexes]
    latentys = [position_array[i][2] for i in latentindexes]
    latentzs = [position_array[i][3] for i in latentindexes]

    trace1 = scatter3d(;x=xs, y=ys, z=zs,
        marker = attr(color="#0000bb", size=12),
        mode = "markers")
    trace2 = scatter3d(;x=latentxs, y=latentys, z=latentzs,
        marker = attr(color="#bb0000", size=12),
        mode = "markers")
    layout = Layout(width=600, height=600,
        autosize=false,
        showlegend=false,
        scene=attr(aspectmode = "cube",
#        camera=attr(projection = attr(type = "orthographic"),
#                   eye = attr(x=0.0, y=0.0, z=2.5),
#                    up = attr(x=0.0, y=1.0, z=0.0)),
        xaxis=attr(range=[-5,5]),
        yaxis=attr(range=[-5,5]),
        zaxis=attr(range=[-5,5])),
        margin=attr(l=0, r=0, b=0, t=0))
    plot([trace1, trace2], layout)
end 
function compute_isospectral_reduction(lambda, S, C)
    N = size(C)[1]
    T = [i for i in 1:N if !(i in S)]
    U = [C[i,j] for i in S, j in S]
    V = [C[i,j] for i in S, j in T]
    W = [C[i,j] for i in T, j in S]
    X = [C[i,j] for i in T, j in T]
    U - V*inv(lambda*I - X)*W
end
function degree_of_commutativity(R, M)
    opnorm(R*M - M*R, 2)
end
function perturb_scatterer!(scattererarray::Array{Scatterer}, index, scale)
    r = scale*randn(3)
    scattererarray[index] = Scatterer(scattererarray[index].position + r, scattererarray[index].cap)
    norm(r,2)
end
function direction_perturb_scatterer!(scattererarray::Array{Scatterer}, index, scale, direction)
    n = normalize(direction)
    r = scale*n
    scattererarray[index] = Scatterer(scattererarray[index].position + r, scattererarray[index].cap)
    norm(r,2)
end
function plot_stability(scattererarray::Array{Scatterer}, index, latent_indices)
    xs = [] 
    ys = []
    v = randn(3)
    for i in 1:1000
        S = copy(scattererarray)
        r = direction_perturb_scatterer!(S, index, 0.3*i/1000, v)
        if r == 0.0
            continue
        end
        push!(xs, r)
        C = compute_capacitance_matrix(S)
        R = compute_isospectral_reduction(rand(), latent_indices, C)
        P = [0 1; 1 0]
        gamma = degree_of_commutativity(R,P)
        push!(ys, gamma)
    end
    trace = scatter(; x = xs, y = ys, mode="markers")
    println(minimum(xs))
    println(minimum(ys))
    linearxs = 0.3/1000*collect(1:1000)
    linearys = 1e-3*linearxs
    # quadraticys = linearxs .^ 2
    trace1 = scatter(; x = linearxs, y = linearys)
    # trace2 = scatter(; x = linearxs, y = quadraticys)
    layout = Layout(
        xaxis = attr(type = "log"),
        yaxis = attr(type = "log")
    )
    plot([trace, trace1], layout)
    # plot(trace)
end
function plot_eigenmodes(scattererarray::Array{Scatterer}, latentindices)
    N = length(scattererarray)
    traces = Array{AbstractTrace}(undef, N)
    v = randn(3)
    for q in 1:N
        xs = []
        ys = []
        for i in 0:1000
            S = copy(scattererarray)
            r = i*3/500 - 3
            direction_perturb_scatterer!(S, latentindices[1], r, v)
            push!(xs, r)
            C = compute_capacitance_matrix(S)
            vecs = eigvecs(C)
            push!(ys, vecs[latentindices[1],q]/vecs[latentindices[2],q]) 
        end
        traces[q] = scatter(; x = xs, y = ys, mode="markers")
    end
    layout = Layout(
        xaxis = attr(range=[-3, 3]),
        yaxis = attr(range=[-1.25, 1.25]),
        showlegend=false,
    )
    plot(traces, layout)
end
function plot_latent_curve()
    xs1 = [1, 0, 0]
    ys1 = [0, 1, 0]
    zs1 = [0, 0, 1]
    trace1 = scatter3d(; x = xs1, y = ys1, z = zs1, 
                       mode = "markers")
    a(t) = 1 + cos(t)/(2*sqrt(2))
    b(t) = 2 + sin(t)/2
    c(t) = 1 + cos(t)/(2*sqrt(2))
    xs2 = [a(t*pi/500) for t in 0:1000]
    ys2 = [b(t*pi/500) for t in 0:1000]
    zs2 = [c(t*pi/500) for t in 0:1000]
    trace2 = scatter3d(; x = xs2, y = ys2, z = zs2,
                       mode = "markers",
                       marker = attr(color = collect(0:1000)))
    xs3 = (xs2 .- (2 .* ys2) .- (2 .* zs2) .+ 2)./3 
    ys3 = ((-2 .* xs2) .- (2 .* ys2) .+ zs2 .+ 2)./3
    zs3 = ((-2 .* xs2) .+ ys2 .- (2 .* zs2) .+ 2)./3
    trace3 = scatter3d(; x = xs3, y = ys3, z = zs3,
                       mode = "markers",
                       marker = attr(color = collect(0:1000)))
    layout = Layout(
        autosize = false,
        showlegend = false,
        scene = attr(aspectmode = "cube",
                xaxis = attr(range = [-2, 3]),
                yaxis = attr(range = [-2, 3]),
                zaxis = attr(range = [-2, 3]),
#                camera = attr(projection = attr(type = "orthographic"),
#                              eye = attr(x = 2.5, y = 0.0, z = 0.0),
#                              up = attr(x = 0.0, y = 0.0, z = 1.0),)
            )
        )

    plot([trace1, trace2, trace3], layout)
end
function commute_space(A)
    n = size(A)[1]
    map(v -> reshape(v, n,n),
    eachcol(nullspace(kron(I(n), A) - kron(A', I(n)))))
end
function matrix_shape(M)
    a  = unique(M)
    Q = Array{Char}(undef, size(M))
    for i in eachindex(a)
        for j in eachindex(M)
            if a[i] == M[j]
                Q[j] = i - 1 + 'a' 
            end
        end
    end
    return Q
end
end
