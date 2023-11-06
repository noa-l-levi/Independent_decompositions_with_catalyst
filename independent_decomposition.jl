# loading the required packages
using RowEchelon, Catalyst, Graphs, Latexify, Combinatorics, MetaGraphs, GraphPlot, DelimitedFiles

# independent decomposition function
function independent_decomposition(reaction_network)
    # create a reaction diagram
    # complexgraph(reaction_network)

    # compute the deficiency
    reactioncomplexes(reaction_network)
    linkageclasses(reaction_network)
    numspecies(reaction_network)
    N = netstoichmat(reaction_network)
    δ = deficiency(reaction_network)

    # get the transpose of the stoichiometric matrix
    R = transpose(N)

    # compute a maximal linearly independent set of vectors in R
    pivots = rref_with_pivots(transpose(R))

    basis = R[pivots[2],:]

    #initialise an undirected graph G
    G = SimpleGraph()
    G = MetaGraph(G)

    # add vertices to G which are the vectors is basis
    for i = 1:length(pivots[2])
        add_vertex!(G)
        st = string("R",pivots[2][i])
        set_prop!(G, i, :name, st)
    end

    # initialise matrix of linear combinations (#rxns x #basis elements)
    linear_combo = zeros(size(N)[2], length(pivots[2]))

    # get matrix which shows the linear combinations of basis rxns that 
    # give the nonbasis rxns (basis rxns will have a row of zeros)
    for i = 1:size(N)[2]
        if !in.(i, Ref(pivots[2]))

            linear_combo[i,:] = transpose(basis)\R[i,:]
        end
    end

    linear_combo = round.(linear_combo, digits=1) 

    # Get the reactions that are linear combinations of at least 2 basis reactions
    # These are the reactions where we'll get the edges
    get_edges = first.(Tuple.(findall(x -> x > 1, abs.(sum(linear_combo, dims=2))))) # problem is in here I think- for S11, r7 is missing because r7 = r5-r6 and 1-1 gets zero

    #  Initialize an array for sets of vertices that will form the edges
    vertex_set = Array{}[];
        
    #  Identify which vertices form edges in each reaction: 
    # get those with non-zero coefficients in the linear combinations
    for i = 1:length(get_edges)
        push!(vertex_set, [first.(Tuple.(findall(!iszero, linear_combo[get_edges[i], :])))]);
    end
        
    # Initialize the edge set
    edges = Vector{}[];
        
    #  Get all possible combinations (not permutations) of the reactions 
    # involved in the linear combinations
    for i = 1:length(vertex_set)
        append!(edges, [collect(combinations(vertex_set[i][1], 2))]);
    end
        
    #  Get just the unique edges (if edges is non-empty)
    if isempty(edges) != true
        edges = unique(reduce(vcat,edges));
    end
        
    # Add these edges to graph G
    for i = 1:length(edges)
        add_edge!(G, edges[i][1], edges[i][2]); 
    end

    # gplot(G)

    # Determine whether a finer decomposition exists, if not, end here
    if is_connected(G) == true
        println("No finer decomposition exists")
        return
    end
    
    if is_connected(G) != true
        println("A finer network decomposition exists")
    end

    # Determine which vertices are connected
    components = connected_components(G);
        
    # Determine the number of connected components of G: this is 
    # the number of partitions R will be decomposed to
    println("The number of independent subnetworks is:")
    num_components = size(components,1)
    display(num_components)
    # Initialize the vector of partitions
    P = [];
        
    # Basis vectors: assign them first into their respective partition 
    for i = 1:size(components)[1]
        push!(P,pivots[2][components[i][1:size(components[i])[1]]]);
    end
    
    # Nonbasis vectors: they go to the same partition as the basis vectors 
    # that form their linear combination
    for i = 1:size(P,1)
        for j = 1:size(P[i],1)
            
            # Get the column number representing the basis vectors in 'linear_combo'
            col = findall(x -> x == P[i][j], pivots[2])[1];
            
            # Check which reactions used a particular basis vector and assign them 
            # to their respective partition
            append!(P[i], findall(x -> x != 0, linear_combo[:, col])); 
        end
    end
    
    # Get only unique elements in each partition
    for i = 1:size(P,1)
        P[i] = unique(P[i]);
    end
    println("The partition of reactions is:")
    display(P)

    # # Send list of partitions to a csv file
    # s1 = "reaction_network"
    # s2 = ".csv"
    # filename = "s1,s2" 
    # writedlm(filename,  P, ',')

    # Check that all reactions have been assigned a partition
    if length(reduce(vcat,P)) == size(R, 1) 
        println("No reactions have been omitted in the partition") 
        # return
    end
    
    #If reactions are missing, repeat from the linear combination step
    if length(reduce(vcat,P)) != size(R, 1) 
        println("Warning: some reactions have been omitted in the partition") 
        missing = size(R, 1) - length(reduce(vcat,P))
        display(missing)
   \
        # #Do not round
        # linear_combination = linear_combo
        # get_edges = first.(Tuple.(findall(x -> x > 1, abs.(sum(linear_combination, dims=2))))) 
        # vertex_set = Array{}[];
        # for i = 1:length(get_edges)
        #     push!(vertex_set, [first.(Tuple.(findall(!iszero, linear_combination[get_edges[i], :])))]);
        # end
        # edges = Vector{}[];
        # for i = 1:length(vertex_set)
        #     append!(edges, [collect(combinations(vertex_set[i][1], 2))]);
        # end
        # if isempty(edges) != true
        #     edges = unique(reduce(vcat,edges));
        # end
        # for i = 1:length(edges)
        #     add_edge!(G, edges[i][1], edges[i][2]); 
        # end
        # if is_connected(G) == true
        #     println("No finer decomposition exists")
        #     return
        # end
        
        # if is_connected(G) != true
        #     println("A finer network decomposition exists")
        # end
        # components = connected_components(G);
        # println("The number of independent subnetworks is:")
        # num_components = size(components,1)
        # display(num_components)
        # P = [];
        # for i = 1:size(components)[1]
        #     push!(P,pivots[2][components[i][1:size(components[i])[1]]]);
        # end
        # for i = 1:size(P,1)
        #     for j = 1:size(P[i],1)
        #         col = findall(x -> x == P[i][j], pivots[2])[1];
        #         append!(P[i], findall(x -> x != 0, linear_combination[:, col])); 
        #     end
        # end
        # for i = 1:size(P,1)
        #     P[i] = unique(P[i]);
        # end
        # println("The partition of reactions is:")
        # display(P)
        # if length(reduce(vcat,P)) == size(R, 1) 
        #     println("No reactions have been omitted in the partition") 
        #     # return
        # end
    end

    
end 

# define the reaction network
test_cap = @reaction_network begin
    k1, A + B --> B + C
    (k2,k3), B + C <--> 2B
    (k4,k5), B <--> 2E
    k6, 2E --> 2D
    k7, C --> A
    k8, D --> E
end

example8 = @reaction_network begin
    k1, X2 --> X1 + X2
    k2, X1 + X5 --> X2 + X5
    k3, 2*X5 + X1 --> X5 + X1
    k4, X2 + X5 --> X3 + X5
    k5, 2*X5 + X2 --> X5 + X2
    k6, X2 + X5 --> X5
    k7, X2 + X5 --> X2
    k8, X3 + X5 --> X4 + X5
    k9, X3 + X5 --> X3 + 2*X5
    k10, X3 + X4 + X5 --> X4 + X5
    k11, X3 + X4 + X5 --> X3 + X5
    k12, X3 + X4 + X5 --> X3 + X4 + 2*X5
    k13, 2*X5 --> X5
end

example9 = @reaction_network begin
    k1, 0 --> X1
    k2, X1 + X3 --> X3 + X2
    k3, X2 --> X3
    k4, X1 + X2 --> X1 + X4
    k5, X3 --> 0
    k6, X4 --> 0
end

EnvZOmpR = @reaction_network begin
    (a5,d5), XD <--> X
    (a1,d1), X <--> XT
    k1, XT --> Xp
    (a2,d2), Xp + Y <--> XpY
    k2, XpY --> X + Yp
    (a3,d3), XT + Yp <--> XTYp
    k3, XTYp --> XT + Y
    (a4,d4), XD + Yp <--> XDYp
    k4, XDYp --> XD + Y
end

# Run the decomposition on a reaction network
independent_decomposition(test_cap)
    
full_hierarchy = @reaction_network begin
    (a1,d1), E + I <--> EI
    (a2,d2), EI + S <--> ESI
    (a3,d3), ESI <--> ES + I
    (a4), ESI --> EI + P
end

newEO = @reaction_network begin
    (k1,k2), E1 <--> E2
    k3, E2 --> E3
    (k5,k6), E1 + S1 <--> C1
    k7, C1 --> E1 + S18
    (k8,k9), E2 + S18 <--> C2 
    k10, C2 --> E2 + S1 
    (k11,k12), E3 + S2 <--> C3
    k13, C3 --> E3 + S28
end

Sred = [[0 1 -1 -1];[1 -1 0 1];[-1 0 1 0];[-1 0 0 0];[0 -1 0 0];[0 0 1 0];[0 0 0 1]]

S11_opposer = @reaction_network begin
    k1, R --> R + X1
    k2, X1 + X2 --> 0
    k3, 0 --> X2
    k4, X3 --> X3 + X2
    k6, X2 + X3 --> X2
    k7, O1 + X3 --> O1 + O1
    k8, O1 --> 0
    k5, 0 --> X3
end