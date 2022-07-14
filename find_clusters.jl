using Plots
using Combinatorics
using Spacey
using LinearAlgebra

""" Check and see if the all the entries in each column of m are equal """
function isCommonShift(m)
    return all(y->y≈m[:,1],eachcol(m)) 
end

""" Get the atom positions of all structures. Return cell sizes as well """
function get_dsites_cell_sizes()
    cd("/Users/glh43/home/fortranCodes/uncle/cluster_tests/")
    f = readlines("structures.in")
    s = findall(x->occursin("Direct",x),f).+1
    e = findall(x->occursin("Energy",x),f).-2
    idx = hcat(s,e)
    cell_sizes = [j[2]-j[1]+1 for j in eachrow(idx)]
    # Collect all atom positions for each structure
    dsites = []
    for j in 1:size(idx,1)
        nextSites = [parse.(Float64,split(i)[1:3]) for i in f[idx[j,1]:idx[j,2]]]
        # Get lattice vectors and convert from direct to Cartesian coordinates
        S = f[idx[j,1].-collect(5:-1:3)]
        S = hcat([parse.(Float64,split(s)) for s ∈ S]...)
        nextSites = S*hcat(nextSites...)
        nextSites = [nextSites[:,i] for i ∈ 1:size(nextSites,2)]
        push!(dsites,nextSites)
    end
    return dsites, cell_sizes
end

""" Get the (integer) symmetry operations of the lattice. Don't include identity """
function get_lattice_symops()
    cd("/Users/glh43/home/fortranCodes/uncle/cluster_tests/")
    g = readlines("lat.in")
    A = hcat([parse.(Float64,split(i)) for i ∈ g[6:8]]...)
    u,v,w = eachcol(A)
    rops,ops = pointGroup_robust(u,v,w)
    filter!(i->i!=I(3),ops) # Take the identity out of the list
    return ops
end

""" Expand a list of positions/structure to cluster combinations """
function generateFiguresFromSites(dsites)
    # For each structure, generate all candidate clusters of all vertex orders
    candClust = [collect(combinations(jSt)) for jSt ∈ dsites]
    # Undo one layer of nesting (structure-by-structure nesting not needed now)
    # candClust will now be a vector of clusters (a cluster is a vector of vectors)
    candClust = vcat(candClust...)
    # Toss out single-vertex cluster candidates
    candClust = filter(x->length(x)> 1,candClust) 
    # Convert from vectors of vectors to matrices 
    candClust = [hcat(i...) for i ∈ candClust]
    # Sort the columns, put vertices in a canonical order (handy for comparisons)
    candClust = [sortslices(i,dims=2) for i ∈ candClust]
    # Shift all clusters so that vertex 1 is the origin
    candClust = [i.-i[:,1] for i ∈ candClust]
    # Toss out exact duplicates (tossing out other kinds of duplicates will happen later)
    candClust = unique(candClust)
    # Sort/group clusters by number of vertices (not necessary but aesthetically pleasing) 
    candClust = candClust[sortperm(size.(candClust,2))]
    return candClust
end

""" Attach s-vectors to figures """
function attachSvectors(clusters)
    sV = [permutations(vcat([i*ones(Int,j) for i in 1:2]...),j)|>collect|>unique for j in 2:6]
    svec = []
    scl = []
    for iCl ∈ clusters
        l = size(iCl,2)
        for iSvec ∈ 1:2^l # Loop over all the s-vectors for this vertex order 
            push!(scl,iCl) # Just another copy of the geometric figure
            push!(svec,sV[l-1][iSvec])
        end
    end
    return scl, svec
end

symops = get_lattice_symops()
d, Nc = get_dsites_cell_sizes()

## Try it with only 1- & 2-site cells
d2 = unique(d[Nc .< 3])
# Convert to matrices (drop nV < 2 clusters)
cl2 = generateFiguresFromSites(d2)

cl2, svec = attachSvectors(cl2)



mask = trues(size(candClust,1))
for (i,iCl) ∈ enumerate(candClust)
    for (j,jCl) ∈ enumerate(candClust[1:i-1])
        if length(iCl) == length(jCl) && iCl≈jCl
            mask[i] = false
        end
    end
end
# 638 is reduced to 408 when 'isapprox' is used to compare
candClust = candClust[mask]
plot(length.(candClust),title="Unique Candidate figures",legend=:none)

# Remove clusters that are non-zero lattice translations of each other
uqClust = []
for (j,jCl) ∈ enumerate(candClust)
    uq = true
    for (i,iCl) ∈ enumerate(uqClust)
        if length(iCl) == length(jCl)
            r = iCl - jCl # Turn r into a matrix
            # All all elements of each row equal?
            if isCommonShift(r) 
                uq = false
                break
            end
        end
    end
    if uq
        push!(uqClust,jCl)
    end
end
size(uqClust) # 408->401

# Finally we need to eliminate clusters that are equivalent under rotation (and rotation+lattice shift)
mask = trues(size(uqClust,1))
for (i,iCl) ∈ enumerate(uqClust)
    for (j,jCl) ∈ enumerate(uqClust[1:i-1])
        if length(iCl) != length(jCl) 
            continue
        end
        for iRot ∈ ops
            t = sort(transpose(transpose(iCl)*iRot),dims=2)
            r = t - jCl # Compute translation difference between the vertices
            if t ≈ jCl # then iCl is a rotation duplicate
                mask[i] = false
                break
            #All all elements of each row equal?
            elseif isCommonShift(r) # Then iCl is a rotated translation duplicate
                mask[i] = false
                println("shift dup")
                break
            end
        end
        if !mask[i] break end
    end
end
rotuqClust = uqClust[mask]
# 277 of these found July 11 8:50 pm
# 267 at 9:40 (added sorting of vertices to 't')
# 108 at 9:30 July 12 (shifted all clusters to have vertex 1 at origin)...but now they are unsorted...
# 314 at 2 am Juli 13

nFigPerOrder = [sum(length.(rotuqClust).==i*3) for i in 2:6]
nClustPerOrder = [2^(i+1)*nFigPerOrder[i] for i ∈ 1:length(nFigPerOrder)]
println("Number of figures and clusters at each order (2-6): \n",nFigPerOrder,"\n",nClustPerOrder)
println("Total number of clusters (including s-vectors): ",sum(nClustPerOrder))

##Debugging
rotuqClust = copy(candClust)


## Now that we have all the geometric clusters we need to attach all possible s-vectors
sV = [permutations(vcat([i*ones(Int,j) for i in 1:2]...),j)|>collect|>unique for j in 2:6]
sCl = [] # Clusters with s-vectors
sVlist = []
for iCl ∈ rotuqClust
    l = size(iCl,2)
    for iSvec ∈ 1:2^l # Loop over all the s-vectors for this vertex order 
        push!(sCl,iCl) # Just another copy of the geometric figure
        push!(sVlist,sV[l-1][iSvec])
    end
end
size(sCl)
# 4784 of these on July 11 9:30 pm
# 4664 at 9:45 pm
# 2344 at 9:30 July 12
# 5700 at 2 am July 13

## Now that we have the s-vectors, we need to eliminate s-clusters that are symmetrically equivalent
mask = trues(size(sCl,1))
for (i,iCl) ∈ enumerate(sCl[1:end])
    println(i)
    for (j,jCl) ∈ enumerate(sCl[1:i-1])
        if length(iCl) != length(jCl) 
            continue
        end
        for iRot ∈ ops
            t = transpose(transpose(iCl)*iRot) # Rotate the i-th sCl
            p = sortperm(collect(eachslice(t,dims=2))) # how t is permuted by the sort
            t = t[:,p] # sort t
            svec = sVlist[i][p] # Permute the labels on the vertices
            if t ≈ jCl && svec==sVlist[j]
                mask[i] = false
                break
            end
            r = t - jCl
            if isCommonShift(r) && svec==sVlist[j]
                mask[i] = false
                println("breaking")
                break
            end
        end
        if !mask[i] break end
    end
end

## Debugging
finalCl = copy(sCl)
finalSV = copy(sVlist)

finalCl = sCl[mask]
#finalCl = [i.-i[:,1] for i ∈ finalCl]
finalSV = sVlist[mask]
#finalCl = finalCl[1:10]
#finalSV = sVlist[1:10]
#4361 on July 12 morning

## Write a clusters.out-type file with these clusters

# Copy the header and point clusters

""" Write out clusters & s-vectors to clusters.out-type file """
function(clusters, svectors)
    str = """
    #--------------------------------------------------------------------------------
    # Cluster number: 0
    # This is the empty cluster = constant term for the cluster expansion.
    #
    # - The empty cluster is not explicitely listed in this file. 
    #   (This comment here only reminds the user of the fact that
    #   the empty cluster does exist in UNCLE.)
    # - The emtpy cluster does not contain any vertices.
    # - The empty cluster is not read in from this file.
    # - The emtpy cluster is always used in the fitting procedure.
    # - The emtpy cluster counts as a separate cluster for the fitting procedure.
    #   (E.g., a simple fit with 1 cluster is a fit with the empty cluster only.)
    #--------------------------------------------------------------------------------
    # Cluster number (total/of this kind):     1 /      1
        1
    # Number of vertices:
    1
    # Avg. Distance:
    0.000000
    # Vertices: (x,y,z) | d-vector label | s-vector
        0.000000        0.000000        0.000000         1       1
    #
    #--------------------------------------------------------------------------------
    # Cluster number (total/of this kind):     2 /      2
        2
    # Number of vertices:
    1
    # Avg. Distance:
    0.000000
    # Vertices: (x,y,z) | d-vector label | s-vector
        0.000000        0.000000        0.000000         1       2
    #
    """
    nCl = 2
    for iFig ∈ 1:length(clusters)
        println("iFig ",iFig)
        nV = size(clusters[iFig],2)
        thiskind = count(size.(clusters[1:iFig],2).==nV)
        total = nCl + iFig 
        fig = clusters[iFig] # Grab the vertex coordinates
        com = sum(fig,dims=2)/nV # Compute center of mass
        radius = √sum((fig.-com).^2)/nV
        str *= "#--------------------------------------------------------------------------------\n"
        str *= "# Cluster number (total/of this kind):   $total / $thiskind\n"
        str *= " $iFig\n"
        str *= "# Number of vertices:\n  $nV\n# Avg. Distance:\n"
        
        str *= "  $radius \n"
        str *= "# Vertices: (x,y,z) | d-vector label | s-vector \n"
        for iV ∈ 1:nV
            a,b,c = fig[:,iV]
            sVlabel = svectors[iFig][iV]
            str *= " $a    $b    $c    1  $sVlabel \n"
        end
        str *= "#\n"
    end

    open("fromStructures.clusters.out","w") do io 
        println(io,str)
    end
    return true
end

#using DelimitedFiles
m = readdlm("enum_PI_matrix.out")
rank(m)

# for >= pairs, rank = 6 (full rank)

# t1 = copy(finalCl[2])
# sv1 = copy(finalSV[2])
# rot
# t2 = transpose(transpose(t1)*rot)
# p = sortperm(collect(eachslice(t2,dims=2)))
# sort!(t2,dims=2)
# sv2 = sv1[p]
# isCommonShift(t1-t2) 
# t3 = copy(finalCl[3])
# sv3 = copy(finalSV[3])
# isCommonShift(t2-t3) && sv3==sv2