module GmshReader

import gmsh_jll
include(gmsh_jll.gmsh_api)

export
    readgmsh

struct NodeSet
    nodetags::Vector{Int}
    dim::Int
    coord::Vector{Vector{Float64}}
end

struct ElementSet
    elementtags::Vector{Int}
    elementname::String
    dim::Int
    order::Int
    numnodes::Int
    localnodecoord::Vector{Vector{Float64}}
    numprimarynodes::Int
    connectivities::Vector{Vector{Int}}
end

struct GmshFile
    nodeset::Dict{String, NodeSet}
    elementset::Dict{String, Vector{ElementSet}}
end

function Base.show(io::IO, mime::MIME"text/plain", gmsh::GmshFile)
    println(io, summary(gmsh), ":")
    show(io, mime, gmsh.nodeset)
    println(io)
    show(io, mime, gmsh.elementset)
end

function collectwithstep(x::AbstractVector, step::Int)
    if step == 0
        [[only(x)]]
    else
        [x[i:i+step-1] for i in 1:step:length(x)]
    end
end

function readgmsh_elementset()
    dimtags::Vector{Tuple{Int, Int}} = gmsh.model.getPhysicalGroups()
    Dict{String, Vector{ElementSet}}(map(dimtags) do (dim, tag)
        tags = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        name = gmsh.model.getPhysicalName(dim, tag)
        group = map(tags) do tag′ # each PhysicalGroup having several entities
            elementtypes, elementtags::Vector{Vector{Int}}, nodetags_all::Vector{Vector{Int}} = gmsh.model.mesh.getElements(dim, tag′)
            sets = map(zip(elementtypes, nodetags_all, elementtags)) do (elttype, nodetags, elttags)
                elementname::String, _dim::Int, order::Int, numnodes::Int, localnodecoord::Vector{Float64}, numprimarynodes::Int = gmsh.model.mesh.getElementProperties(elttype) 
                @assert dim == _dim
                conns = collectwithstep(nodetags, numnodes)
                lcoord = collectwithstep(localnodecoord, dim)
                ElementSet(elttags, elementname, dim, order, numnodes, lcoord, numprimarynodes, conns)
            end
            reduce(vcat, sets) # all `dim`s are always the same in each group (maybe?)
        end
        name => group
    end)
end

function readgmsh_nodeset()
    dimtags::Vector{Tuple{Int, Int}} = gmsh.model.getPhysicalGroups()
    Dict{String, NodeSet}(map(dimtags) do (dim, tag)
        name = gmsh.model.getPhysicalName(dim, tag)
        nodetags::Vector{Int}, coord::Vector{Float64} = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)
        nodeset = NodeSet(nodetags, dim, collectwithstep(coord, 3))
        name => nodeset
    end)
end

function readgmsh(filename::String)
    @assert endswith(filename, ".msh")
    @assert isfile(filename)

    gmsh.initialize()
    gmsh.open(filename)

    gmsh.model.mesh.renumberNodes()
    gmsh.model.mesh.renumberElements()

    nodeset = readgmsh_nodeset()
    elementset = readgmsh_elementset()

    gmsh.finalize()
    GmshFile(nodeset, elementset)
end

end # module
