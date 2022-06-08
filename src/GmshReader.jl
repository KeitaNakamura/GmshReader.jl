module GmshReader

using StaticArrays

import gmsh_jll
include(gmsh_jll.gmsh_api)

export
    readgmsh

struct ElementSet{dim}
    elementtags::Vector{Int}
    elementname::String
    order::Int
    numnodes::Int
    localnodecoord::Vector{SVector{dim, Float64}}
    numprimarynodes::Int
    connectivities::Vector{Vector{Int}}
end
num_elements(eltset::ElementSet) = length(eltset.connectivities)

function readgmsh_elements()
    dimtags::Vector{Tuple{Int, Int}} = gmsh.model.getPhysicalGroups()
    Dict{String, Vector{ElementSet}}(map(dimtags) do (dim, tag)
        tags = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        name = gmsh.model.getPhysicalName(dim, tag)
        group = map(tags) do tag′ # each PhysicalGroup having several entities
            elementtypes, elementtags::Vector{Vector{Int}}, nodetags_all::Vector{Vector{Int}} = gmsh.model.mesh.getElements(dim, tag′)
            sets = map(zip(elementtypes, nodetags_all, elementtags)) do (elttype, nodetags, elttags)
                elementname::String, _dim::Int, order::Int, numnodes::Int, localnodecoord::Vector{Float64}, numprimarynodes::Int = gmsh.model.mesh.getElementProperties(elttype) 
                @assert dim == _dim
                conns = [nodetags[i:i + (numnodes - 1)] for i in 1:numnodes:length(nodetags)]
                lcoord = SVector{dim, Float64}[(localnodecoord[i:i + (dim - 1)]) for i in 1:dim:length(localnodecoord)]
                ElementSet{dim}(elttags, elementname, order, numnodes, lcoord, numprimarynodes, conns)
            end
            reduce(vcat, sets) # all `dim`s are always the same in each group (maybe?)
        end
        name => group
    end)
end

function readgmsh_nodes()
    nodeid, nodes = gmsh.model.mesh.getNodes()
    dim::Int = gmsh.model.getDimension()
    [SVector{dim}(nodes[i:i + (dim - 1)]) for i in 1:3:length(nodes)]
end

function readgmsh(filename::String)
    @assert endswith(filename, ".msh")
    @assert isfile(filename)

    gmsh.initialize()
    gmsh.open(filename)

    gmsh.model.mesh.renumberNodes()
    gmsh.model.mesh.renumberElements()

    nodes = readgmsh_nodes()
    elements = readgmsh_elements()

    gmsh.finalize()
    Dict{String, Any}("nodes" => nodes, "elements" => elements)
end

end # module
