module GmshReader

import gmsh_jll
include(gmsh_jll.gmsh_api)

export
    readgmsh,
    NodeSet,
    ElementSet,
    PhysicalGroup,
    GmshFile

struct NodeSet
    nodetags::Vector{Int}
    dim::Int
    coord::Vector{Vector{Float64}}
end

struct ElementSet
    elementname::String
    elementtags::Vector{Int}
    dim::Int
    order::Int
    numnodes::Int
    localnodecoord::Vector{Vector{Float64}}
    numprimarynodes::Int
    connectivities::Vector{Vector{Int}}
end

struct Entity <: AbstractVector{ElementSet}
    data::Vector{ElementSet}
end
Base.size(e::Entity) = size(e.data)
Base.getindex(e::Entity, i::Int) = e.data[i]
Base.convert(::Type{Entity}, data::Vector{ElementSet}) = Entity(data)
Base.summary(io::IO, e::Entity) = print(io, string(Entity, " vector with ", length(e), " element ", ifelse(length(e)==1, "type", "types")))

struct PhysicalGroup
    nodeset::NodeSet
    entities::Vector{Entity}
end
function Base.show(io::IO, ::PhysicalGroup)
    print(io, "PhysicalGroup(nodeset, entities)")
end

struct GmshFile
    name::String
    nodeset::NodeSet
    physicalgroups::Dict{String, PhysicalGroup}
end
function Base.show(io::IO, file::GmshFile)
    print(io, "GmshFile(\"", file.name, "\", nodeset, physicalgroups)")
end

function collectwithstep(x::AbstractVector, step::Int)
    if step == 0
        [[only(x)]]
    else
        [x[i:i+step-1] for i in 1:step:length(x)]
    end
end

function readgmsh_physicalgroup_elements()
    dimtags::Vector{Tuple{Int, Int}} = gmsh.model.getPhysicalGroups()
    Dict{String, Vector{Entity}}(map(dimtags) do (dim, tag)
        # loop over PhysicalGroups
        # PhysicalGroup have several entities
        # all entities always have the same dimension (maybe?)
        tags = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        name = gmsh.model.getPhysicalName(dim, tag)
        group = map(tags) do tag′ # each entity
            elementtypes, elementtags::Vector{Vector{Int}}, nodetags_all::Vector{Vector{Int}} = gmsh.model.mesh.getElements(dim, tag′)
            # elements in an entity are grouped into element types
            # and all types should be unique (maybe...)
            # so let's make a Dict having elementnames as `keys`
            map(zip(elementtypes, nodetags_all, elementtags)) do (elttype, nodetags, elttags)
                elementname::String, _dim::Int, order::Int, numnodes::Int, localnodecoord::Vector{Float64}, numprimarynodes::Int = gmsh.model.mesh.getElementProperties(elttype) 
                @assert dim == _dim
                conns = collectwithstep(nodetags, numnodes)
                lcoord = collectwithstep(localnodecoord, dim)
                ElementSet(elementname, elttags, dim, order, numnodes, lcoord, numprimarynodes, conns)
            end
        end
        name => group
    end)
end

function readgmsh_physicalgroup_nodes()
    dimtags::Vector{Tuple{Int, Int}} = gmsh.model.getPhysicalGroups()
    Dict{String, NodeSet}(map(dimtags) do (dim, tag)
        name = gmsh.model.getPhysicalName(dim, tag)
        nodetags::Vector{Int}, coord::Vector{Float64} = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)
        nodeset = NodeSet(nodetags, dim, collectwithstep(coord, 3))
        name => nodeset
    end)
end

function readgmsh_nodeset()
    nodetags::Vector{Int}, coord::Vector{Float64} = gmsh.model.mesh.getNodes()
    dim::Int = gmsh.model.getDimension()
    NodeSet(nodetags, dim, collectwithstep(coord, 3))
end

function readgmsh(filename::String; fixsurface::Bool = false)
    @assert endswith(filename, ".msh")
    @assert isfile(filename)

    gmsh.initialize()
    gmsh.open(filename)

    gmsh.model.mesh.renumberNodes()
    gmsh.model.mesh.renumberElements()

    nodeset = readgmsh_nodeset()
    physicalgroups = Dict(map(readgmsh_physicalgroup_nodes(), readgmsh_physicalgroup_elements()) do (key1, val1), (key2, val2)
        @assert key1 == key2
        key1 => PhysicalGroup(val1, val2)
    end)
    file = GmshFile(filename, nodeset, physicalgroups)

    fixsurface && fix_surfacemesh!(file)

    gmsh.finalize()
    file
end

#####################
# fixsurface option #
#####################

SURFACE_LIST = Dict{String, Vector{Vector{Int}}}(
    "Line 2" => [[1], [2]],
    "Line 3" => [[1], [2]],
    "Triangle 3" => [[1,2], [2,3], [3,1]],
    "Triangle 6" => [[1,4,2], [2,5,3], [3,6,1]],
    "Quadrilateral 4" => [[1,2], [2,3], [3,4], [4,1]],
    "Quadrilateral 9" => [[1,5,2], [2,6,3], [3,7,4], [4,8,1]],
    "Tetrahedron 4"  => [[1,4,3], [4,2,3], [2,1,3], [1,2,4]],
    "Tetrahedron 10" => [[1,8,4,9,3,7], [4,10,2,6,3,9], [2,5,1,7,3,6], [1,5,2,10,4,8]],
    "Hexahedron 8"  => [[5,6,7,8], [2,1,4,3], [1,5,8,4], [6,2,3,7], [1,2,6,5], [8,7,3,4]],
    "Hexahedron 20" => [[5,17,6,19,7,20,8,18], [2,9,1,10,4,14,3,12], [1,11,5,18,8,16,4,10], [6,13,2,12,3,15,7,19], [1,9,2,13,6,17,5,11], [8,20,7,15,3,14,4,16]],
)

function fix_surfacemesh!(file::GmshFile)
    # grouping with dimension
    dim_elementset = Dict{Int, Vector{ElementSet}}()
    for (name, physicalgroup) in file.physicalgroups
        for entity in physicalgroup.entities
            for elementset in entity
                list = get!(dim_elementset, elementset.dim, ElementSet[])
                push!(list, elementset)
            end
        end
    end
    # solid elementset
    dim = maximum(keys(dim_elementset))
    solid_elementset_vector = dim_elementset[dim]
    # loop over surface elementset
    for (dim′, elementset_vector) in dim_elementset
        dim′ == dim && continue
        for elementset in elementset_vector
            for conn in elementset.connectivities # loop over elements
                fix_surfacemesh!(conn, solid_elementset_vector)
            end
        end
    end
end

function fix_surfacemesh!(surface_conn::Vector{Int}, solid_elementset_vector::Vector{ElementSet})
    for elementset in solid_elementset_vector
        solid_name = elementset.elementname
        surface_list = SURFACE_LIST[solid_name]
        for solid_conn in elementset.connectivities
            for inds in surface_list
                if Set(solid_conn[inds]) == Set(surface_conn)
                    @. surface_conn = solid_conn[inds]
                    return surface_conn
                end
            end
        end
    end
end

end # module
