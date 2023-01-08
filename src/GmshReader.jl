module GmshReader

import gmsh_jll
include(gmsh_jll.gmsh_api)

const IS_ACTIVE = Ref(false)
gmsh_initialize() = !IS_ACTIVE[] && (gmsh.initialize(); IS_ACTIVE[]=true)
gmsh_finalize() = IS_ACTIVE[] && (gmsh.finalize(); IS_ACTIVE[]=false)

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

"""
# GmshFile structure

```
GmshFile
├── <name> :: String
├── <nodeset> :: NodeSet
└── <physicalgroups> :: Dict{String, PhysicalGroup}
    │
    ├── key1 => value1
    │           ├── <nodeset> :: NodeSet
    │           └── <entities> :: Vector{Entity}
    │               │
    │               ├── entity1 :: Entity
    │               │   ├── elementset1 :: ElementSet
    │               │   └── elementset2 :: ElementSet
    │               │
    │               └── entity2 :: Entity
    │                   ├ ...
    │
    └── key2 => value2
                ├ ...
```
"""
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

function readgmsh_nodeset()
    nodetags::Vector{Int}, coord::Vector{Float64} = gmsh.model.mesh.getNodes()
    dim::Int = gmsh.model.getDimension()
    NodeSet(nodetags, dim, collectwithstep(coord, 3))
end

function readgmsh_physicalgroups()
    dimtags::Vector{Tuple{Int, Int}} = gmsh.model.getPhysicalGroups()
    Dict{String, PhysicalGroup}(map(dimtags) do (dim, tag) # loop over PhysicalGroups
        # nodes
        nodetags::Vector{Int}, coord::Vector{Float64} = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)
        nodeset = NodeSet(nodetags, dim, collectwithstep(coord, 3))

        # PhysicalGroup have entities
        # all entities always have the same dimension (maybe?)
        tags = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        entities = map(tags) do tag′ # each entity
            elementtypes, elementtags::Vector{Vector{Int}}, nodetags_all::Vector{Vector{Int}} = gmsh.model.mesh.getElements(dim, tag′)
            # elements in an entity are grouped into element types
            # and all types should be unique (maybe...)
            map(zip(elementtypes, nodetags_all, elementtags)) do (elttype, nodetags, elttags)
                elementname::String, _dim::Int, order::Int, numnodes::Int, localnodecoord::Vector{Float64}, numprimarynodes::Int = gmsh.model.mesh.getElementProperties(elttype) 
                @assert dim == _dim
                conns = collectwithstep(nodetags, numnodes)
                lcoord = collectwithstep(localnodecoord, dim)
                ElementSet(elementname, elttags, dim, order, numnodes, lcoord, numprimarynodes, conns)
            end
        end

        name = gmsh.model.getPhysicalName(dim, tag)
        name => PhysicalGroup(nodeset, entities)
    end)
end

function readgmsh(filename::String; fixsurface::Bool = false)
    @assert endswith(filename, ".msh")
    @assert isfile(filename)

    gmsh_initialize()
    gmsh.open(filename)

    gmsh.model.mesh.renumberNodes()
    gmsh.model.mesh.renumberElements()

    nodeset = readgmsh_nodeset()
    physicalgroups = readgmsh_physicalgroups()
    file = GmshFile(filename, nodeset, physicalgroups)

    fixsurface && fix_surfacemesh!(file)

    gmsh_finalize()
    file
end

#####################
# fixsurface option #
#####################

const FACE_LIST = Dict{String, Vector{Vector{Int}}}(
    # Line
    "Line 2" => [[1], [2]],
    "Line 3" => [[1], [2]],
    # Triangle
    "Triangle 3" => [[1,2],   [2,3],   [3,1]  ],
    "Triangle 6" => [[1,2,4], [2,3,5], [3,1,6]],
    # Quadrangle
    "Quadrilateral 4" => [[1,2],   [2,3],   [3,4],   [4,1]  ],
    "Quadrilateral 8" => [[1,2,5], [2,3,6], [3,4,7], [4,1,8]],
    "Quadrilateral 9" => [[1,2,5], [2,3,6], [3,4,7], [4,1,8]],
    # Tetrahedron
    "Tetrahedron 4"  => [[1,4,3],       [4,2,3],        [2,1,3],       [1,2,4]       ],
    "Tetrahedron 10" => [[1,4,3,8,9,7], [4,2,3,10,6,9], [2,1,3,5,7,6], [1,2,4,5,10,8]],
    # Hexahedron
    "Hexahedron 8"  => [[5,6,7,8],                [2,1,4,3],               [1,5,8,4],                [6,2,3,7],                [1,2,6,5],               [8,7,3,4]               ],
    "Hexahedron 20" => [[5,6,7,8,17,19,20,18],    [2,1,4,3,9,10,14,12],    [1,5,8,4,11,18,16,10],    [6,2,3,7,13,12,15,19],    [1,2,6,5,9,13,17,11],    [8,7,3,4,20,15,14,16]   ],
    "Hexahedron 27" => [[5,6,7,8,17,19,20,18,26], [2,1,4,3,9,10,14,12,21], [1,5,8,4,11,18,16,10,23], [6,2,3,7,13,12,15,19,24], [1,2,6,5,9,13,17,11,22], [8,7,3,4,20,15,14,16,25]],
)

function fix_surfacemesh!(file::GmshFile)
    # grouping with dimension
    dim2elementsets = Dict{Int, Vector{ElementSet}}()
    for (name, physicalgroup) in file.physicalgroups
        for entity in physicalgroup.entities
            for elementset in entity
                list = get!(dim2elementsets, elementset.dim, ElementSet[])
                push!(list, elementset)
            end
        end
    end

    # solid elementset
    dim = maximum(keys(dim2elementsets))
    solid_eltsets = dim2elementsets[dim]
    face_eltsets_vec = [eltsets for (d,eltsets) in dim2elementsets if d != dim]

    # loop over surface elementset and sort! all connectivities
    for face_eltsets in face_eltsets_vec, elementset in face_eltsets
        foreach(sort!, elementset.connectivities)
    end

    # loop over all solid elements
    @inbounds for solid_eltset in solid_eltsets
        solid_name = solid_eltset.elementname
        facelist = FACE_LIST[solid_name]
        tmp = similar(first(facelist))
        for conn in solid_eltset.connectivities
            for inds in facelist
                tmp .= @view conn[inds]
                sort!(tmp)
                # check all surface elements
                for face_eltsets in face_eltsets_vec, face_eltset in face_eltsets, face_conn in face_eltset.connectivities
                    if tmp == face_conn
                        face_conn .= @view conn[inds]
                        break
                    end
                end
            end
        end
    end
end

#############
# Utilities #
#############

"""
    element_properties(familyname, order, serendip = false)

Return element properties given its family name `familyname` ("Point", "Line", "Triangle",
"Quadrangle", "Tetrahedron", "Pyramid", "Prism", "Hexahedron") and polynomial `order`.
If `serendip` is `true`, return the corresponding serendip element type (element without interior nodes).
"""
function element_properties(familyname::String, order::Int, serendip::Bool = false)
    gmsh_initialize()
    elementtype = gmsh.model.mesh.getElementType(familyname, order, serendip)
    elementname::String, dim::Int, order::Int, numnodes::Int, localnodecoord′::Vector{Float64}, numprimarynodes::Int = gmsh.model.mesh.getElementProperties(elementtype)
    gmsh_finalize()
    localnodecoord = collectwithstep(localnodecoord′, dim)
    (; elementname, dim, order, numnodes, localnodecoord, numprimarynodes)
end

end # module
