# GmshReader

[![Build Status](https://github.com/KeitaNakamura/GmshReader.jl/workflows/CI/badge.svg)](https://github.com/KeitaNakamura/GmshReader.jl/actions)
[![codecov](https://codecov.io/gh/KeitaNakamura/GmshReader.jl/branch/main/graph/badge.svg?token=KCLNS7IKOM)](https://codecov.io/gh/KeitaNakamura/GmshReader.jl)

## Installation

```julia
pkg> add https://github.com/KeitaNakamura/GmshReader.jl.git
```

## Usage

`readgmsh` reads `.msh` file and return `GmshReader.GmshFile`.

```julia
julia> readgmsh("test/cube.msh")
Info    : Reading 'test/cube.msh'...
Info    : 27 entities
Info    : 14 nodes
Info    : 48 elements
Info    : Done reading 'test/cube.msh'
GmshFile("test/cube.msh", nodeset, physicalgroups)
```

`GmshReader.GmshFile` structure is as follows:

```
GmshReader.GmshFile
├── <name> :: String
├── <nodeset> :: GmshReader.NodeSet
└── <physicalgroups> :: Dict{String, GmshReader.PhysicalGroup}
    │
    ├── key1 => value1
    │           ├── <nodeset> :: GmshReader.NodeSet
    │           └── <entities> :: Vector{GmshReader.Entity}
    │               │
    │               ├── entity1 :: GmshReader.Entity
    │               │   ├── elementset1 :: GmshReader.ElementSet
    │               │   └── elementset2 :: GmshReader.ElementSet
    │               │
    │               └── entity2 :: GmshReader.Entity
    │                   ├ ...
    │
    └── key2 => value2
                ├ ...
```
