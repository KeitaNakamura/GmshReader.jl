using GmshReader
using Test

using LinearAlgebra
using Statistics

@testset "readgmsh" begin
    @test_throws Exception readgmsh("square")    # no extension
    @test_throws Exception readgmsh("square.ms") # wrong extension

    gmsh = (@inferred readgmsh("square.msh"))::GmshFile

    phygroups = gmsh.physicalgroups
    @testset "Nodes" begin
        for boundname in ("left", "right", "top", "bottom")
            set = phygroups[boundname].nodeset
            for (tag, x) in zip(set.nodetags, set.coord)
                i = only(findall(==(tag), phygroups["main"].nodeset.nodetags))
                @test phygroups["main"].nodeset.coord[i] == x
                @test gmsh.nodeset.coord[i] == x
            end
        end
        @test phygroups["main"].nodeset.coord == gmsh.nodeset.coord
    end

    @testset "Elements" begin
        # left
        for set in only(phygroups["left"].entities)
            @test set.elementname == "Line 2"
            @test length(set.connectivities) == 2
            @test gmsh.nodeset.coord[set.connectivities[1]] ≈ [[0.0,1.0,0.0], [0.0,0.5,0.0]]
            @test gmsh.nodeset.coord[set.connectivities[2]] ≈ [[0.0,0.5,0.0], [0.0,0.0,0.0]]
        end
        # right
        for set in only(phygroups["right"].entities)
            @test set.elementname == "Line 2"
            @test length(set.connectivities) == 2
            @test gmsh.nodeset.coord[set.connectivities[1]] ≈ [[1.0,0.0,0.0], [1.0,0.5,0.0]]
            @test gmsh.nodeset.coord[set.connectivities[2]] ≈ [[1.0,0.5,0.0], [1.0,1.0,0.0]]
        end
        # top
        for set in only(phygroups["top"].entities)
            @test set.elementname == "Line 2"
            @test length(set.connectivities) == 2
            @test gmsh.nodeset.coord[set.connectivities[1]] ≈ [[1.0,1.0,0.0], [0.5,1.0,0.0]]
            @test gmsh.nodeset.coord[set.connectivities[2]] ≈ [[0.5,1.0,0.0], [0.0,1.0,0.0]]
        end
        # bottom
        for set in only(phygroups["bottom"].entities)
            @test set.elementname == "Line 2"
            @test length(set.connectivities) == 2
            @test gmsh.nodeset.coord[set.connectivities[1]] ≈ [[0.0,0.0,0.0], [0.5,0.0,0.0]]
            @test gmsh.nodeset.coord[set.connectivities[2]] ≈ [[0.5,0.0,0.0], [1.0,0.0,0.0]]
        end
        # main
        for set in only(phygroups["main"].entities)
            @test set.elementtags == [9,10,11,12]
            @test set.elementname == "Quadrilateral 4"
            @test set.connectivities == [[3,7,9,6], [6,9,5,2], [7,4,8,9], [9,8,1,5]]
        end
    end
end

function compute_normal(solidfamily::String, facenodes::Vector{Vector{Float64}})
    if solidfamily == "Triangle" || solidfamily == "Quadrangle"
        v = facenodes[2] - facenodes[1]
        n = [v[2], -v[1]]
        return normalize(n)
    end
    if solidfamily == "Tetrahedron"
        v1 = facenodes[2] - facenodes[1]
        v2 = facenodes[3] - facenodes[1]
        return normalize(v1 × v2)
    end
    if solidfamily == "Hexahedron"
        v1 = facenodes[2] - facenodes[1]
        v2 = facenodes[4] - facenodes[1]
        return normalize(v1 × v2)
    end
    error()
end

@testset "fixsurface option" begin
    @testset "$(GmshReader.element_properties(familyname, args...).elementname)" for (familyname, args...) in (
            ("Triangle",    1),
            ("Triangle",    2),
            ("Quadrangle",  1),
            ("Quadrangle",  2, true),
            ("Quadrangle",  2),
            ("Tetrahedron", 1),
            ("Tetrahedron", 2),
            ("Hexahedron",  1),
            ("Hexahedron",  2, true),
            ("Hexahedron",  2),
        )
        prop = GmshReader.element_properties(familyname, args...)
        surface_list = GmshReader.SURFACE_LIST[prop.elementname]
        for conn in surface_list
            surface_nodes = prop.localnodecoord[conn]
            n = compute_normal(familyname, surface_nodes)
            xc = mean(prop.localnodecoord)
            xc_surface = mean(surface_nodes)
            @test (xc_surface - xc) ⋅ n > 0 # normal vector is regarded as outer direction when both vectors are the same direction
        end
    end
    @testset "readgmsh" begin
        for filename in ("cube", "cube2")
            for fixsurface in (true, false)
                gmsh = readgmsh(joinpath(@__DIR__, "$filename.msh"); fixsurface)
                coord = gmsh.nodeset.coord
                phygroups = gmsh.physicalgroups
                for (name, n) in (("left", [-1,0,0]), ("right", [1,0,0]), ("bottom", [0,-1,0]), ("top", [0,1,0]), ("back", [0,0,-1]), ("front", [0,0,1]))
                    for elementset in only(phygroups[name].entities)
                        for conn in elementset.connectivities
                            if filename == "cube"
                                @assert length(conn) == 3
                            elseif filename == "cube2"
                                @assert length(conn) == 6
                            else
                                error("unreachable")
                            end
                            x1 = coord[conn[2]] - coord[conn[1]]
                            x2 = coord[conn[end]] - coord[conn[1]]
                            if !fixsurface && name == "back"
                                @test normalize(x1 × x2) ≈ -n
                            else
                                @test normalize(x1 × x2) ≈ n
                            end
                        end
                    end
                end
            end
        end
    end
end

@testset "Utilities" begin
    @testset "element_properties" begin
        # quad 4
        prop = GmshReader.element_properties("Quadrangle", 1)
        @test prop.elementname == "Quadrilateral 4"
        @test prop.dim == 2
        @test prop.order == 1
        @test prop.numnodes == 4
        @test prop.localnodecoord == [[-1.0,-1.0], [1.0,-1.0], [1.0,1.0], [-1.0,1.0]]
        @test prop.numprimarynodes == 4
        # quad 8
        prop = GmshReader.element_properties("Quadrangle", 2, true)
        @test prop.elementname == "Quadrilateral 8"
        @test prop.dim == 2
        @test prop.order == 2
        @test prop.numnodes == 8
        @test prop.localnodecoord == [[-1.0,-1.0], [1.0,-1.0], [1.0,1.0], [-1.0,1.0], [0.0,-1.0], [1.0,0.0], [0.0,1.0], [-1.0,0.0]]
        @test prop.numprimarynodes == 4
        # quad 9
        prop = GmshReader.element_properties("Quadrangle", 2)
        @test prop.elementname == "Quadrilateral 9"
        @test prop.dim == 2
        @test prop.order == 2
        @test prop.numnodes == 9
        @test prop.localnodecoord == [[-1.0,-1.0], [1.0,-1.0], [1.0,1.0], [-1.0,1.0], [0.0,-1.0], [1.0,0.0], [0.0,1.0], [-1.0,0.0], [0.0,0.0]]
        @test prop.numprimarynodes == 4
    end
end
