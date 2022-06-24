using GmshReader
using Test

using LinearAlgebra

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

@testset "fixsurface option" begin
    for filename in ("cube", "cube2")
        for fixsurface in (true, false)
            gmsh = readgmsh("$filename.msh"; fixsurface)
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
