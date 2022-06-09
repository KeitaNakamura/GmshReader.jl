using GmshReader
using Test

@testset "readgmsh" begin
    @test_throws Exception readgmsh("square")    # no extension
    @test_throws Exception readgmsh("square.ms") # wrong extension

    gmsh = (@inferred readgmsh("square.msh"))::GmshReader.GmshFile
    nodeset = gmsh.nodeset
    elementset = gmsh.elementset

    @testset "NodeSet" begin
        for boundname in ("left", "right", "top", "bottom")
            set = nodeset[boundname]
            for (tag, x) in zip(set.nodetags, set.coord)
                i = only(findall(==(tag), nodeset["main"].nodetags))
                @test nodeset["main"].coord[i] == x
            end
        end
    end

    @testset "ElementSet" begin
        # left
        for set in elementset["left"]
            @test set.elementname == "Line 2"
            @test length(set.connectivities) == 2
            @test nodeset["main"].coord[set.connectivities[1]] ≈ [[0.0,1.0,0.0], [0.0,0.5,0.0]]
            @test nodeset["main"].coord[set.connectivities[2]] ≈ [[0.0,0.5,0.0], [0.0,0.0,0.0]]
        end
        # right
        for set in elementset["right"]
            @test set.elementname == "Line 2"
            @test length(set.connectivities) == 2
            @test nodeset["main"].coord[set.connectivities[1]] ≈ [[1.0,0.0,0.0], [1.0,0.5,0.0]]
            @test nodeset["main"].coord[set.connectivities[2]] ≈ [[1.0,0.5,0.0], [1.0,1.0,0.0]]
        end
        # top
        for set in elementset["top"]
            @test set.elementname == "Line 2"
            @test length(set.connectivities) == 2
            @test nodeset["main"].coord[set.connectivities[1]] ≈ [[1.0,1.0,0.0], [0.5,1.0,0.0]]
            @test nodeset["main"].coord[set.connectivities[2]] ≈ [[0.5,1.0,0.0], [0.0,1.0,0.0]]
        end
        # bottom
        for set in elementset["bottom"]
            @test set.elementname == "Line 2"
            @test length(set.connectivities) == 2
            @test nodeset["main"].coord[set.connectivities[1]] ≈ [[0.0,0.0,0.0], [0.5,0.0,0.0]]
            @test nodeset["main"].coord[set.connectivities[2]] ≈ [[0.5,0.0,0.0], [1.0,0.0,0.0]]
        end
        # main
        for set in elementset["main"]
            @test set.elementtags == [9,10,11,12]
            @test set.elementname == "Quadrilateral 4"
            @test set.connectivities == [[3,7,9,6], [6,9,5,2], [7,4,8,9], [9,8,1,5]]
        end
    end
end
