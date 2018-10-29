#pragma once

#include <gtest\gtest.h>
#include <geogram\mesh\mesh_io.h>
#include <geogram\mesh\mesh_reorder.h>

#include <MeshGenerator.hpp>
#include <HalfedgeMeshWrapper.hpp>
#include <Algorithm\MeshCutAlgorithm.hpp>

TEST(MeshCutAlgorithm, Function)
{
	std::string Root = "..\\Mesh\\";
	std::vector<std::string> BenchMarkModels = 
	{
		"plane_hole1.obj",
		"plane_hole2.obj",
		"donut_genus1.obj",
		"simple_donut_genus1.obj",
		"sphere_hole1.obj",
		"sphere_hole2.obj",
		"kitten.obj",
		"kitten_holes.obj",
		"king.obj"
	};

	for (std::string & Filename : BenchMarkModels)
	{
		GEO::Mesh Mesh;
		GEO::mesh_load(Root + Filename, Mesh);
		Mesh.vertices.set_double_precision();

		GEO::mesh_reorder(Mesh, GEO::MESH_ORDER_MORTON);
		
		std::unique_ptr<MeshAlgorithm::IMeshAlgorithm> MeshCut(new MeshAlgorithm::MeshCutAlgorithm());
		GEO::index_t StartFacet = 0;
		MeshCut->PutArg(MeshAlgorithm::MeshCutAlgorithm::PARAMS_KEY_START_FACET, StartFacet);
		MeshCut->Execute(&Mesh);
		ASSERT_EQ(MeshCut->GetBoundaryNumber(&Mesh), 1);
		Mesh.clear(false, false);
	}
}