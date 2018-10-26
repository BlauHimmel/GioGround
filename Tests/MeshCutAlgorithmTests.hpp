#pragma once

#include <gtest\gtest.h>
#include <geogram\mesh\mesh_io.h>

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
		std::unique_ptr<MeshAlgorithm::MeshCutAlgorithm> MeshCut(new MeshAlgorithm::MeshCutAlgorithm());
		GEO::index_t StartFacet = 0;
		MeshCut->PutArg(MeshAlgorithm::MeshCutAlgorithm::PARAMS_KEY_START_FACET, StartFacet);
		MeshCut->Execute(&Mesh);

		HalfedgeMeshWrapper Wrapper(&Mesh);

		GEO::index_t iBeginCorner = GEO::NO_CORNER;
		for (GEO::index_t iCorner = 0; iCorner < Mesh.facet_corners.nb(); ++iCorner)
		{
			if (Mesh.facet_corners.adjacent_facet(iCorner) == GEO::NO_FACET)
			{
				iBeginCorner = iCorner;
				break;
			}
		}

		GEO::index_t iCorner = iBeginCorner;
		GEO::index_t iIter = 0;
		GEO::index_t nMaxIterationTime = Mesh.facets.nb() * 3;

		do
		{
			iCorner = Wrapper.Next(iCorner);
			do
			{
				iCorner = Wrapper.Corner2Corner[iCorner];
			} while (Mesh.facet_corners.adjacent_facet(iCorner) != GEO::NO_FACET);
			iIter++;
		} while (iCorner != iBeginCorner && iIter < nMaxIterationTime);
		ASSERT_EQ(iCorner, iBeginCorner);

		Mesh.clear(false, false);
	}
}