#pragma once

#include <gtest\gtest.h>

#include <MeshGenerator.hpp>
#include <HalfedgeMeshWrapper.hpp>
#include <Algorithm\MeshCutAlgorithm.hpp>

TEST(MeshCutAlgorithm, Function)
{
	GEO::Mesh Mesh;
	MeshGenerator::MeshGenHexagon(&Mesh);

	std::unique_ptr<MeshAlgorithm::MeshCutAlgorithm> MeshCut(new MeshAlgorithm::MeshCutAlgorithm());
	GEO::index_t StartFacet = 0;
	MeshCut->PutArg(MeshAlgorithm::MeshCutAlgorithm::PARAMS_KEY_START_FACET, StartFacet);
	MeshCut->Execute(&Mesh);

	GEO::index_t nBoundary = 0;
	for (GEO::index_t iCorner = 0; iCorner < Mesh.facet_corners.nb(); ++iCorner) 
	{
		if (Mesh.facet_corners.adjacent_facet(iCorner) == GEO::NO_FACET) 
		{
			nBoundary++;
		}
	}

	ASSERT_EQ(nBoundary, 1);
}