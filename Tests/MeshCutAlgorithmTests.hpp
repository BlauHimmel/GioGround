#pragma once

#include <gtest\gtest.h>

#include <MeshGenerator.hpp>
#include <HalfedgeMeshWrapper.hpp>
#include <Algorithm\MeshCutAlgorithm.hpp>

TEST(MeshCutAlgorithm, Function)
{
	GEO::Mesh Mesh;
	MeshGenerator::MeshGenHexagon(&Mesh);

	MeshAlgorithm::MeshCutAlgorithm MeshCut;
	GEO::index_t StartFacet = 0;
	MeshCut.PutArg(MeshAlgorithm::MeshCutAlgorithm::PARAMS_KEY_START_FACET, StartFacet);
	MeshCut.Execute(&Mesh);

	GEO::Attribute<bool> AttriIsCornerMarked;
	AttriIsCornerMarked.bind(Mesh.facet_corners.attributes(), "IsCornerMarked");

	ASSERT_TRUE(AttriIsCornerMarked[0]);
	ASSERT_TRUE(AttriIsCornerMarked[1]);
	ASSERT_TRUE(AttriIsCornerMarked[2]);
	ASSERT_TRUE(AttriIsCornerMarked[3]);
	ASSERT_FALSE(AttriIsCornerMarked[4]);
	ASSERT_FALSE(AttriIsCornerMarked[5]);
	ASSERT_TRUE(AttriIsCornerMarked[6]);
	ASSERT_TRUE(AttriIsCornerMarked[7]);
	ASSERT_TRUE(AttriIsCornerMarked[8]);
	ASSERT_FALSE(AttriIsCornerMarked[9]);
	ASSERT_TRUE(AttriIsCornerMarked[10]);
	ASSERT_FALSE(AttriIsCornerMarked[11]);
	ASSERT_FALSE(AttriIsCornerMarked[12]);
	ASSERT_TRUE(AttriIsCornerMarked[13]);
	ASSERT_TRUE(AttriIsCornerMarked[14]);
	ASSERT_TRUE(AttriIsCornerMarked[15]);
	ASSERT_FALSE(AttriIsCornerMarked[16]);
	ASSERT_TRUE(AttriIsCornerMarked[17]);
}