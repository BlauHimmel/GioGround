#pragma once

#include <gtest\gtest.h>

#include <MeshGenerator.hpp>
#include <HalfedgeMeshWrapper.hpp>

TEST(HalfedgeMeshWrapper, Function_Next)
{
	GEO::Mesh Mesh;
	ASSERT_TRUE(MeshGenerator::MeshGenHexagon(&Mesh));
	HalfedgeMeshWrapper Halfedge(&Mesh);

	ASSERT_EQ(Halfedge.Next(0), 1);
	ASSERT_EQ(Halfedge.Next(1), 2);
	ASSERT_EQ(Halfedge.Next(2), 3);
	ASSERT_EQ(Halfedge.Next(3), 0);
	
	ASSERT_EQ(Halfedge.Next(4), 5);
	ASSERT_EQ(Halfedge.Next(5), 6);
	ASSERT_EQ(Halfedge.Next(6), 7);
	ASSERT_EQ(Halfedge.Next(7), 4);

	ASSERT_EQ(Halfedge.Next(8), 9);
	ASSERT_EQ(Halfedge.Next(9), 10);
	ASSERT_EQ(Halfedge.Next(10), 8);

	ASSERT_EQ(Halfedge.Next(11), 12);
	ASSERT_EQ(Halfedge.Next(12), 13);
	ASSERT_EQ(Halfedge.Next(13), 14);
	ASSERT_EQ(Halfedge.Next(14), 11);

	ASSERT_EQ(Halfedge.Next(15), 16);
	ASSERT_EQ(Halfedge.Next(16), 17);
	ASSERT_EQ(Halfedge.Next(17), 15);
}

TEST(HalfedgeMeshWrapper, Function_Prev)
{
	GEO::Mesh Mesh;
	ASSERT_TRUE(MeshGenerator::MeshGenHexagon(&Mesh));
	HalfedgeMeshWrapper Halfedge(&Mesh);

	ASSERT_EQ(Halfedge.Prev(0), 3);
	ASSERT_EQ(Halfedge.Prev(1), 0);
	ASSERT_EQ(Halfedge.Prev(2), 1);
	ASSERT_EQ(Halfedge.Prev(3), 2);

	ASSERT_EQ(Halfedge.Prev(4), 7);
	ASSERT_EQ(Halfedge.Prev(5), 4);
	ASSERT_EQ(Halfedge.Prev(6), 5);
	ASSERT_EQ(Halfedge.Prev(7), 6);

	ASSERT_EQ(Halfedge.Prev(8), 10);
	ASSERT_EQ(Halfedge.Prev(9), 8);
	ASSERT_EQ(Halfedge.Prev(10), 9);

	ASSERT_EQ(Halfedge.Prev(11), 14);
	ASSERT_EQ(Halfedge.Prev(12), 11);
	ASSERT_EQ(Halfedge.Prev(13), 12);
	ASSERT_EQ(Halfedge.Prev(14), 13);

	ASSERT_EQ(Halfedge.Prev(15), 17);
	ASSERT_EQ(Halfedge.Prev(16), 15);
	ASSERT_EQ(Halfedge.Prev(17), 16);
}

TEST(HalfedgeMeshWrapper, Function_Opposite)
{
	GEO::Mesh Mesh;
	ASSERT_TRUE(MeshGenerator::MeshGenHexagon(&Mesh));
	HalfedgeMeshWrapper Halfedge(&Mesh);

	ASSERT_EQ(Halfedge.Opposite(0), 10);
	ASSERT_EQ(Halfedge.Opposite(1), 14);
	ASSERT_EQ(Halfedge.Opposite(2), 17);
	ASSERT_EQ(Halfedge.Opposite(3), 7);

	ASSERT_EQ(Halfedge.Opposite(4), GEO::NO_CORNER);
	ASSERT_EQ(Halfedge.Opposite(5), GEO::NO_CORNER);
	ASSERT_EQ(Halfedge.Opposite(6), 8);
	ASSERT_EQ(Halfedge.Opposite(7), 3);

	ASSERT_EQ(Halfedge.Opposite(8), 6);
	ASSERT_EQ(Halfedge.Opposite(9), GEO::NO_CORNER);
	ASSERT_EQ(Halfedge.Opposite(10), 0);

	ASSERT_EQ(Halfedge.Opposite(11), GEO::NO_CORNER);
	ASSERT_EQ(Halfedge.Opposite(12), GEO::NO_CORNER);
	ASSERT_EQ(Halfedge.Opposite(13), 15);
	ASSERT_EQ(Halfedge.Opposite(14), 1);

	ASSERT_EQ(Halfedge.Opposite(15), 13);
	ASSERT_EQ(Halfedge.Opposite(16), GEO::NO_CORNER);
	ASSERT_EQ(Halfedge.Opposite(17), 2);
}

TEST(HalfedgeMeshWrapper, Function_Facet)
{
	GEO::Mesh Mesh;
	ASSERT_TRUE(MeshGenerator::MeshGenHexagon(&Mesh));
	HalfedgeMeshWrapper Halfedge(&Mesh);

	ASSERT_EQ(Halfedge.Facet(0), 0);
	ASSERT_EQ(Halfedge.Facet(1), 0);
	ASSERT_EQ(Halfedge.Facet(2), 0);
	ASSERT_EQ(Halfedge.Facet(3), 0);

	ASSERT_EQ(Halfedge.Facet(4), 1);
	ASSERT_EQ(Halfedge.Facet(5), 1);
	ASSERT_EQ(Halfedge.Facet(6), 1);
	ASSERT_EQ(Halfedge.Facet(7), 1);

	ASSERT_EQ(Halfedge.Facet(8), 2);
	ASSERT_EQ(Halfedge.Facet(9), 2);
	ASSERT_EQ(Halfedge.Facet(10), 2);

	ASSERT_EQ(Halfedge.Facet(11), 3);
	ASSERT_EQ(Halfedge.Facet(12), 3);
	ASSERT_EQ(Halfedge.Facet(13), 3);
	ASSERT_EQ(Halfedge.Facet(14), 3);

	ASSERT_EQ(Halfedge.Facet(15), 4);
	ASSERT_EQ(Halfedge.Facet(16), 4);
	ASSERT_EQ(Halfedge.Facet(17), 4);
}

TEST(HalfedgeMeshWrapper, Function_NextAroundVertex)
{
	GEO::Mesh Mesh;
	ASSERT_TRUE(MeshGenerator::MeshGenHexagon(&Mesh));
	HalfedgeMeshWrapper Halfedge(&Mesh);

	ASSERT_EQ(Halfedge.NextAroundVertex(10), GEO::NO_CORNER);
	ASSERT_EQ(Halfedge.NextAroundVertex(1), 10);
	ASSERT_EQ(Halfedge.NextAroundVertex(11), 1);

	ASSERT_EQ(Halfedge.NextAroundVertex(12), GEO::NO_CORNER);

	ASSERT_EQ(Halfedge.NextAroundVertex(13), GEO::NO_CORNER);
	ASSERT_EQ(Halfedge.NextAroundVertex(16), 13);

	ASSERT_EQ(Halfedge.NextAroundVertex(17), GEO::NO_CORNER);
	ASSERT_EQ(Halfedge.NextAroundVertex(3), 17);
	ASSERT_EQ(Halfedge.NextAroundVertex(4), 3);

	ASSERT_EQ(Halfedge.NextAroundVertex(5), GEO::NO_CORNER);

	ASSERT_EQ(Halfedge.NextAroundVertex(6), GEO::NO_CORNER);
	ASSERT_EQ(Halfedge.NextAroundVertex(9), 6);

	ASSERT_EQ(Halfedge.NextAroundVertex(0), 7);
	ASSERT_EQ(Halfedge.NextAroundVertex(7), 8);
	ASSERT_EQ(Halfedge.NextAroundVertex(8), 0);

	ASSERT_EQ(Halfedge.NextAroundVertex(2), 14);
	ASSERT_EQ(Halfedge.NextAroundVertex(15), 2);
	ASSERT_EQ(Halfedge.NextAroundVertex(14), 15);
}

TEST(HalfedgeMeshWrapper, Function_OriginVertex)
{
	GEO::Mesh Mesh;
	ASSERT_TRUE(MeshGenerator::MeshGenHexagon(&Mesh));
	HalfedgeMeshWrapper Halfedge(&Mesh);

	ASSERT_EQ(Halfedge.Origin(10), 0);
	ASSERT_EQ(Halfedge.Origin(1), 0);
	ASSERT_EQ(Halfedge.Origin(11), 0);

	ASSERT_EQ(Halfedge.Origin(12), 1);

	ASSERT_EQ(Halfedge.Origin(13), 2);
	ASSERT_EQ(Halfedge.Origin(16), 2);

	ASSERT_EQ(Halfedge.Origin(17), 3);
	ASSERT_EQ(Halfedge.Origin(3), 3);
	ASSERT_EQ(Halfedge.Origin(4), 3);

	ASSERT_EQ(Halfedge.Origin(5), 4);

	ASSERT_EQ(Halfedge.Origin(6), 5);
	ASSERT_EQ(Halfedge.Origin(9), 5);

	ASSERT_EQ(Halfedge.Origin(0), 6);
	ASSERT_EQ(Halfedge.Origin(7), 6);
	ASSERT_EQ(Halfedge.Origin(8), 6);

	ASSERT_EQ(Halfedge.Origin(2), 7);
	ASSERT_EQ(Halfedge.Origin(15), 7);
	ASSERT_EQ(Halfedge.Origin(14), 7);
}

TEST(HalfedgeMeshWrapper, Function_DestVertex)
{
	GEO::Mesh Mesh;
	ASSERT_TRUE(MeshGenerator::MeshGenHexagon(&Mesh));
	HalfedgeMeshWrapper Halfedge(&Mesh);

	ASSERT_EQ(Halfedge.Dest(10), 6);
	ASSERT_EQ(Halfedge.Dest(1), 7);
	ASSERT_EQ(Halfedge.Dest(11), 1);

	ASSERT_EQ(Halfedge.Dest(12), 2);

	ASSERT_EQ(Halfedge.Dest(13), 7);
	ASSERT_EQ(Halfedge.Dest(16), 3);

	ASSERT_EQ(Halfedge.Dest(17), 7);
	ASSERT_EQ(Halfedge.Dest(3), 6);
	ASSERT_EQ(Halfedge.Dest(4), 4);

	ASSERT_EQ(Halfedge.Dest(5), 5);

	ASSERT_EQ(Halfedge.Dest(6), 6);
	ASSERT_EQ(Halfedge.Dest(9), 0);

	ASSERT_EQ(Halfedge.Dest(0), 0);
	ASSERT_EQ(Halfedge.Dest(7), 3);
	ASSERT_EQ(Halfedge.Dest(8), 5);

	ASSERT_EQ(Halfedge.Dest(2), 3);
	ASSERT_EQ(Halfedge.Dest(15), 2);
	ASSERT_EQ(Halfedge.Dest(14), 0);
}