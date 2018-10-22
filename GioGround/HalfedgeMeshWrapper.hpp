#pragma once

#include <geogram\mesh\mesh.h>

#include <memory>

#include "Macro.hpp"

// Note : the sequence of iteration is CCW
struct HalfedgeMeshWrapper
{
	GEO::Mesh * pMesh = nullptr;
	GEO::vector<GEO::index_t> Vertex2Corner;
	GEO::vector<GEO::index_t> Corner2Facet;
	GEO::vector<GEO::index_t> Corner2Corner;

	HalfedgeMeshWrapper(Ref GEO::Mesh * pMesh);

	void Reset();
	GEO::index_t Next(In GEO::index_t iCorner);
	GEO::index_t Prev(In GEO::index_t iCorner);
	GEO::index_t Opposite(In GEO::index_t iCorner);
	GEO::index_t Facet(In GEO::index_t iCorner);
	GEO::index_t NextAroundVertex(In GEO::index_t iCorner);
	GEO::index_t Origin(In GEO::index_t iCorner);
	GEO::index_t Dest(In GEO::index_t iCorner);
};