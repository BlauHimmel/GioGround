#pragma once

#include <geogram\mesh\mesh.h>

#include "Macro.hpp"

namespace MeshGenerator
{
	bool MeshGenHexagon(InOut GEO::Mesh * pMesh);

	bool MeshGenPyramid(InOut GEO::Mesh * pMesh);

	bool MeshUVSphere(InOut GEO::Mesh * pMesh, GEO::index_t nSegments, GEO::index_t nRings);
}