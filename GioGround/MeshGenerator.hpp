#pragma once

#include <geogram\mesh\mesh.h>

#include "Macro.hpp"

namespace MeshGenerator
{
	bool MeshGenHexagon(InOut GEO::Mesh * pMesh);

	bool MeshGenPyramid(InOut GEO::Mesh * pMesh);
}