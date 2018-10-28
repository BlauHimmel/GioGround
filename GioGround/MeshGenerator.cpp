#include "MeshGenerator.hpp"

#include <geogram\mesh\mesh_reorder.h>

namespace MeshGenerator
{
	bool MeshGenHexagon(InOut GEO::Mesh * pMesh)
	{
		if (pMesh == nullptr)
		{
			return false;
		}

		static double pPoints[] =
		{
			-0.5,   0.0,  0.0,
			-0.25, -0.5,  0.0,
			 0.25, -0.5,  0.0,
			 0.5,   0.0,  0.0,
			 0.25,  0.5,  0.0,
			-0.25,  0.5,  0.0,
			 0.0,   0.25, 0.0,
			 0.0,  -0.25, 0.0
		};

		static GEO::index_t nDim = 3;
		static GEO::index_t nPoints = 8;

		pMesh->clear(false, false);

		pMesh->vertices.create_vertices(8);
		pMesh->vertices.assign_points(pPoints, nDim, nPoints);

		pMesh->facets.create_quad(6, 0, 7, 3);
		pMesh->facets.create_quad(3, 4, 5, 6);
		pMesh->facets.create_triangle(6, 5, 0);
		pMesh->facets.create_quad(0, 1, 2, 7);
		pMesh->facets.create_triangle(7, 2, 3);

		pMesh->facets.connect();
		pMesh->vertices.set_double_precision();

		return true;
	}
}