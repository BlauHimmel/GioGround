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

		pMesh->vertices.create_vertices(nPoints);
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

	bool MeshGenPyramid(InOut GEO::Mesh * pMesh)
	{
		if (pMesh == nullptr)
		{
			return false;
		}

		static double pPoints[] =
		{
			 3.0,  0.0, 1.0,
			-1.0, -1.0, 0.0,
			 1.0, -1.0, 0.0,
			 1.0,  1.0, 0.0,
			-1.0,  1.0, 0.0
		};

		static GEO::index_t nDim = 3;
		static GEO::index_t nPoints = 5;

		pMesh->clear(false, false);

		pMesh->vertices.create_vertices(nPoints);
		pMesh->vertices.assign_points(pPoints, nDim, nPoints);

		pMesh->facets.create_triangle(0, 1, 2);
		pMesh->facets.create_triangle(0, 2, 3);
		pMesh->facets.create_triangle(0, 3, 4);
		pMesh->facets.create_triangle(0, 4, 1);

		pMesh->facets.connect();
		pMesh->vertices.set_double_precision();

		return true;
	}

	bool MeshUVSphere(InOut GEO::Mesh * pMesh, GEO::index_t nSegments, GEO::index_t nRings)
	{
		if (pMesh == nullptr)
		{
			return false;
		}

		pMesh->clear(false, false);

		GEO::index_t nVertex = (nSegments + 1) * (nRings + 1);
		GEO::index_t nFacet = nSegments * nRings * 2;
		
		pMesh->vertices.create_vertices(nVertex);
		pMesh->facets.create_facets(nFacet, 3);

		GEO::Attribute<double> AttriVertexUV;
		AttriVertexUV.bind(pMesh->vertices.attributes(), "tex_coord");
		AttriVertexUV.redim(2);

		double DeltaTheta = M_PI / nSegments;
		double DeltaPhi = 2.0 * M_PI / nRings;

		GEO::index_t iVertex = 0;

		for (double Theta = 0.0; Theta <= M_PI; Theta += DeltaTheta)
		{
			for (double Phi = 0.0; Phi <= 2.0 * M_PI; Phi += DeltaPhi)
			{
				double X = std::sin(Theta) * std::cos(Phi);
				double Y = std::sin(Theta) * std::sin(Phi);
				double Z = std::cos(Theta);

				double U = Phi;
				double V = Theta;

				double * pV = pMesh->vertices.point_ptr(iVertex);
				pV[0] = X;
				pV[1] = Y;
				pV[2] = Z;

				AttriVertexUV[2 * iVertex + 0] = U;
				AttriVertexUV[2 * iVertex + 1] = V;

				iVertex++;
			}
		}

		return true;
	}
}