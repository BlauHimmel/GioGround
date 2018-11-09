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

	bool MeshUVSphere(InOut GEO::Mesh * pMesh, GEO::index_t nSegments, GEO::index_t nRings, double Radius, double PoleEpsilon)
	{
		if (pMesh == nullptr)
		{
			return false;
		}

		pMesh->clear(false, false);

		GEO::index_t nVertex = (nSegments + 1) * (nRings + 1);
		GEO::index_t nFacet = nSegments * nRings * 2;
		
		pMesh->vertices.create_vertices(nVertex);

		GEO::Attribute<double> AttriVertexUV;
		AttriVertexUV.bind(pMesh->vertices.attributes(), "tex_coord");
		AttriVertexUV.redim(2);

		GEO::index_t iVertex = 0;

		for (GEO::index_t i = 0; i <= nRings; i++)
		{
			double Theta = (double(i) / double(nRings)) * (M_PI - 2.0 * PoleEpsilon) + PoleEpsilon;
			for (GEO::index_t j = 0; j <= nSegments; j++)
			{
				double Phi = (double(j) / double(nSegments)) * 2.0 * M_PI;
		
				double X = Radius * std::sin(Theta) * std::cos(Phi);
				double Y = Radius * std::sin(Theta) * std::sin(Phi);
				double Z = Radius * std::cos(Theta);

				double U = Phi / (2.0 * M_PI);
				double V = 1.0 - Theta / M_PI;

				double * pV = pMesh->vertices.point_ptr(iVertex);
				pV[0] = X;
				pV[1] = Y;
				pV[2] = Z;

				AttriVertexUV[2 * iVertex + 0] = U;
				AttriVertexUV[2 * iVertex + 1] = V;

				iVertex++;
			}
		}

		for (GEO::index_t i = 0; i < nRings; i++)
		{
			for (GEO::index_t j = 0; j < nSegments; j++)
			{
				GEO::index_t iLB = i * (nSegments + 1) + j;
				GEO::index_t iRB = iLB + 1;
				GEO::index_t iLT = iLB + (nSegments + 1);
				GEO::index_t iRT = iLT + 1;

				pMesh->facets.create_triangle(iRT, iRB, iLB);
				pMesh->facets.create_triangle(iLB, iLT, iRT);
			}
		}

		pMesh->facets.connect();
		pMesh->vertices.set_double_precision();

		return true;
	}
}