#include "LSCMAlgorithm.hpp"

#include "..\HalfedgeMeshWrapper.hpp"

#include <geogram\NL\nl.h>
#include <mkl.h>

namespace MeshAlgorithm
{
	bool LSCMAlgorithm::Execute(InOut GEO::Mesh * pMesh)
	{
		if (pMesh == nullptr)
		{
			return false;
		}

		if (!Reset())
		{
			return false;
		}

		if (!CheckAndGetArgs(pMesh))
		{
			return false;
		}

		TIMER_START(TimeTotal);

		/* TODO */

		TIMER_END(TimeTotal);

		GEO::Logger::out("LSCM") << "LSCM finished! (Total : " << TimeTotal << " sec)" << std::endl;
		return true;
	}

	bool LSCMAlgorithm::Visualize(In GEO::Mesh * const pMesh) const
	{
		return true;
	}

	bool LSCMAlgorithm::CheckAndGetArgs(In GEO::Mesh * const pMesh)
	{
		pMesh->assert_is_valid();
		assert(pMesh->vertices.double_precision());

		if (GetBoundaryNumber(pMesh) != 1)
		{
			return false;
		}

		return true;
	}

	bool LSCMAlgorithm::Reset()
	{
		return true;
	}

	void LSCMAlgorithm::ProjectTriangle(
		In const GEO::vec3 & P0,
		In const GEO::vec3 & P1,
		In const GEO::vec3 & P2,
		Out GEO::vec2 & Z0,
		Out GEO::vec2 & Z1,
		Out GEO::vec2 & Z2
	) const
	{
		GEO::vec3 X = GEO::normalize(P1 - P0);
		GEO::vec3 Z = GEO::normalize(GEO::cross(X, P2 - P0));
		GEO::vec3 Y = GEO::cross(Z, X);

		double X0 = 0;
		double Y0 = 0;
		double X1 = (P1 - P0).length();
		double Y1 = 0;
		double X2 = GEO::dot((P2 - P0), X);
		double Y2 = GEO::dot((P2 - P0), Y);

		Z0 = GEO::vec2(X0, Y0);
		Z1 = GEO::vec2(X1, Y1);
		Z2 = GEO::vec2(X2, Y2);
	}

	void LSCMAlgorithm::SolveLeastSquareEquation(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper)
	{
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);

		GEO::Mesh * pMesh = pHalfedgeMeshWrapper->pMesh;

		GEO::index_t nFacet = pMesh->facets.nb();
		GEO::index_t nVertex = pMesh->vertices.nb();
		const GEO::index_t nFixedPoint = 2;

		int m = 2 * int(nFacet);
		int n = 2 * int(nVertex - nFixedPoint);

		double * A/*Size = m * n*/ = reinterpret_cast<double*>(mkl_malloc(sizeof(double) * m * n, 64));
		double * B/*Size = m * 1*/ = reinterpret_cast<double*>(mkl_malloc(sizeof(double) * m, 64));

		LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', m, n, 1, A, n, B, 1);

		mkl_free(A);
		mkl_free(B);
	}

}
