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

		GEO::Logger::out("LSCM") << "Buiding and solving the linear system...";
		TIMER_START(TimeSolveLeastSquareEquation);
		SolveLeastSquareEquation(pMesh);
		TIMER_END(TimeSolveLeastSquareEquation);
		GEO::Logger::out("LSCM") << "(" << TimeSolveLeastSquareEquation << " sec)" << std::endl;

		GEO::Logger::out("LSCM") << "Normalizing UV...";
		TIMER_START(TimeNormalizeUV);
		NormalizeUV(pMesh);
		TIMER_END(TimeNormalizeUV);
		GEO::Logger::out("LSCM") << "(" << TimeNormalizeUV << " sec)" << std::endl;

		GEO::Logger::out("LSCM") << "Generating output...";
		TIMER_START(TimeGenerateOutput);
		GenerateOutput(pMesh);
		TIMER_END(TimeGenerateOutput);
		GEO::Logger::out("LSCM") << "(" << TimeGenerateOutput << " sec)" << std::endl;

		GEO::Logger::out("LSCM") << "Clearing attributes...";
		TIMER_START(TimeClearAttribute);
		ClearAttribute(pMesh);
		TIMER_END(TimeClearAttribute);
		GEO::Logger::out("LSCM") << "(" << TimeClearAttribute << " sec)" << std::endl;

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

		double X0 = 0.0;
		double Y0 = 0.0;
		double X1 = (P1 - P0).length();
		double Y1 = 0.0;
		double X2 = GEO::dot((P2 - P0), X);
		double Y2 = GEO::dot((P2 - P0), Y);

		Z0 = GEO::vec2(X0, Y0);
		Z1 = GEO::vec2(X1, Y1);
		Z2 = GEO::vec2(X2, Y2);
	}

	void LSCMAlgorithm::GetBoundingBox(
		In GEO::Mesh * const pMesh,
		Out double & XMin,
		Out double & YMin,
		Out double & ZMin,
		Out double & XMax,
		Out double & YMax,
		Out double & ZMax
	) const
	{
		assert(pMesh != nullptr);

		XMin = 1e30;
		YMin = 1e30;
		ZMin = 1e30;
		XMax = -1e30;
		YMax = -1e30;
		ZMax = -1e30;

		for (GEO::index_t i = 0; i < pMesh->vertices.nb(); i++)
		{
			GEO::vec3 Pt = pMesh->vertices.point(i);

			XMin = std::min(XMin, Pt.x);
			YMin = std::min(YMin, Pt.y);
			ZMin = std::min(ZMin, Pt.z);
			XMax = std::max(XMax, Pt.x);
			YMax = std::max(YMax, Pt.y);
			ZMax = std::max(ZMax, Pt.z);
		}
	}

	void LSCMAlgorithm::GetLargestDimension(
		In double XMin,
		In double YMin,
		In double ZMin,
		In double XMax,
		In double YMax,
		In double ZMax,
		Out GEO::vec3 & V1,
		Out GEO::vec3 & V2
	) const
	{
		double Dx = XMax - XMin;
		double Dy = YMax - YMin;
		double Dz = ZMax - ZMin;

		if (Dx <= Dy && Dx <= Dz)
		{
			if (Dy > Dz)
			{
				V1 = GEO::vec3(0.0, 1.0, 0.0);
				V2 = GEO::vec3(0.0, 0.0, 1.0);
			}
			else
			{
				V2 = GEO::vec3(0.0, 1.0, 0.0);
				V1 = GEO::vec3(0.0, 0.0, 1.0);
			}
		}
		else if (Dy <= Dx && Dy <= Dz)
		{
			if (Dx > Dz)
			{
				V1 = GEO::vec3(1.0, 0.0, 0.0);
				V2 = GEO::vec3(0.0, 0.0, 1.0);
			}
			else
			{
				V2 = GEO::vec3(1.0, 0.0, 0.0);
				V1 = GEO::vec3(0.0, 0.0, 1.0);
			}
		}
		else if (Dz <= Dx && Dz <= Dy)
		{
			if (Dx > Dy)
			{
				V1 = GEO::vec3(1.0, 0.0, 0.0);
				V2 = GEO::vec3(0.0, 1.0, 0.0);
			}
			else
			{
				V2 = GEO::vec3(1.0, 0.0, 0.0);
				V1 = GEO::vec3(0.0, 1.0, 0.0);
			}
		}
	}

	void LSCMAlgorithm::FixedPoint(
		In GEO::Mesh * const pMesh,
		In GEO::vec3 V1,
		In GEO::vec3 V2,
		Out GEO::index_t & iFixedPointMin,
		Out GEO::index_t & iFixedPointMax
	) const
	{
		assert(pMesh != nullptr);

		double UMin = 1e30;
		double UMax = -1e30;

		GEO::Attribute<double> AttriVertexUV;
		AttriVertexUV.bind(pMesh->vertices.attributes(), "UV");
		AttriVertexUV.redim(2);

		for (GEO::index_t i = 0; i < pMesh->vertices.nb(); i++)
		{
			GEO::vec3 Pt = pMesh->vertices.point(i);
			double U = GEO::dot(Pt, V1);
			double V = GEO::dot(Pt, V2);
			
			AttriVertexUV[2 * i + 0] = U;
			AttriVertexUV[2 * i + 1] = V;

			if (U < UMin)
			{
				iFixedPointMin = i;
				UMin = U;
			}

			if (U > UMax)
			{
				iFixedPointMax = i;
				UMax = U;
			}
		}
	}

	void LSCMAlgorithm::SolveLeastSquareEquation(In GEO::Mesh * const pMesh)
	{
		assert(pMesh != nullptr);

		GEO::index_t nFacet = pMesh->facets.nb();
		GEO::index_t nVertex = pMesh->vertices.nb();

		double XMin, YMin, ZMin;
		double XMax, YMax, ZMax;
		GEO::vec3 V1, V2;
		GEO::index_t iFixedPointMax, iFixedPointMin;

		GetBoundingBox(pMesh, XMin, YMin, ZMin, XMax, YMax, ZMax);
		GetLargestDimension(XMin, YMin, ZMin, XMax, YMax, ZMax, V1, V2);
		FixedPoint(pMesh, V1, V2, iFixedPointMin, iFixedPointMax);

		nlNewContext();

		nlSolverParameteri(NL_NB_VARIABLES, NLint(2 * nVertex));
		nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
		nlSolverParameteri(NL_MAX_ITERATIONS, NLint(5 * nVertex));
		nlSolverParameterd(NL_THRESHOLD, 1e-6);

		nlBegin(NL_SYSTEM);
		GEO::Attribute<double> AttriVertexUV(pMesh->vertices.attributes(), "UV");
		for (GEO::index_t i = 0; i < pMesh->vertices.nb(); i++)
		{
			double U = AttriVertexUV[2 * i + 0];
			double V = AttriVertexUV[2 * i + 1];

			nlSetVariable(2 * i + 0, U);
			nlSetVariable(2 * i + 1, V);
		}
		nlLockVariable(2 * iFixedPointMin + 0);
		nlLockVariable(2 * iFixedPointMin + 1);
		nlLockVariable(2 * iFixedPointMax + 0);
		nlLockVariable(2 * iFixedPointMax + 1);

		nlBegin(NL_MATRIX);
		for (GEO::index_t iFacet = 0; iFacet < pMesh->facets.nb(); iFacet++)
		{
			GEO::index_t nVertex = pMesh->facets.nb_vertices(iFacet);
			assert(nVertex == 3);

			GEO::index_t iP0 = pMesh->facets.vertex(iFacet, 0);
			GEO::index_t iP1 = pMesh->facets.vertex(iFacet, 1);
			GEO::index_t iP2 = pMesh->facets.vertex(iFacet, 2);
			
			const GEO::vec3 & P0 = pMesh->vertices.point(iP0);
			const GEO::vec3 & P1 = pMesh->vertices.point(iP1);
			const GEO::vec3 & P2 = pMesh->vertices.point(iP2);

			GEO::vec2 Z0, Z1, Z2;
			ProjectTriangle(P0, P1, P2, Z0, Z1, Z2);
			
			GEO::vec2 Z01 = Z1 - Z0;
			GEO::vec2 Z02 = Z2 - Z0;

			double A = Z01.x;
			double B = Z01.y;
			double C = Z02.x;
			double D = Z02.y;
			assert(B == 0.0);

			GEO::index_t iU0 = 2 * iP0 + 0;
			GEO::index_t iV0 = 2 * iP0 + 1;
			GEO::index_t iU1 = 2 * iP1 + 0;
			GEO::index_t iV1 = 2 * iP1 + 1;
			GEO::index_t iU2 = 2 * iP2 + 0;
			GEO::index_t iV2 = 2 * iP2 + 1;

			/* Real part */
			nlBegin(NL_ROW);
			nlCoefficient(iU0, -A + C);
			nlCoefficient(iV0,  B - D);
			nlCoefficient(iU1, -C    );
			nlCoefficient(iV1,  D    );
			nlCoefficient(iU2,  A    );
			nlCoefficient(iV2,  B    );
			nlEnd(NL_ROW);

			/* Imaginary part */
			nlBegin(NL_ROW);
			nlCoefficient(iU0, -B + D);
			nlCoefficient(iV0, -A + C);
			nlCoefficient(iU1, -D    );
			nlCoefficient(iV1, -C    );
			nlCoefficient(iU2, -B    );
			nlCoefficient(iV2,  A    );
			nlEnd(NL_ROW);
		}
		nlEnd(NL_MATRIX);
		nlEnd(NL_SYSTEM);

		nlSolve();

		m_Vertices.resize(pMesh->vertices.nb() * 2);
		for (GEO::index_t i = 0; i < pMesh->vertices.nb(); i++)
		{
			double U = nlGetVariable(2 * i + 0);
			double V = nlGetVariable(2 * i + 1);
			m_Vertices[2 * i + 0] = U;
			m_Vertices[2 * i + 1] = V;
		}

		nlDeleteContext(nlGetCurrent());
	}

	void LSCMAlgorithm::NormalizeUV(In GEO::Mesh * const pMesh)
	{
		assert(pMesh != nullptr);

		double UMin = 1e30;
		double VMin = 1e30;
		double UMax = -1e30;
		double VMax = -1e30;

		for (GEO::index_t i = 0; i < pMesh->vertices.nb(); i++)
		{
			UMin = std::min(UMin, m_Vertices[2 * i + 0]);
			VMin = std::min(VMin, m_Vertices[2 * i + 1]);
			UMax = std::max(UMax, m_Vertices[2 * i + 0]);
			VMax = std::max(VMax, m_Vertices[2 * i + 1]);
		}

		double Inv = 1.0 / std::max(UMax - UMin, VMax - VMin);

		for (GEO::index_t i = 0; i < pMesh->vertices.nb(); i++)
		{
			m_Vertices[2 * i + 0] -= UMin;
			m_Vertices[2 * i + 0] *= Inv;
			m_Vertices[2 * i + 1] -= VMin;
			m_Vertices[2 * i + 1] *= Inv;
		}
	}

	void LSCMAlgorithm::GenerateOutput(In GEO::Mesh * const pMesh)
	{
		assert(pMesh != nullptr);

		pMesh->vertices.set_dimension(6);

		GEO::AttributesManager & CornerAttributeManager = pMesh->facet_corners.attributes();
		if (CornerAttributeManager.is_defined("tex_coord"))
		{
			CornerAttributeManager.delete_attribute_store("tex_coord");
		}

		GEO::Attribute<double> AttriVertexTexCoord(pMesh->vertices.attributes(), "tex_coord");
		AttriVertexTexCoord.redim(2);

		double XYZMin[3];
		double XYZMax[3];
		GetBBox(pMesh, XYZMin, XYZMax, false);
		double DeltaX = XYZMax[0] - XYZMin[0];
		double DeltaY = XYZMax[1] - XYZMin[1];
		double Delta = std::max(DeltaX, DeltaY);
		double HalfDelta = 0.5  * Delta;

		for (GEO::index_t iVertex = 0; iVertex < pMesh->vertices.nb(); iVertex++)
		{
			double * pVertex = pMesh->vertices.point_ptr(iVertex);

			pVertex[3] = m_Vertices[iVertex * 2 + 0] * Delta - HalfDelta;
			pVertex[4] = m_Vertices[iVertex * 2 + 1] * Delta - HalfDelta;

			AttriVertexTexCoord[iVertex * 2 + 0] = m_Vertices[iVertex * 2 + 0];
			AttriVertexTexCoord[iVertex * 2 + 1] = m_Vertices[iVertex * 2 + 1];
		}
	}

	void LSCMAlgorithm::ClearAttribute(In GEO::Mesh * const pMesh)
	{
		assert(pMesh != nullptr);

		GEO::AttributesManager & VertexAttributeManager = pMesh->vertices.attributes();

		GEO::vector<std::string> AttributeNames;

		VertexAttributeManager.list_attribute_names(AttributeNames);
		for (std::string & Name : AttributeNames)
		{
			if (Name == "UV")
			{
				VertexAttributeManager.delete_attribute_store(Name);
			}
		}
	}
}
