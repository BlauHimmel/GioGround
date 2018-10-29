#include "BarycentricMappingAlgorithm.hpp"

#include "..\HalfedgeMeshWrapper.hpp"

#include <mkl.h>

#include <unordered_set>
#include <functional>

namespace MeshAlgorithm
{
	std::string BarycentricMappingAlgorithm::PARAMS_KEY_COEFFICIENT_TYPE = "CoefficientType";
	std::string BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[3] = 
	{
		"Mean Value Coordinates",
		"Discrete Harmonic Coordinates",
		"Wachspress Coordinates"
	};

	std::string BarycentricMappingAlgorithm::PARAMS_KEY_DOMAIN_SHAPE = "DomainShape";
	std::string BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_DOMAIN_SHAPE[2] = 
	{
		"Square",
		"Circle"
	};

	std::string BarycentricMappingAlgorithm::PARAMS_KEY_BOUNDARY_FIX_WEIGHT = "BoundaryFixWeight";
	std::string BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_BOUNDARY_FIX_WEIGHT[2] = 
	{
		"Mean",
		"Length"
	};

	bool BarycentricMappingAlgorithm::Execute(InOut GEO::Mesh * pMesh)
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

		GEO::Logger::out("Barycentric") << "Create Halfedge data structure...";
		TIMER_START(TimeCreateHalfedge);
		HalfedgeMeshWrapper Wrapper(pMesh);
		TIMER_END(TimeCreateHalfedge);
		GEO::Logger::out("Barycentric") << "(" << TimeCreateHalfedge << " sec)" << std::endl;

		GEO::Logger::out("Barycentric") << "Executing vertices partition...";
		TIMER_START(TimeFindInteriorBoundary);
		FindInteriorBoundary(&Wrapper);
		TIMER_END(TimeFindInteriorBoundary);
		GEO::Logger::out("Barycentric") << "(" << TimeFindInteriorBoundary << " sec)" << std::endl;

		GEO::Logger::out("Barycentric") << "Fixing boundary vertices...";
		TIMER_START(TimeFixBoundaryVertices);
		FixBoundaryVertices(&Wrapper);
		TIMER_END(TimeFixBoundaryVertices);
		GEO::Logger::out("Barycentric") << "(" << TimeFixBoundaryVertices << " sec)" << std::endl;

		GEO::Logger::out("Barycentric") << "Building and solving linear equation...";
		TIMER_START(TimeSoloveLinearEquation);
		SoloveLinearEquation(&Wrapper);
		TIMER_END(TimeSoloveLinearEquation);
		GEO::Logger::out("Barycentric") << "(" << TimeSoloveLinearEquation << " sec)" << std::endl;

		GEO::Logger::out("Barycentric") << "Generating output...";
		TIMER_START(TimeGenerateOutput);
		GenerateOutput(&Wrapper);
		TIMER_END(TimeGenerateOutput);
		GEO::Logger::out("Barycentric") << "(" << TimeGenerateOutput << " sec)" << std::endl;

		TIMER_END(TimeTotal);

		GEO::Logger::out("Barycentric") << "Barycentric mapping finished! (Total : " << TimeTotal << " sec)" << std::endl;
		return true;
	}

	bool BarycentricMappingAlgorithm::Visualize(In GEO::Mesh * const pMesh) const
	{
		return true;
	}

	bool BarycentricMappingAlgorithm::CheckAndGetArgs(In GEO::Mesh * const pMesh)
	{
		pMesh->assert_is_valid();
		assert(pMesh->vertices.double_precision());

		if (!GetArg(PARAMS_KEY_COEFFICIENT_TYPE, &m_CoefficientType))
		{
			return false;
		}

		if (m_CoefficientType != PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[0] &&
			m_CoefficientType != PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[1] &&
			m_CoefficientType != PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[2])
		{
			return false;
		}

		if (!GetArg(PARAMS_KEY_DOMAIN_SHAPE, &m_DomainShape))
		{
			return false;
		}

		if (m_DomainShape != PARAMS_VALUE_SUPPORTED_DOMAIN_SHAPE[0] &&
			m_DomainShape != PARAMS_VALUE_SUPPORTED_DOMAIN_SHAPE[1])
		{
			return false;
		}

		if (!GetArg(PARAMS_KEY_BOUNDARY_FIX_WEIGHT, &m_BoundaryFixWeight))
		{
			return false;
		}

		if (m_BoundaryFixWeight != PARAMS_VALUE_SUPPORTED_BOUNDARY_FIX_WEIGHT[0] &&
			m_BoundaryFixWeight != PARAMS_VALUE_SUPPORTED_BOUNDARY_FIX_WEIGHT[1])
		{
			return false;
		}

		if (GetBoundaryNumber(pMesh) != 1)
		{
			return false;
		}

		return true;
	}

	bool BarycentricMappingAlgorithm::Reset()
	{
		m_CoefficientType = PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[0];
		m_DomainShape = PARAMS_VALUE_SUPPORTED_DOMAIN_SHAPE[0];
		m_BoundaryFixWeight = PARAMS_VALUE_SUPPORTED_BOUNDARY_FIX_WEIGHT[0];
		m_nInteriorVertices = GEO::index_t(-1);
		m_nBoundaryVertices = GEO::index_t(-1);
		m_iInteriorVertices.clear();
		m_iBoundaryVertices.clear();
		m_InteriorVertices.clear();
		m_BoundaryVertices.clear();
		return true;
	}

	void BarycentricMappingAlgorithm::FindInteriorBoundary(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper)
	{
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);

		GEO::Mesh * pMesh = pHalfedgeMeshWrapper->pMesh;
		
		GEO::vector<GEO::index_t> BoundaryCornersNoOrder;

		for (GEO::index_t iVertex = 0; iVertex < pMesh->vertices.nb(); iVertex++)
		{
			bool bIsBoundaryVertex = false;
			GEO::index_t iBeginCorner = pHalfedgeMeshWrapper->Vertex2Corner[iVertex];
			GEO::index_t iCorner = iBeginCorner;
			do
			{
				if (pMesh->facet_corners.adjacent_facet(iCorner) == GEO::NO_CORNER)
				{
					bIsBoundaryVertex = true;
					break;
				}
				iCorner = pHalfedgeMeshWrapper->Corner2Corner[iCorner];
			} while (iCorner != iBeginCorner);

			if (bIsBoundaryVertex)
			{
				BoundaryCornersNoOrder.push_back(iCorner);
			}
			else
			{
				m_iInteriorVertices.push_back(iVertex);
			}
		}

		std::unordered_set<GEO::index_t> BoundaryCornersSet(BoundaryCornersNoOrder.begin(), BoundaryCornersNoOrder.end());

		GEO::index_t iCornerBegin = BoundaryCornersNoOrder[0];
		GEO::index_t iCorner = iCornerBegin;

		do
		{
			m_iBoundaryVertices.push_back(pMesh->facet_corners.vertex(iCorner));
			
			iCorner = pHalfedgeMeshWrapper->Next(iCorner);
			while (BoundaryCornersSet.find(iCorner) == BoundaryCornersSet.end())
			{
				iCorner = pHalfedgeMeshWrapper->Corner2Corner[iCorner];
			}
		} while (iCorner != iCornerBegin);

		m_nInteriorVertices = m_iInteriorVertices.size();
		m_nBoundaryVertices = m_iBoundaryVertices.size();
	}

	void BarycentricMappingAlgorithm::FixBoundaryVertices(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper)
	{
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);

		GEO::Mesh * pMesh = pHalfedgeMeshWrapper->pMesh;

		GEO::vector<double> Weight;
		if (m_BoundaryFixWeight == PARAMS_VALUE_SUPPORTED_BOUNDARY_FIX_WEIGHT[0])
		{
			double MeanWeight = 1.0 / double(m_nBoundaryVertices);
			Weight.resize(m_nBoundaryVertices, MeanWeight);
		}
		else if (m_BoundaryFixWeight == PARAMS_VALUE_SUPPORTED_BOUNDARY_FIX_WEIGHT[1])
		{
			Weight.resize(m_nBoundaryVertices);
			double Total = 0.0;
			for (GEO::index_t i = 0; i < m_nBoundaryVertices - 1; i++)
			{
				GEO::vec3 Vertex1 = pMesh->vertices.point(m_iInteriorVertices[i]);
				GEO::vec3 Vertex2 = pMesh->vertices.point(m_iInteriorVertices[i + 1]);
				double Distance = std::sqrt(std::pow(Vertex1.x - Vertex2.x, 2.0) + std::pow(Vertex1.y - Vertex2.y, 2.0) + std::pow(Vertex1.z - Vertex2.z, 2.0));
				Total += Distance;
				Weight[i] = Distance;
			}

			GEO::vec3 Vertex1 = pMesh->vertices.point(m_iInteriorVertices[m_nBoundaryVertices - 1]);
			GEO::vec3 Vertex2 = pMesh->vertices.point(m_iInteriorVertices[0]);
			double Distance = std::sqrt(std::pow(Vertex1.x - Vertex2.x, 2.0) + std::pow(Vertex1.y - Vertex2.y, 2.0) + std::pow(Vertex1.z - Vertex2.z, 2.0));
			Total += Distance;
			Weight[m_nBoundaryVertices - 1] = Distance;

			double InvTotal = 1.0 / Total;
			for (GEO::index_t i = 0; i < m_nBoundaryVertices; i++)
			{
				Weight[i] *= InvTotal;
			}
		}

		if (m_DomainShape == PARAMS_VALUE_SUPPORTED_DOMAIN_SHAPE[0])
		{
			FixSquareBoundaryVertices(Weight);
		}
		else if (m_DomainShape == PARAMS_VALUE_SUPPORTED_DOMAIN_SHAPE[1])
		{
			FixSphereBoundaryVertices(Weight);
		}
		else
		{
			/*m_DomainShape must be in PARAMS_VALUE_SUPPORTED_DOMAIN_SHAPE*/
			assert(false);
		}
	}

	void BarycentricMappingAlgorithm::FixSquareBoundaryVertices(In GEO::vector<double> & Weight)
	{
		double MaxX = 1.0, MinX = 0.0;
		double MaxY = 1.0, MinY = 0.0;
		double X = MinX, Y = MinY, Z = 0.0;
		double BoundaryLength = (MaxX - MinX + MaxY - MinY) * 2.0;
		GEO::index_t iDirection = 0;
		double Compensate = 0.0;

		m_BoundaryVertices.resize(m_nBoundaryVertices * 3, 0.0);

		for (GEO::index_t i = 0; i < Weight.size(); i++)
		{
			double Step = BoundaryLength * Weight[i];

			m_BoundaryVertices[3 * i + 0] = X;
			m_BoundaryVertices[3 * i + 1] = Y;
			m_BoundaryVertices[3 * i + 2] = Z;

			if (iDirection == 0)
			{
				if (X + Step > MaxX)
				{
					Compensate = Step - (MaxX - X);
					X = MaxX;
					iDirection++;
				}
				else
				{
					X += Step;
					if (X == MaxX)
					{
						iDirection++;
					}
				}
			}
			else if (iDirection == 1)
			{
				if (Y + Step > MaxY)
				{
					Compensate = Step - (MaxY - Y);
					Y = MaxY;
					iDirection++;
				}
				else
				{
					Y += Step;
					if (Compensate != 0.0)
					{
						Y += Compensate;
						Compensate = 0.0;
					}
					if (Y == MaxY)
					{
						iDirection++;
					}
				}
			}
			else if (iDirection == 2)
			{
				if (X - Step < MinX)
				{
					Compensate = Step - (X - MinX);
					X = MinX;
					iDirection++;
				}
				else
				{
					X -= Step;
					if (Compensate != 0.0)
					{
						X -= Compensate;
						Compensate = 0.0;
					}
					if (X == MinX)
					{
						iDirection++;
					}
				}
			}
			else if (iDirection == 3)
			{
				if (Y - Step < MinY)
				{
					Compensate = Step - (Y - MinY);
					Y = MinY;
					iDirection++;
				}
				else
				{
					Y -= Step;
					if (Compensate != 0.0)
					{
						Y -= Compensate;
						Compensate = 0.0;
					}
					if (Y == MinY)
					{
						iDirection++;
					}
				}
			}
		}
	}

	void BarycentricMappingAlgorithm::FixSphereBoundaryVertices(In GEO::vector<double> & Weight)
	{
		double Radius = 0.5;
		double Theta = 0.0;
		double Parameters = 2.0 * M_PI;

		m_BoundaryVertices.resize(m_nBoundaryVertices * 3, 0.0);

		for (GEO::index_t i = 0; i < Weight.size(); i++)
		{
			double X = Radius * std::cos(Theta);
			double Y = Radius * std::sin(Theta);
			double Z = 0.0;

			m_BoundaryVertices[3 * i + 0] = X + Radius;
			m_BoundaryVertices[3 * i + 1] = Y + Radius;
			m_BoundaryVertices[3 * i + 2] = Z;

			double Step = Parameters * Weight[i];
			Theta += Step;
		}
	}

	void BarycentricMappingAlgorithm::SoloveLinearEquation(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper)
	{
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);

		GEO::Mesh * pMesh = pHalfedgeMeshWrapper->pMesh;

		int n = int(m_nInteriorVertices);
		int b = int(m_nBoundaryVertices);

		double * A/*Size = n * n*/ = reinterpret_cast<double*>(mkl_malloc(sizeof(double) * n * n, 64));
		double * B/*Size = n * 2*/ = reinterpret_cast<double*>(mkl_malloc(sizeof(double) * n * 2, 64));

		memset(A, 0, sizeof(double) * n * n);
		memset(B, 0, sizeof(double) * n * 2);

		for (GEO::index_t i = 0; i < GEO::index_t(n); i++)
		{
			A[i * n + i] = 1.0;
			GEO::vector<GEO::index_t> iAdjVertices = GetAdjacentVertices(pHalfedgeMeshWrapper, m_iInteriorVertices[i]);
			for (GEO::index_t j = 0; j < iAdjVertices.size(); j++)
			{
				auto Iter = std::find(m_iInteriorVertices.begin(), m_iInteriorVertices.end(), iAdjVertices[j]);
				if (Iter != m_iInteriorVertices.end())
				{
					double Lambda_ij = Lambda_ij_BarycentricCoordinates(pHalfedgeMeshWrapper, m_iInteriorVertices[i], iAdjVertices[j]);
					A[i * n + (Iter - m_iInteriorVertices.begin())] = -1.0 * Lambda_ij;
				}
			}
		}

		for (GEO::index_t i = 0; i < GEO::index_t(n); i++)
		{
			GEO::vector<GEO::index_t> iAdjVertices = GetAdjacentVertices(pHalfedgeMeshWrapper, m_iInteriorVertices[i]);
			double u = 0.0, v = 0.0;
			for (GEO::index_t j = 0; j < iAdjVertices.size(); j++)
			{
				if (IsBoundaryVertex(pHalfedgeMeshWrapper, iAdjVertices[j]))
				{
					double Lambda_ij = Lambda_ij_BarycentricCoordinates(pHalfedgeMeshWrapper, m_iInteriorVertices[i], iAdjVertices[j]);

					auto Iter = std::find(m_iBoundaryVertices.begin(), m_iBoundaryVertices.end(), iAdjVertices[j]);
					GEO::index_t iBoundaryIndex = GEO::index_t(Iter - m_iBoundaryVertices.begin());
					u += Lambda_ij * m_BoundaryVertices[iBoundaryIndex * 3 + 0];
					v += Lambda_ij * m_BoundaryVertices[iBoundaryIndex * 3 + 1];
				}
			}

			B[i * 2 + 0] = u;
			B[i * 2 + 1] = v;
		}

		int * ipiv = reinterpret_cast<int*>(mkl_malloc(sizeof(int) * n, 64));

		// Solve linear equation using LU decomposition
		int Info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, 2, A, n, ipiv, B, 2);

		if (Info < 0)
		{
			GEO::Logger::err("MKL") << "[" << -Info << "] param had an illegal value." << std::endl;
		}

		if (Info > 0)
		{
			GEO::Logger::err("MKL") << "U_[" << Info << ", " << Info << "] is exactly zero."
				"The factorization has been completed, but the factor U is exactly singular,"
				" so the solution could not be computed." << std::endl;
		}

		assert(Info == 0);

		m_InteriorVertices.resize(m_nInteriorVertices * 3, 0.0);
		for (GEO::index_t i = 0; i < m_nInteriorVertices; i++)
		{
			m_InteriorVertices[i * 3 + 0] = B[i * 2 + 0];
			m_InteriorVertices[i * 3 + 1] = B[i * 2 + 1];
		}
	}

	double BarycentricMappingAlgorithm::Lambda_ij_BarycentricCoordinates(
		In HalfedgeMeshWrapper * pHalfedgeMeshWrapper, 
		In GEO::index_t i, /*i = 0....n-1*/
		In GEO::index_t j
	) const
	{
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);
		assert(std::find(m_iInteriorVertices.begin(), m_iInteriorVertices.end(), i) != m_iInteriorVertices.end());

		GEO::Mesh * pMesh = pHalfedgeMeshWrapper->pMesh;

		std::function<double(GEO::index_t, GEO::index_t, GEO::index_t, GEO::index_t)> w_ij_Func;

		if (m_CoefficientType == PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[0])
		{
			w_ij_Func = std::bind(
				&BarycentricMappingAlgorithm::w_ij_MeanValueCoordinates, 
				this, 
				pHalfedgeMeshWrapper,
				std::placeholders::_1,
				std::placeholders::_2,
				std::placeholders::_3,
				std::placeholders::_4
			);
		}
		else if (m_CoefficientType == PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[1])
		{
			w_ij_Func = std::bind(
				&BarycentricMappingAlgorithm::w_ij_DiscreteHarmonicCoordinates, 
				this, pHalfedgeMeshWrapper,
				std::placeholders::_1,
				std::placeholders::_2,
				std::placeholders::_3,
				std::placeholders::_4
			);

		}
		else if (m_CoefficientType == PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[2])
		{
			w_ij_Func = std::bind(
				&BarycentricMappingAlgorithm::w_ij_WachspressCoordinates, 
				this, pHalfedgeMeshWrapper,
				std::placeholders::_1,
				std::placeholders::_2,
				std::placeholders::_3,
				std::placeholders::_4
			);
		}
		else
		{
			/*m_CoefficientType must be in PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE*/
			assert(false);
		}

		GEO::index_t iIBeginCorner = pHalfedgeMeshWrapper->Vertex2Corner[i];
		GEO::index_t iICorner = iIBeginCorner;

		while (pHalfedgeMeshWrapper->Dest(iICorner) != j)
		{
			iICorner = pHalfedgeMeshWrapper->Corner2Corner[iICorner];
			/*i and j must be adjacent*/
			assert(iICorner != iIBeginCorner);
		}

		GEO::index_t iAlpha_ij_Corner = iICorner;
		GEO::index_t iAlpha_ji_Corner = pHalfedgeMeshWrapper->Opposite(iICorner);

		double w_ij = w_ij_Func(i, j, iAlpha_ij_Corner, iAlpha_ji_Corner);
		double Sigma_k_Ni_w_ik = w_ij;

		GEO::index_t iBeginCorner = iICorner;
		GEO::index_t iCorner = pHalfedgeMeshWrapper->Corner2Corner[iICorner];
		while (iCorner != iBeginCorner)
		{
			GEO::index_t iAlpha_ik_Corner = iCorner;
			GEO::index_t iAlpha_ki_Corner = pHalfedgeMeshWrapper->Opposite(iCorner);
			GEO::index_t k = pHalfedgeMeshWrapper->Dest(iCorner);

			Sigma_k_Ni_w_ik += w_ij_Func(i, k, iAlpha_ik_Corner, iAlpha_ki_Corner);

			iCorner = pHalfedgeMeshWrapper->Corner2Corner[iCorner];
		}

		return w_ij / Sigma_k_Ni_w_ik;
	}

	double BarycentricMappingAlgorithm::w_ij_WachspressCoordinates(
		In HalfedgeMeshWrapper * pHalfedgeMeshWrapper,
		In GEO::index_t i, /*i = 0....n-1*/
		In GEO::index_t j,
		In GEO::index_t iAlpha_ij_Corner,
		In GEO::index_t iAlpha_ji_Corner
	) const
	{
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);
		assert(std::find(m_iInteriorVertices.begin(), m_iInteriorVertices.end(), i) != m_iInteriorVertices.end());
		assert(iAlpha_ij_Corner != GEO::NO_CORNER);
		assert(iAlpha_ji_Corner != GEO::NO_CORNER);

		GEO::Mesh * pMesh = pHalfedgeMeshWrapper->pMesh;

		GEO::index_t iAlpha_ij_Beta_ij_Gamma_ij_Facet = pHalfedgeMeshWrapper->Facet(iAlpha_ij_Corner);
		GEO::index_t iAlpha_ji_Beta_ji_Gamma_ji_Facet = pHalfedgeMeshWrapper->Facet(iAlpha_ji_Corner);

		GEO::index_t iAplha_ij_Vertex = i;
		GEO::index_t iBeta_ij_Vertex = j;
		GEO::index_t iGamma_ij_Vertex = pHalfedgeMeshWrapper->Origin(pHalfedgeMeshWrapper->Prev(iAlpha_ij_Corner));

		GEO::index_t iAplha_ji_Vertex = j;
		GEO::index_t iBeta_ji_Vertex = i;
		GEO::index_t iGamma_ji_Vertex = pHalfedgeMeshWrapper->Origin(pHalfedgeMeshWrapper->Prev(iAlpha_ji_Corner));

		GEO::vec3 Alpha_ij_Vertex = pMesh->vertices.point(iAplha_ij_Vertex);
		GEO::vec3 Beta_ij_Vertex = pMesh->vertices.point(iBeta_ij_Vertex);
		GEO::vec3 Gamma_ij_Vertex = pMesh->vertices.point(iGamma_ij_Vertex);

		GEO::vec3 Alpha_ji_Vertex = pMesh->vertices.point(iAplha_ji_Vertex);
		GEO::vec3 Beta_ji_Vertex = pMesh->vertices.point(iBeta_ji_Vertex);
		GEO::vec3 Gamma_ji_Vertex = pMesh->vertices.point(iGamma_ji_Vertex);

		GEO::vec3 Beta_ij_Alpha_ij_Vector = Alpha_ij_Vertex - Beta_ij_Vertex;
		GEO::vec3 Beta_ij_Gamma_ij_Vector = Gamma_ij_Vertex - Beta_ij_Vertex;

		GEO::vec3 Alpha_ji_Beta_ji_Vector = Beta_ij_Alpha_ij_Vector;
		GEO::vec3 Alpha_ji_Gamma_ji_Vector = Gamma_ji_Vertex - Alpha_ji_Vertex;

		double r_ij = Beta_ij_Alpha_ij_Vector.length();

		double Beta_ij = std::acos(GEO::dot(Beta_ij_Alpha_ij_Vector, Beta_ij_Gamma_ij_Vector) /
			(r_ij * Beta_ij_Gamma_ij_Vector.length()));
		double Alpha_ji = std::acos(GEO::dot(Alpha_ji_Beta_ji_Vector, Alpha_ji_Gamma_ji_Vector) /
			(r_ij * Alpha_ji_Gamma_ji_Vector.length()));

		double Tan_Beta_ij = std::tan(Beta_ij);
		double Tan_Alpha_ji = std::tan(Alpha_ji);

		assert(Tan_Beta_ij != 0.0 && Tan_Alpha_ji != 0.0);

		return ((1.0 / Tan_Beta_ij) + (1.0 / Tan_Alpha_ji)) / (r_ij * r_ij);
	}

	double BarycentricMappingAlgorithm::w_ij_DiscreteHarmonicCoordinates(
		In HalfedgeMeshWrapper * pHalfedgeMeshWrapper,
		In GEO::index_t i, /*i = 0....n-1*/
		In GEO::index_t j,
		In GEO::index_t iAlpha_ij_Corner,
		In GEO::index_t iAlpha_ji_Corner
	) const
	{
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);
		assert(std::find(m_iInteriorVertices.begin(), m_iInteriorVertices.end(), i) != m_iInteriorVertices.end());
		assert(iAlpha_ij_Corner != GEO::NO_CORNER);
		assert(iAlpha_ji_Corner != GEO::NO_CORNER);

		GEO::Mesh * pMesh = pHalfedgeMeshWrapper->pMesh;

		GEO::index_t iAlpha_ij_Beta_ij_Gamma_ij_Facet = pHalfedgeMeshWrapper->Facet(iAlpha_ij_Corner);
		GEO::index_t iAlpha_ji_Beta_ji_Gamma_ji_Facet = pHalfedgeMeshWrapper->Facet(iAlpha_ji_Corner);

		GEO::index_t iAplha_ij_Vertex = i;
		GEO::index_t iBeta_ij_Vertex = j;
		GEO::index_t iGamma_ij_Vertex = pHalfedgeMeshWrapper->Origin(pHalfedgeMeshWrapper->Prev(iAlpha_ij_Corner));

		GEO::index_t iAplha_ji_Vertex = j;
		GEO::index_t iBeta_ji_Vertex = i;
		GEO::index_t iGamma_ji_Vertex = pHalfedgeMeshWrapper->Origin(pHalfedgeMeshWrapper->Prev(iAlpha_ji_Corner));

		GEO::vec3 Alpha_ij_Vertex = pMesh->vertices.point(iAplha_ij_Vertex);
		GEO::vec3 Beta_ij_Vertex = pMesh->vertices.point(iBeta_ij_Vertex);
		GEO::vec3 Gamma_ij_Vertex = pMesh->vertices.point(iGamma_ij_Vertex);

		GEO::vec3 Alpha_ji_Vertex = pMesh->vertices.point(iAplha_ji_Vertex);
		GEO::vec3 Beta_ji_Vertex = pMesh->vertices.point(iBeta_ji_Vertex);
		GEO::vec3 Gamma_ji_Vertex = pMesh->vertices.point(iGamma_ji_Vertex);
		
		GEO::vec3 Gamma_ij_Aplha_ij_Vector = Alpha_ij_Vertex - Gamma_ij_Vertex;
		GEO::vec3 Gamma_ij_Beta_ij_Vector = Beta_ij_Vertex - Gamma_ij_Vertex;

		GEO::vec3 Gamma_ji_Aplha_ji_Vector = Alpha_ji_Vertex - Gamma_ji_Vertex;
		GEO::vec3 Gamma_ji_Beta_ji_Vector = Beta_ji_Vertex - Gamma_ji_Vertex;

		double Gamma_ij = std::acos(GEO::dot(Gamma_ij_Aplha_ij_Vector, Gamma_ij_Beta_ij_Vector) /
			(Gamma_ij_Aplha_ij_Vector.length() * Gamma_ij_Beta_ij_Vector.length()));
		double Gamma_ji = std::acos(GEO::dot(Gamma_ji_Aplha_ji_Vector, Gamma_ji_Beta_ji_Vector) /
			(Gamma_ji_Aplha_ji_Vector.length() * Gamma_ji_Beta_ji_Vector.length()));

		double Tan_Gamma_ij = std::tan(Gamma_ij);
		double Tan_Gamma_ji = std::tan(Gamma_ji);

		assert(Tan_Gamma_ij != 0.0 && Tan_Gamma_ji != 0.0);

		return 1.0 / Tan_Gamma_ij + 1.0 / Tan_Gamma_ji;
	}

	double BarycentricMappingAlgorithm::w_ij_MeanValueCoordinates(
		In HalfedgeMeshWrapper * pHalfedgeMeshWrapper,
		In GEO::index_t i, /*i = 0....n-1*/
		In GEO::index_t j,
		In GEO::index_t iAlpha_ij_Corner,
		In GEO::index_t iAlpha_ji_Corner
	) const
	{
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);
		assert(std::find(m_iInteriorVertices.begin(), m_iInteriorVertices.end(), i) != m_iInteriorVertices.end());
		assert(iAlpha_ij_Corner != GEO::NO_CORNER);
		assert(iAlpha_ji_Corner != GEO::NO_CORNER);

		GEO::Mesh * pMesh = pHalfedgeMeshWrapper->pMesh;

		GEO::index_t iAlpha_ij_Beta_ij_Gamma_ij_Facet = pHalfedgeMeshWrapper->Facet(iAlpha_ij_Corner);
		GEO::index_t iAlpha_ji_Beta_ji_Gamma_ji_Facet = pHalfedgeMeshWrapper->Facet(iAlpha_ji_Corner);

		GEO::index_t iAplha_ij_Vertex = i;
		GEO::index_t iBeta_ij_Vertex = j;
		GEO::index_t iGamma_ij_Vertex = pHalfedgeMeshWrapper->Origin(pHalfedgeMeshWrapper->Prev(iAlpha_ij_Corner));

		GEO::index_t iAplha_ji_Vertex = j;
		GEO::index_t iBeta_ji_Vertex = i;
		GEO::index_t iGamma_ji_Vertex = pHalfedgeMeshWrapper->Origin(pHalfedgeMeshWrapper->Prev(iAlpha_ji_Corner));

		GEO::vec3 Alpha_ij_Vertex = pMesh->vertices.point(iAplha_ij_Vertex);
		GEO::vec3 Beta_ij_Vertex = pMesh->vertices.point(iBeta_ij_Vertex);
		GEO::vec3 Gamma_ij_Vertex = pMesh->vertices.point(iGamma_ij_Vertex);

		GEO::vec3 Alpha_ji_Vertex = pMesh->vertices.point(iAplha_ji_Vertex);
		GEO::vec3 Beta_ji_Vertex = pMesh->vertices.point(iBeta_ji_Vertex);
		GEO::vec3 Gamma_ji_Vertex = pMesh->vertices.point(iGamma_ji_Vertex);

		GEO::vec3 Alpha_ij_Beta_ij_Vector = Beta_ij_Vertex - Alpha_ij_Vertex;
		GEO::vec3 Alpha_ij_Gamma_ij_Vector = Gamma_ij_Vertex - Alpha_ij_Vertex;

		GEO::vec3 Beta_ji_Alpha_ji_Vector = Alpha_ij_Beta_ij_Vector;
		GEO::vec3 Beta_ji_Gamma_ji_Vector = Gamma_ji_Vertex - Beta_ji_Vertex;

		double r_ij = Alpha_ij_Beta_ij_Vector.length();

		double Alpha_ij = std::acos(GEO::dot(Alpha_ij_Beta_ij_Vector, Alpha_ij_Gamma_ij_Vector) /
			(r_ij * Alpha_ij_Gamma_ij_Vector.length()));
		double Beta_ji = std::acos(GEO::dot(Beta_ji_Alpha_ji_Vector, Beta_ji_Gamma_ji_Vector) /
			(r_ij * Beta_ji_Gamma_ji_Vector.length()));

		double Tan_Half_Alpha_ij = std::tan(0.5 * Alpha_ij);
		double Tan_Half_Beta_ji = std::tan(0.5 * Beta_ji);

		return (Tan_Half_Alpha_ij + Tan_Half_Beta_ji) / r_ij;
	}

	void BarycentricMappingAlgorithm::GenerateOutput(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper)
	{
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);
		assert(m_InteriorVertices.size() + m_BoundaryVertices.size() == pHalfedgeMeshWrapper->pMesh->vertices.nb() * 3);

		GEO::Mesh * pMesh = pHalfedgeMeshWrapper->pMesh;

		m_ParameterizationMesh.copy(*pMesh);

		GEO::AttributesManager & CornerAttributeManager = pMesh->facet_corners.attributes();
		if (CornerAttributeManager.is_defined("tex_coord"))
		{
			CornerAttributeManager.delete_attribute_store("tex_coord");
		}

		GEO::Attribute<double> AttriVertexTexCoord(pMesh->vertices.attributes(), "tex_coord");
		AttriVertexTexCoord.redim(2);

		for (GEO::index_t i = 0; i < m_nInteriorVertices; i++)
		{
			GEO::index_t iVertex = m_iInteriorVertices[i];

			double * pVertex = m_ParameterizationMesh.vertices.point_ptr(iVertex);
			pVertex[0] = m_InteriorVertices[i * 3 + 0];
			pVertex[1] = m_InteriorVertices[i * 3 + 1];
			pVertex[2] = m_InteriorVertices[i * 3 + 2];

			AttriVertexTexCoord[iVertex * 2 + 0] = m_InteriorVertices[i * 3 + 0];
			AttriVertexTexCoord[iVertex * 2 + 1] = m_InteriorVertices[i * 3 + 1];
		}

		for (GEO::index_t i = 0; i < m_nBoundaryVertices; i++)
		{
			GEO::index_t iVertex = m_iBoundaryVertices[i];

			double * pVertex = m_ParameterizationMesh.vertices.point_ptr(iVertex);
			pVertex[0] = m_BoundaryVertices[i * 3 + 0];
			pVertex[1] = m_BoundaryVertices[i * 3 + 1];
			pVertex[2] = m_BoundaryVertices[i * 3 + 2];

			AttriVertexTexCoord[iVertex * 2 + 0] = m_BoundaryVertices[i * 3 + 0];
			AttriVertexTexCoord[iVertex * 2 + 1] = m_BoundaryVertices[i * 3 + 1];
		}
	}
}