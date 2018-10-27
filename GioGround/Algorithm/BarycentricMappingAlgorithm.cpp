#include "BarycentricMappingAlgorithm.hpp"

#include "..\HalfedgeMeshWrapper.hpp"

#include <unordered_set>

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

		GEO::Logger::out("Barycentric") << "Executing vertices permutation...";
		TIMER_START(TimePermutation);
		Permutation(&Wrapper);
		TIMER_END(TimePermutation);
		GEO::Logger::out("Barycentric") << "(" << TimePermutation << " sec)" << std::endl;

		GEO::Logger::out("Barycentric") << "Fixing boundary vertices...";
		TIMER_START(TimeFixBoundaryVertices);
		FixBoundaryVertices(&Wrapper);
		TIMER_END(TimeFixBoundaryVertices);
		GEO::Logger::out("Barycentric") << "(" << TimeFixBoundaryVertices << " sec)" << std::endl;

		TIMER_END(TimeTotal);

		GEO::Logger::out("Barycentric") << "Cut mesh finished! (Total : " << TimeTotal << " sec)" << std::endl;
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

		return true;
	}

	bool BarycentricMappingAlgorithm::Reset()
	{
		m_CoefficientType = PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[0];
		return true;
	}

	void BarycentricMappingAlgorithm::Permutation(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper)
	{
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);

		GEO::Mesh * pMesh = pHalfedgeMeshWrapper->pMesh;
		
		GEO::vector<GEO::index_t> InteriorVertices;
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
				InteriorVertices.push_back(iVertex);
			}
		}

		std::unordered_set<GEO::index_t> BoundaryCornersSet(BoundaryCornersNoOrder.begin(), BoundaryCornersNoOrder.end());
		GEO::vector<GEO::index_t> BoundaryVertices;

		GEO::index_t iCornerBegin = BoundaryCornersNoOrder[0];
		GEO::index_t iCorner = iCornerBegin;

		do
		{
			BoundaryVertices.push_back(pMesh->facet_corners.vertex(iCorner));
			
			iCorner = pHalfedgeMeshWrapper->Next(iCorner);
			while (BoundaryCornersSet.find(iCorner) == BoundaryCornersSet.end())
			{
				iCorner = pHalfedgeMeshWrapper->Corner2Corner[iCorner];
			}
		} while (iCorner != iCornerBegin);

		GEO::vector<GEO::index_t> Permutation(pMesh->vertices.nb(), GEO::index_t(-1));

		for (GEO::index_t iInteriorVertex = 0; iInteriorVertex < InteriorVertices.size(); iInteriorVertex++)
		{
			Permutation[InteriorVertices[iInteriorVertex]] = iInteriorVertex;
		}

		for (GEO::index_t iBoundaryVertex = 0; iBoundaryVertex < BoundaryVertices.size(); iBoundaryVertex++)
		{
			Permutation[BoundaryVertices[iBoundaryVertex]] = iBoundaryVertex + InteriorVertices.size();
		}

		pMesh->vertices.permute_elements(Permutation);

		m_nInteriorVertices = InteriorVertices.size();
		m_nBoundaryVertices = BoundaryVertices.size();
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
				GEO::vec3 Vertex1 = pMesh->vertices.point(i);
				GEO::vec3 Vertex2 = pMesh->vertices.point(i + 1);
				double Distance = std::sqrt(std::pow(Vertex1.x - Vertex2.x, 2.0) + std::pow(Vertex1.y - Vertex2.y, 2.0) + std::pow(Vertex1.z - Vertex2.z, 2.0));
				Total += Distance;
				Weight[i] = Distance;
			}

			GEO::vec3 Vertex1 = pMesh->vertices.point(m_nBoundaryVertices - 1);
			GEO::vec3 Vertex2 = pMesh->vertices.point(0);
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
			FixSquareBoundaryVertices(pHalfedgeMeshWrapper, Weight);
		}
		else if (m_DomainShape == PARAMS_VALUE_SUPPORTED_DOMAIN_SHAPE[1])
		{
			FixSphereBoundaryVertices(pHalfedgeMeshWrapper, Weight);
		}
	}

	void BarycentricMappingAlgorithm::FixSquareBoundaryVertices(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper, In GEO::vector<double> & Weight)
	{
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);

		double MaxX = 1.0, MinX = 0.0;
		double MaxY = 1.0, MinY = 0.0;
		double X = MinX, Y = MinY, Z = 0.0;
		double BoundaryLength = (MaxX - MinX + MaxY - MinY) * 2.0;
		GEO::index_t iDirection = 0;

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
					X = MaxX;
					Y += Step - (MaxX - X);
					iDirection++;
				}
				else
				{
					X += Step;
				}
			}
			else if (iDirection == 1)
			{
				if (Y + Step > MaxY)
				{
					Y = MaxY;
					X -= Step - (MaxY - Y);
					iDirection++;
				}
				else
				{
					Y += Step;
				}
			}
			else if (iDirection == 2)
			{
				if (X - Step < MinX)
				{
					X = MinX;
					Y -= Step - (X - MinX);
					iDirection++;
				}
				else
				{
					X -= Step;
				}
			}
			else if (iDirection == 3)
			{
				if (Y - Step < MinY)
				{
					Y = MinY;
					X += Step - (Y - MinY);
					iDirection++;
				}
				else
				{
					Y -= Step;
				}
			}
		}
	}

	void BarycentricMappingAlgorithm::FixSphereBoundaryVertices(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper, In GEO::vector<double> & Weight)
	{
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);

		double Radius = 0.5;
		double Theta = 0.0;
		double Parameters = 2.0 * M_PI;

		m_BoundaryVertices.resize(m_nBoundaryVertices * 3, 0.0);

		for (GEO::index_t i = 0; i < Weight.size(); i++)
		{
			double X = Radius * std::cos(Theta);
			double Y = Radius * std::sin(Theta);
			double Z = 0.0;

			m_BoundaryVertices[3 * i + 0] = X;
			m_BoundaryVertices[3 * i + 1] = Y;
			m_BoundaryVertices[3 * i + 2] = Z;

			double Step = Parameters * Weight[i];
			Theta += Step;
		}
	}
}