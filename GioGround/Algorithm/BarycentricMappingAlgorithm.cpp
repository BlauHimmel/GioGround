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

		GEO::Logger::out("Barycentric Mapping") << "Create Halfedge data structure...";
		TIMER_START(TimeCreateHalfedge);
		HalfedgeMeshWrapper Wrapper(pMesh);
		TIMER_END(TimeCreateHalfedge);
		GEO::Logger::out("Barycentric Mapping") << "(" << TimeCreateHalfedge << " sec)" << std::endl;

		GEO::Logger::out("Barycentric Mapping") << "Executing vertices permutation...";
		TIMER_START(TimePermutation);
		Permutation(&Wrapper);
		TIMER_END(TimePermutation);
		GEO::Logger::out("Barycentric Mapping") << "(" << TimePermutation << " sec)" << std::endl;

		TIMER_END(TimeTotal);

		GEO::Logger::out("Barycentric Mapping") << "Cut mesh finished! (Total : " << TimeTotal << " sec)" << std::endl;
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
		GEO::vector<GEO::index_t> BoundaryVerticesNoOrder;

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
				BoundaryVerticesNoOrder.push_back(iVertex);
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
	}

	void BarycentricMappingAlgorithm::FixBoundaryVertices(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper)
	{
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);

		if (m_DomainShape == PARAMS_VALUE_SUPPORTED_DOMAIN_SHAPE[0])
		{
			FixSquareBoundaryVertices(pHalfedgeMeshWrapper);
		}
		else if (m_DomainShape == PARAMS_VALUE_SUPPORTED_DOMAIN_SHAPE[1])
		{
			FixSphereBoundaryVertices(pHalfedgeMeshWrapper);
		}
	}

	void BarycentricMappingAlgorithm::FixSquareBoundaryVertices(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper)
	{
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);

	}

	void BarycentricMappingAlgorithm::FixSphereBoundaryVertices(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper)
	{
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);

	}
}