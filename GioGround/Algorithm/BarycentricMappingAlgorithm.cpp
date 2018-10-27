#include "BarycentricMappingAlgorithm.hpp"

#include "..\HalfedgeMeshWrapper.hpp"

namespace MeshAlgorithm
{
	std::string BarycentricMappingAlgorithm::PARAMS_KEY_COEFFICIENT_TYPE = "CoefficientType";
	std::string BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[3] = 
	{
		"Mean Value Coordinates",
		"Discrete Harmonic Coordinates",
		"Wachspress Coordinates"
	};

	bool BarycentricMappingAlgorithm::Execute(InOut GEO::Mesh * pMesh)
	{
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

		return true;
	}

	bool BarycentricMappingAlgorithm::Reset()
	{
		m_CoefficientType = PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[0];
		return true;
	}
}