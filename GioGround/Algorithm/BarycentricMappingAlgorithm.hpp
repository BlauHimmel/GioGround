#pragma once

#include "IMeshAlgorithm.hpp"

namespace MeshAlgorithm
{
	// Params: [CoefficientType -> "std::string" Required]
	class BarycentricMappingAlgorithm : public IMeshAlgorithm
	{
	private:
		std::string m_CoefficientType = PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[0];

	public:
		virtual bool Execute(InOut GEO::Mesh * pMesh) override;
		virtual bool Visualize(In GEO::Mesh * const pMesh) const override;

	public:
		static std::string PARAMS_KEY_COEFFICIENT_TYPE;
		static std::string PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[3];

	protected:
		virtual bool CheckAndGetArgs(In GEO::Mesh * const pMesh) override;
		virtual bool Reset() override;
	};
}