#pragma once

#include "IMeshAlgorithm.hpp"

namespace MeshAlgorithm
{
	// Params: [CoefficientType -> "std::string" Required]
	class BarycentricMappingAlgorithm : public IMeshAlgorithm
	{
	private:
		std::string m_CoefficientType = PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[0];
		std::string m_DomainShape = PARAMS_VALUE_SUPPORTED_DOMAIN_SHAPE[0];

	public:
		virtual bool Execute(InOut GEO::Mesh * pMesh) override;
		virtual bool Visualize(In GEO::Mesh * const pMesh) const override;

	public:
		static std::string PARAMS_KEY_COEFFICIENT_TYPE;
		static std::string PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[3];

		static std::string PARAMS_KEY_DOMAIN_SHAPE;
		static std::string PARAMS_VALUE_SUPPORTED_DOMAIN_SHAPE[2];

		static std::string PARAMS_KEY_BOUNDARY_FIX_WEIGHT;
		static std::string PARAMS_VALUE_SUPPORTED_BOUNDARY_FIX_WEIGHT[2];

	protected:
		virtual bool CheckAndGetArgs(In GEO::Mesh * const pMesh) override;
		virtual bool Reset() override;

	private:
		// Step1: 1....n interior vertices, n+1....n+b boundary vertices;
		void Permutation(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper);

		// Step2
		void FixBoundaryVertices(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper);
		void FixSquareBoundaryVertices(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper);
		void FixSphereBoundaryVertices(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper);

	};
}