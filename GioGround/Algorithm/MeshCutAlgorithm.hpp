#pragma once

#include "IMeshAlgorithm.hpp"

namespace MeshAlgorithm
{
	// Params: [StartFacet -> GEO::index_t Required]
	class MeshCutAlgorithm : public IMeshAlgorithm
	{
	private:
		GEO::index_t m_StartFacet = GEO::NO_FACET;
		GEO::vector<GEO::index_t> m_NotMarkedCornersIdx;

		GEO::vector<double> m_CuttingEdgePoints;
		GEO::vector<GEO::index_t> m_CuttingEdgePointsIdx;

	public:
		static std::string PARAMS_KEY_START_FACET;

	public:
		virtual bool Execute(InOut GEO::Mesh * pMesh) override;
		virtual bool Visualize(In GEO::Mesh * const pMesh) const override;

	protected:
		virtual bool CheckAndGetArgs(In GEO::Mesh * const pMesh) override;
		virtual bool Reset() override;

	private:
		// Step1
		void MarkMeshByBFS(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper);

		// Step2
		void CutBranchEdge(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper);
		[[deprecated("This function has been deprecated for the reason of efficiency")]]
		void FindLoopByDFS(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper);
		
		// Step3
		void ComputeCuttingEdgePoints(In GEO::Mesh * const pMesh);

		// Step4
		void CutMeshByEdge(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper);

		// Step5
		void ClearAttribute(In GEO::Mesh * const pMesh);
	};
}