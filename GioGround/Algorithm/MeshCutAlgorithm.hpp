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

		GEO::vector<float> m_CuttingEdgePointsF;
		GEO::vector<double> m_CuttingEdgePointsD;
		GEO::vector<GEO::index_t> m_CuttingEdgePointsIdx;

	public:
		static std::string PARAMS_KEY_START_FACET;

	public:
		virtual bool Execute(InOut GEO::Mesh * pMesh) override;
		virtual bool Visualize(In GEO::Mesh * const pMesh) const override;

	protected:
		virtual bool CheckAndGetArgs(In GEO::Mesh * const pMesh) override;

	private:
		void MarkMeshByBFS(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper);
		void CutBranchEdge(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper);
		
		void CutMeshByEdge(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper);
		void FindLoopByDFS(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper);

		void ComputeCuttingEdgePoints(In GEO::Mesh * const pMesh);
	};
}