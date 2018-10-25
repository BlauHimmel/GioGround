#include "MeshCutAlgorithm.hpp"

#include "..\HalfedgeMeshWrapper.hpp"

#include <list>
#include <queue>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <assert.h>

namespace MeshAlgorithm
{
	std::string MeshCutAlgorithm::PARAMS_KEY_START_FACET = "StartFacet";

	bool MeshCutAlgorithm::Execute(InOut GEO::Mesh * pMesh)
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

		GEO::Logger::out("MeshCutAlgorithm") << "Create Halfedge data structure...";
		TIMER_START(TimeCreateHalfedge);
		HalfedgeMeshWrapper Wrapper(pMesh);
		TIMER_END(TimeCreateHalfedge);
		GEO::Logger::out("MeshCutAlgorithm") << "(" << TimeCreateHalfedge << " sec)" << std::endl;

		GEO::Logger::out("MeshCutAlgorithm") << "Marking facets of mesh using BFS...";
		TIMER_START(TimeMarkMeshByBFS);
		MarkMeshByBFS(&Wrapper);
		TIMER_END(TimeMarkMeshByBFS);
		GEO::Logger::out("MeshCutAlgorithm") << "(" << TimeMarkMeshByBFS << " sec)" << std::endl;

		GEO::Logger::out("MeshCutAlgorithm") << "Cutting branch of unmarked corners...";
		TIMER_START(TimeCutBranchEdge);
		CutBranchEdge(&Wrapper);
		TIMER_END(TimeCutBranchEdge);
		GEO::Logger::out("MeshCutAlgorithm") << "(" << TimeCutBranchEdge << " sec)" << std::endl;

		GEO::Logger::out("MeshCutAlgorithm") << "Cutting mesh by left edges...";
		TIMER_START(TimeCutMeshByEdge);
		CutMeshByEdge(&Wrapper);
		TIMER_END(TimeCutMeshByEdge);
		GEO::Logger::out("MeshCutAlgorithm") << "(" << TimeCutMeshByEdge << " sec)" << std::endl;

		GEO::Logger::out("MeshCutAlgorithm") << "Computing the points of cutting edge used for visualization...";
		TIMER_START(TimeComputeCuttingEdgePoints);
		ComputeCuttingEdgePoints(pMesh);
		TIMER_END(TimeComputeCuttingEdgePoints);
		GEO::Logger::out("MeshCutAlgorithm") << "(" << TimeComputeCuttingEdgePoints << " sec)" << std::endl;

		TIMER_END(TimeTotal);

		GEO::Logger::out("MeshCutAlgorithm") << "Cut mesh finished! (Total : " << TimeTotal << " sec)" << std::endl;
		return true;
	}

	bool MeshCutAlgorithm::Visualize(In GEO::Mesh * const pMesh) const
	{
		if (pMesh == nullptr)
		{
			return false;
		}
		
		GEO::index_t EdgeWidth = 3;
		GEO::vec4f EdgeColor(0.0f, 0.0f, 1.0f, 1.0f);

		glupSetMeshWidth(EdgeWidth);
		glupSetColor4fv(GLUP_FRONT_AND_BACK_COLOR, EdgeColor.data());
		glupBegin(GLUP_LINES);

		for (GEO::index_t iElement = 0; iElement < m_CuttingEdgePoints.size(); iElement += 3)
		{
			glupVertex3d(m_CuttingEdgePoints[iElement + 0], m_CuttingEdgePoints[iElement + 1], m_CuttingEdgePoints[iElement + 2]);
		}

		glupEnd();
		return true;
	}

	bool MeshCutAlgorithm::CheckAndGetArgs(In GEO::Mesh * const pMesh)
	{
		pMesh->assert_is_valid();
		assert(pMesh->vertices.double_precision());

		if (!GetArg(PARAMS_KEY_START_FACET, &m_StartFacet))
		{
			return false;
		}

		if (m_StartFacet < 0 || m_StartFacet >= pMesh->facets.nb())
		{
			return false;
		}

		return true;
	}

	bool MeshCutAlgorithm::Reset()
	{
		m_StartFacet = GEO::NO_FACET;
		m_NotMarkedCornersIdx.clear();

		m_CuttingEdgePoints.clear();
		m_CuttingEdgePointsIdx.clear();

		return true;
	}

	void MeshCutAlgorithm::MarkMeshByBFS(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper)
	{
		assert(m_StartFacet != GEO::NO_FACET);
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);

		GEO::Mesh * pMesh = pHalfedgeMeshWrapper->pMesh;

		GEO::Attribute<bool> AttriIsFacetVisited;
		GEO::Attribute<bool> AttriIsCornerMarked;

		AttriIsFacetVisited.bind(pMesh->facets.attributes(), "IsFacetVisited");
		AttriIsCornerMarked.bind(pMesh->facet_corners.attributes(), "IsCornerMarked");

		AttriIsFacetVisited.fill(false);
		AttriIsCornerMarked.fill(false);

		std::queue<GEO::index_t/*iFacet*/> Queue;
		Queue.push(m_StartFacet);

		while (!Queue.empty())
		{
			GEO::index_t iFacet = Queue.front();
			Queue.pop();

			for (GEO::index_t iCorner = pMesh->facets.corners_begin(iFacet); iCorner != pMesh->facets.corners_end(iFacet); iCorner++)
			{
				GEO::index_t iOppsiteCorner = pHalfedgeMeshWrapper->Opposite(iCorner);
				if (iOppsiteCorner == GEO::NO_CORNER)
				{
					continue;
				}

				GEO::index_t iAdjFacet = pHalfedgeMeshWrapper->Facet(iOppsiteCorner);
				if (AttriIsFacetVisited[iAdjFacet])
				{
					continue;
				}

				AttriIsFacetVisited[iAdjFacet] = true;
				AttriIsCornerMarked[iCorner] = true;
				AttriIsCornerMarked[iOppsiteCorner] = true;

				Queue.push(iAdjFacet);
			}
		}

		for (GEO::index_t iCorner = 0; iCorner < pMesh->facet_corners.nb(); iCorner++)
		{
			if (!AttriIsCornerMarked[iCorner])
			{
				m_NotMarkedCornersIdx.push_back(iCorner);
			}
		}
	}

	void MeshCutAlgorithm::CutBranchEdge(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper)
	{
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);

		GEO::Mesh * pMesh = pHalfedgeMeshWrapper->pMesh;

		GEO::Attribute<bool> AttriIsCornerMarked(pMesh->facet_corners.attributes(), "IsCornerMarked");

		bool bIsFindBranch = true;
		while (bIsFindBranch)
		{
			bIsFindBranch = false;

			for (GEO::index_t i = 0; i < m_NotMarkedCornersIdx.size(); i++)
			{
				GEO::index_t iCornerBegin = m_NotMarkedCornersIdx[i];

				if (AttriIsCornerMarked[iCornerBegin])
				{
					continue;
				}

				GEO::index_t iCorner = iCornerBegin;
				GEO::index_t nOutDegree = 0;
				GEO::index_t iCuttingCorner = GEO::NO_CORNER;
				bool bIsBoundary = false;

				do
				{
					nOutDegree = 0;
					iCuttingCorner = GEO::NO_CORNER;

					do
					{
						if (!AttriIsCornerMarked[iCorner])
						{
							nOutDegree++;
							iCuttingCorner = iCorner;
							if (nOutDegree > 1)
							{
								break;
							}
						}
						iCorner = pHalfedgeMeshWrapper->NextAroundVertex(iCorner);
					} while (iCorner != iCornerBegin && iCorner != GEO::NO_CORNER);

					if (nOutDegree == 1)
					{
						GEO::index_t iOppsiteCorner = pHalfedgeMeshWrapper->Opposite(iCuttingCorner);

						/*iOppsiteCorner == GEO::NO_CORNER while iCorner at boundary*/
						if (iOppsiteCorner != GEO::NO_CORNER)
						{
							bIsFindBranch = true;
							bIsBoundary = false;

							AttriIsCornerMarked[iCuttingCorner] = true;
							AttriIsCornerMarked[iOppsiteCorner] = true;
						}
						else
						{
							bIsBoundary = true;
						}

						iCornerBegin = pHalfedgeMeshWrapper->Next(iCuttingCorner);
						iCorner = iCornerBegin;
					}

				} while (nOutDegree == 1 && !bIsBoundary);
			}
		}

		/*Update the m_NotMarkedCornersIdx after cutting branch*/
		m_NotMarkedCornersIdx.clear();
		for (GEO::index_t iCorner = 0; iCorner < pMesh->facet_corners.nb(); iCorner++)
		{
			if (!AttriIsCornerMarked[iCorner])
			{
				m_NotMarkedCornersIdx.push_back(iCorner);
			}
		}
	}

	void MeshCutAlgorithm::FindLoopByDFS(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper)
	{
		GEO::Mesh * pMesh = pHalfedgeMeshWrapper->pMesh;

		GEO::Attribute<bool> AttriIsCornerMarked(pMesh->facet_corners.attributes(), "IsCornerMarked");
		std::unordered_set<GEO::index_t> CuttingCornerIdx;

		for (GEO::index_t i = 0; i < m_NotMarkedCornersIdx.size(); i++)
		{
			std::list<GEO::index_t/*iCorner*/> Stack;
			std::vector<bool> CornerFlag(pMesh->facet_corners.nb(), false);

			GEO::index_t iCorner = m_NotMarkedCornersIdx[i];
			Stack.push_front(iCorner);

			GEO::index_t iOppositeCorner = pHalfedgeMeshWrapper->Opposite(iCorner);
			CornerFlag[iCorner] = true;
			if (iOppositeCorner != GEO::NO_FACET)
			{
				CornerFlag[iOppositeCorner] = true;
			}

			while (!Stack.empty())
			{
				iCorner = Stack.front();

				GEO::index_t iCornerBegin = pHalfedgeMeshWrapper->Next(iCorner);
				iCorner = iCornerBegin;
				bool bNoAdjFound = true;
				do
				{
					if (iCorner == m_NotMarkedCornersIdx[i]) // Loop find
					{
						for (auto Iter = Stack.begin(); Iter != Stack.end(); ++Iter)
						{
							CuttingCornerIdx.insert(*Iter);
						}
						break;
					}

					if (!CornerFlag[iCorner] && std::find(m_NotMarkedCornersIdx.begin(), m_NotMarkedCornersIdx.end(), iCorner) != m_NotMarkedCornersIdx.end())
					{
						Stack.push_front(iCorner);

						iOppositeCorner = pHalfedgeMeshWrapper->Opposite(iCorner);
						CornerFlag[iCorner] = true;
						if (iOppositeCorner != GEO::NO_FACET)
						{
							CornerFlag[iOppositeCorner] = true;
						}

						bNoAdjFound = false;
						break;
					}

					iOppositeCorner = pHalfedgeMeshWrapper->Opposite(iCorner);
					if (iOppositeCorner == GEO::NO_FACET)
					{
						break;
					}

					iCorner = pHalfedgeMeshWrapper->Next(iOppositeCorner);
				} while (iCorner != iCornerBegin && iCorner != GEO::NO_FACET);

				if (bNoAdjFound)
				{
					Stack.pop_front();
				}
			}
		}

		for (GEO::index_t i = 0; i < m_NotMarkedCornersIdx.size(); i++)
		{
			GEO::index_t iCorner = m_NotMarkedCornersIdx[i];
			if (CuttingCornerIdx.find(iCorner) == CuttingCornerIdx.end())
			{
				AttriIsCornerMarked[iCorner] = true;
			}
		}

		/*Update the m_NotMarkedCornersIdx after cutting branch*/
		m_NotMarkedCornersIdx.clear();
		for (GEO::index_t iCorner = 0; iCorner < pMesh->facet_corners.nb(); iCorner++)
		{
			if (!AttriIsCornerMarked[iCorner])
			{
				m_NotMarkedCornersIdx.push_back(iCorner);
			}
		}
	}

	void MeshCutAlgorithm::CutMeshByEdge(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper)
	{
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);

		GEO::Mesh * pMesh = pHalfedgeMeshWrapper->pMesh;

		GEO::Attribute<bool> AttriIsCornerMarked(pMesh->facet_corners.attributes(), "IsCornerMarked");

		for (GEO::index_t i = 0; i < m_NotMarkedCornersIdx.size(); i++)
		{
			GEO::index_t iBeginCorner = m_NotMarkedCornersIdx[i];

			if (AttriIsCornerMarked[iBeginCorner])
			{
				continue;
			}
			
			GEO::index_t iCorner = iBeginCorner;
			GEO::index_t iNewVertex = pMesh->vertices.create_vertex();
			GEO::index_t iOldVertex = pMesh->facet_corners.vertex(iCorner);
			GEO::index_t iCuttingCorner = iBeginCorner;

			double * pOldVertex = pMesh->vertices.point_ptr(iOldVertex);
			double * pNewVertex = pMesh->vertices.point_ptr(iNewVertex);

			pNewVertex[0] = pOldVertex[0];
			pNewVertex[1] = pOldVertex[1];
			pNewVertex[2] = pOldVertex[2];

			bool bIsCirculateEnd = false;

			///
			int nNewVertex = 1;
			int nCirculateStep = 0;
			///

			do
			{
				GEO::index_t iFacet = pHalfedgeMeshWrapper->Facet(iCorner);
				
				///
				auto nV = pMesh->facets.nb_vertices(iFacet);
				auto iV1 = pMesh->facets.vertex(iFacet, 0);
				auto iV2 = pMesh->facets.vertex(iFacet, 1);
				auto iV3 = pMesh->facets.vertex(iFacet, 2);
				nCirculateStep++;
				///

				GEO::index_t iLocalVertex = pMesh->facets.find_vertex(iFacet, iOldVertex);
				pMesh->facets.set_vertex(iFacet, iLocalVertex, iNewVertex);

				iCorner = pHalfedgeMeshWrapper->NextAroundVertex(iCorner);

				bIsCirculateEnd = (iCorner == GEO::NO_CORNER || iCorner == iBeginCorner);

				if (!bIsCirculateEnd && !AttriIsCornerMarked[iCorner])
				{
					iNewVertex = pMesh->vertices.create_vertex();
					pNewVertex = pMesh->vertices.point_ptr(iNewVertex);

					pNewVertex[0] = pOldVertex[0];
					pNewVertex[1] = pOldVertex[1];
					pNewVertex[2] = pOldVertex[2];

					AttriIsCornerMarked[iCuttingCorner] = true;

					iCuttingCorner = iCorner;

					///
					nNewVertex++;
					///
				}

				/*Change the tag of last cutting edge when circulation was finished.*/
				if (bIsCirculateEnd)
				{
					AttriIsCornerMarked[iCuttingCorner] = true;
				}

			} while (!bIsCirculateEnd);

			if (iCorner == GEO::NO_CORNER && pHalfedgeMeshWrapper->Opposite(iBeginCorner) != GEO::NO_CORNER)
			{
				iCorner = pHalfedgeMeshWrapper->Corner2Corner[iBeginCorner];

				do 
				{
					GEO::index_t iFacet = pHalfedgeMeshWrapper->Facet(iCorner);
					GEO::index_t iLocalVertex = pMesh->facets.find_vertex(iFacet, iOldVertex);
					pMesh->facets.set_vertex(iFacet, iLocalVertex, iNewVertex);

					if (!AttriIsCornerMarked[iCorner])
					{
						iNewVertex = pMesh->vertices.create_vertex();
						pNewVertex = pMesh->vertices.point_ptr(iNewVertex);

						pNewVertex[0] = pOldVertex[0];
						pNewVertex[1] = pOldVertex[1];
						pNewVertex[2] = pOldVertex[2];

						AttriIsCornerMarked[iCorner] = true;
					}

					iCorner = pHalfedgeMeshWrapper->Corner2Corner[iCorner];
				} while (iCorner != GEO::NO_CORNER);
			}
		}

		pMesh->vertices.remove_isolated();
	}

	void MeshCutAlgorithm::ComputeCuttingEdgePoints(In GEO::Mesh * const pMesh)
	{
		GEO::Mesh CopyMesh;
		CopyMesh.copy(*pMesh);

		HalfedgeMeshWrapper Wrapper(&CopyMesh);
		GEO::Attribute<bool> AttriIsCornerMarked(CopyMesh.facet_corners.attributes(), "IsCornerMarked");

		for (GEO::index_t iCorner = 0; iCorner < CopyMesh.facet_corners.nb(); iCorner++)
		{
			if (!AttriIsCornerMarked[iCorner])
			{
				AttriIsCornerMarked[iCorner] = true;

				GEO::index_t iOppositeCorner = Wrapper.Opposite(iCorner);
				if (iOppositeCorner == GEO::NO_CORNER)
				{
					iOppositeCorner = Wrapper.Next(iCorner);
				}
				else
				{
					AttriIsCornerMarked[iOppositeCorner] = true;
				}

				GEO::index_t iOrigin = Wrapper.Origin(iCorner);
				GEO::index_t iDest = Wrapper.Dest(iCorner);

				GEO::vec3 vOrigin = CopyMesh.vertices.point(iOrigin);
				GEO::vec3 vDest = CopyMesh.vertices.point(iDest);

				m_CuttingEdgePoints.insert(m_CuttingEdgePoints.end(), { vOrigin[0], vOrigin[1], vOrigin[2] });
				m_CuttingEdgePoints.insert(m_CuttingEdgePoints.end(), { vDest[0], vDest[1], vDest[2] });

				m_CuttingEdgePointsIdx.insert(m_CuttingEdgePointsIdx.end(), { iOrigin, iDest });
			}
		}
	}
}