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

		GEO::Logger::out("Mesh Cut") << "Create Halfedge data structure...";
		TIMER_START(TimeCreateHalfedge);
		HalfedgeMeshWrapper Wrapper(pMesh);
		TIMER_END(TimeCreateHalfedge);
		GEO::Logger::out("Mesh Cut") << "(" << TimeCreateHalfedge << " sec)" << std::endl;

		GEO::Logger::out("Mesh Cut") << "Marking facets of mesh using BFS...";
		TIMER_START(TimeMarkMeshByBFS);
		MarkMeshByBFS(&Wrapper);
		TIMER_END(TimeMarkMeshByBFS);
		GEO::Logger::out("Mesh Cut") << "(" << TimeMarkMeshByBFS << " sec)" << std::endl;

		GEO::Logger::out("Mesh Cut") << "Cutting branch of unmarked corners...";
		TIMER_START(TimeCutBranchEdge);
		CutBranchEdge(&Wrapper);
		TIMER_END(TimeCutBranchEdge);
		GEO::Logger::out("Mesh Cut") << "(" << TimeCutBranchEdge << " sec)" << std::endl;

		GEO::Logger::out("Mesh Cut") << "Connecting boundaries of mesh...";
		TIMER_START(TimeConnectBoundary);
		ConnectBoundary(&Wrapper);
		TIMER_END(TimeConnectBoundary);
		GEO::Logger::out("Mesh Cut") << "(" << TimeConnectBoundary << " sec)" << std::endl;

		GEO::Logger::out("Mesh Cut") << "Computing the points of cutting edge used for visualization...";
		TIMER_START(TimeComputeCuttingEdgePoints);
		ComputeCuttingEdgePoints(pMesh);
		TIMER_END(TimeComputeCuttingEdgePoints);
		GEO::Logger::out("Mesh Cut") << "(" << TimeComputeCuttingEdgePoints << " sec)" << std::endl;

		GEO::Logger::out("Mesh Cut") << "Cutting mesh by left edges...";
		TIMER_START(TimeCutMeshByEdge);
		CutMeshByEdge(&Wrapper);
		TIMER_END(TimeCutMeshByEdge);
		GEO::Logger::out("Mesh Cut") << "(" << TimeCutMeshByEdge << " sec)" << std::endl;

		GEO::Logger::out("Mesh Cut") << "Clearing attributes...";
		TIMER_START(TimeClearAttribute);
		ClearAttribute(pMesh);
		TIMER_END(TimeClearAttribute);
		GEO::Logger::out("Mesh Cut") << "(" << TimeClearAttribute << " sec)" << std::endl;

		TIMER_END(TimeTotal);

		GEO::Logger::out("Mesh Cut") << "Cut mesh finished! (Total : " << TimeTotal << " sec)" << std::endl;
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
	}

	void MeshCutAlgorithm::CutBranchEdge(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper)
	{
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);

		GEO::Mesh * pMesh = pHalfedgeMeshWrapper->pMesh;

		GEO::Attribute<bool> AttriIsCornerMarked(pMesh->facet_corners.attributes(), "IsCornerMarked");

		/*Update the m_NotMarkedCornersIdx after cutting branch*/
		m_NotMarkedCornersIdx.clear();
		for (GEO::index_t iCorner = 0; iCorner < pMesh->facet_corners.nb(); iCorner++)
		{
			if (!AttriIsCornerMarked[iCorner])
			{
				m_NotMarkedCornersIdx.push_back(iCorner);
			}
		}

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
	}

	void MeshCutAlgorithm::FindLoopByDFS(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper)
	{
		GEO::Mesh * pMesh = pHalfedgeMeshWrapper->pMesh;
		
		GEO::Attribute<bool> AttriIsCornerMarked(pMesh->facet_corners.attributes(), "IsCornerMarked");

		/*Update the m_NotMarkedCornersIdx after cutting branch*/
		m_NotMarkedCornersIdx.clear();
		for (GEO::index_t iCorner = 0; iCorner < pMesh->facet_corners.nb(); iCorner++)
		{
			if (!AttriIsCornerMarked[iCorner])
			{
				m_NotMarkedCornersIdx.push_back(iCorner);
			}
		}

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
	}

	void MeshCutAlgorithm::ConnectBoundary(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper)
	{
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);

		GEO::Mesh * pMesh = pHalfedgeMeshWrapper->pMesh;

		GEO::Attribute<GEO::index_t> AttriBoundaryIndex;
		AttriBoundaryIndex.bind(pMesh->facet_corners.attributes(), "BoundaryIndex");
		AttriBoundaryIndex.fill(GEO::index_t(-1));

		GEO::Attribute<bool> AttriIsCornerMarked(pMesh->facet_corners.attributes(), "IsCornerMarked");

		/*Update the m_NotMarkedCornersIdx after cutting branch*/
		m_NotMarkedCornersIdx.clear();
		for (GEO::index_t iCorner = 0; iCorner < pMesh->facet_corners.nb(); iCorner++)
		{
			if (!AttriIsCornerMarked[iCorner])
			{
				m_NotMarkedCornersIdx.push_back(iCorner);
			}
		}

		GEO::vector<GEO::vector<GEO::index_t>> BoundariesCorners;
		BoundariesCorners.push_back(GEO::vector<GEO::index_t>());

		/*Index the boundary, 0 for we generated, others for the original boundaries*/
		bool bOnlyOriginalBoundary = true;
		GEO::index_t iBoundary = 1;
		for (GEO::index_t i = 0; i < m_NotMarkedCornersIdx.size(); i++)
		{
			GEO::index_t iCorner = m_NotMarkedCornersIdx[i];
			GEO::index_t iOppositeIndex = pHalfedgeMeshWrapper->Opposite(iCorner);
			if (iOppositeIndex != GEO::NO_CORNER)
			{
				AttriBoundaryIndex[iCorner] = 0;
				AttriBoundaryIndex[iOppositeIndex] = 0;
				BoundariesCorners[GEO::index_t(0)].insert(BoundariesCorners[GEO::index_t(0)].end(), { iCorner, iOppositeIndex });
				bOnlyOriginalBoundary = false;
				continue;
			}

			if (AttriBoundaryIndex[iCorner] != GEO::index_t(-1))
			{
				continue;
			}

			GEO::vector<GEO::index_t> BoundaryCorners;

			do 
			{
				AttriBoundaryIndex[iCorner] = iBoundary;
				BoundaryCorners.push_back(iCorner);

				GEO::index_t iNextCorner = pHalfedgeMeshWrapper->Next(iCorner);
				iCorner = iNextCorner;
				while (AttriIsCornerMarked[iCorner])
				{
					iCorner = pHalfedgeMeshWrapper->Corner2Corner[iCorner];
				}

			} while (AttriBoundaryIndex[iCorner] == GEO::index_t(-1));

			BoundariesCorners.push_back(BoundaryCorners);
			iBoundary++;
		}

		/*Find the shortest path by BFS*/
		GEO::index_t iMainBoundary = 0;
		if (bOnlyOriginalBoundary)
		{
			GEO::index_t nMaxSize = 0;
			for (GEO::index_t i = 0; i < BoundariesCorners.size(); i++)
			{
				if (nMaxSize < BoundariesCorners[i].size())
				{
					nMaxSize = BoundariesCorners[i].size();
					iMainBoundary = i;
				}
			}
		}

		for (GEO::index_t i = 0; i < BoundariesCorners.size(); i++)
		{
			if (i == iMainBoundary || BoundariesCorners[i].size() == 0)
			{
				continue;
			}

			std::queue<GEO::index_t/*iFacet*/> Queue;
			GEO::Attribute<GEO::index_t> AttriFacetDistance;
			AttriFacetDistance.bind(pMesh->facets.attributes(), "FacetDistance");
			AttriFacetDistance.fill(GEO::index_t(-1));
			GEO::index_t Distance = 1;

			for (GEO::index_t j = 0; j < BoundariesCorners[i].size(); j++)
			{
				GEO::index_t iFacet = pHalfedgeMeshWrapper->Facet(BoundariesCorners[i][j]);
				AttriFacetDistance[iFacet] = Distance;
				Queue.push(iFacet);
			}

			Queue.push(GEO::index_t(-1));

			GEO::index_t iFoundFacet = GEO::NO_FACET;
			GEO::index_t NearDistance = GEO::index_t(-1);

			while (!Queue.empty())
			{
				GEO::index_t iFacet = Queue.front();
				Queue.pop();
				if (iFacet == GEO::index_t(-1) && !Queue.empty())
				{
					Distance++;
					Queue.push(GEO::index_t(-1));
					continue;
				}

				for (GEO::index_t j = 0; j < pMesh->facets.nb_vertices(iFacet); j++)
				{
					GEO::index_t iCorner = pMesh->facets.corner(iFacet, j);
					if (AttriBoundaryIndex[iCorner] == iMainBoundary)
					{
						iFoundFacet = iFacet;
						NearDistance = Distance;
						goto JumpOut; // =====> GOTO
					}
				}

				for (GEO::index_t j = 0; j < pMesh->facets.nb_vertices(iFacet); j++)
				{
					GEO::index_t iAdjFacet = pMesh->facets.adjacent(iFacet, j);
					if (iAdjFacet != GEO::NO_FACET && AttriFacetDistance[iAdjFacet] == GEO::index_t(-1))
					{
						Queue.push(iAdjFacet);
						AttriFacetDistance[iAdjFacet] = Distance + 1;
					}
				}
			}

			JumpOut: // <===== GOTO

			GEO::index_t iFacet = iFoundFacet;
			GEO::index_t CurrentDistance = NearDistance;

			GEO::vector<GEO::index_t> CornerPath;

			/*Trace back to connect the boundary*/
			do 
			{
				for (GEO::index_t j = 0; j < pMesh->facets.nb_vertices(iFacet); j++)
				{
					GEO::index_t iAdjFacet = pMesh->facets.adjacent(iFacet, j);

					if (iAdjFacet != GEO::NO_FACET && AttriFacetDistance[iAdjFacet] == CurrentDistance - 1)
					{
						GEO::index_t iCorner = pMesh->facets.corner(iFacet, j);
						GEO::index_t iOppositeCorner = pHalfedgeMeshWrapper->Opposite(iCorner);

						AttriIsCornerMarked[iCorner] = false;
						AttriIsCornerMarked[iOppositeCorner] = false;

						CornerPath.insert(CornerPath.end(), { iCorner, iOppositeCorner });
						
						iFacet = iAdjFacet;
						CurrentDistance--;

						break;
					}
				}
			} while (AttriFacetDistance[iFacet] != 1);

			for (GEO::index_t j = 0; j < CornerPath.size(); j++)
			{
				GEO::index_t iBeginCorner = CornerPath[j];
				GEO::index_t iOppositeCorner = (j % 2 == 0) ? CornerPath[j + 1] : CornerPath[j - 1];
				GEO::index_t iCorner = iBeginCorner;

				if (AttriIsCornerMarked[iBeginCorner] || AttriIsCornerMarked[iOppositeCorner])
				{
					continue;
				}

				GEO::index_t nOutDegree = 0;

				do
				{
					if (!AttriIsCornerMarked[iCorner])
					{
						nOutDegree++;
						if (nOutDegree > 1)
						{
							break;
						}
					}
					iCorner = pHalfedgeMeshWrapper->Corner2Corner[iCorner];
				} while (iCorner != iBeginCorner);

				if (nOutDegree == 1)
				{
					AttriIsCornerMarked[iBeginCorner] = true;
					AttriIsCornerMarked[iOppositeCorner] = true;
				}

				GEO::index_t iFacet = pHalfedgeMeshWrapper->Facet(iBeginCorner);
				bool bIsAllEdgeMarked = true;

				for (GEO::index_t k = 0; k < pMesh->facets.nb_corners(iFacet); k++)
				{
					if (AttriIsCornerMarked[pMesh->facets.corner(iFacet, k)])
					{
						bIsAllEdgeMarked = false;
						break;
					}
				}

				if (bIsAllEdgeMarked)
				{
					AttriIsCornerMarked[iBeginCorner] = true;
					AttriIsCornerMarked[iOppositeCorner] = true;
				}
			}

			AttriFacetDistance.unbind();
			pMesh->facets.attributes().delete_attribute_store("FacetDistance");
		}
	}

	void MeshCutAlgorithm::ComputeCuttingEdgePoints(In GEO::Mesh * const pMesh)
	{
		GEO::Mesh CopyMesh;
		CopyMesh.copy(*pMesh);

		HalfedgeMeshWrapper Wrapper(&CopyMesh);
		GEO::Attribute<bool> AttriIsCornerMarked(CopyMesh.facet_corners.attributes(), "IsCornerMarked");

		/*Update the m_NotMarkedCornersIdx after cutting branch*/
		m_NotMarkedCornersIdx.clear();
		for (GEO::index_t iCorner = 0; iCorner < pMesh->facet_corners.nb(); iCorner++)
		{
			if (!AttriIsCornerMarked[iCorner])
			{
				m_NotMarkedCornersIdx.push_back(iCorner);
			}
		}

		for (GEO::index_t i = 0; i < m_NotMarkedCornersIdx.size(); i++)
		{
			GEO::index_t iCorner = m_NotMarkedCornersIdx[i];
			if (!AttriIsCornerMarked[iCorner])
			{
				AttriIsCornerMarked[iCorner] = true;

				GEO::index_t iOppositeCorner = Wrapper.Opposite(iCorner);
				if (iOppositeCorner != GEO::NO_CORNER)
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

	void MeshCutAlgorithm::CutMeshByEdge(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper)
	{
		assert(pHalfedgeMeshWrapper != nullptr);
		assert(pHalfedgeMeshWrapper->pMesh != nullptr);

		GEO::Mesh * pMesh = pHalfedgeMeshWrapper->pMesh;

		GEO::Attribute<bool> AttriIsCornerMarked(pMesh->facet_corners.attributes(), "IsCornerMarked");

		/*Update the m_NotMarkedCornersIdx after cutting branch*/
		m_NotMarkedCornersIdx.clear();
		for (GEO::index_t iCorner = 0; iCorner < pMesh->facet_corners.nb(); iCorner++)
		{
			if (!AttriIsCornerMarked[iCorner])
			{
				m_NotMarkedCornersIdx.push_back(iCorner);
			}
		}

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
			GEO::index_t iCurrentVertex = iOldVertex;

			pMesh->facet_corners.set_adjacent_facet(iBeginCorner, GEO::NO_CORNER);

			double * pOldVertex = pMesh->vertices.point_ptr(iOldVertex);
			double * pNewVertex = pMesh->vertices.point_ptr(iNewVertex);

			pNewVertex[0] = pOldVertex[0];
			pNewVertex[1] = pOldVertex[1];
			pNewVertex[2] = pOldVertex[2];

			do
			{
				GEO::index_t iFacet = pHalfedgeMeshWrapper->Facet(iCorner);
				GEO::index_t iLocalVertex = pMesh->facets.find_vertex(iFacet, iOldVertex);
				pMesh->facets.set_vertex(iFacet, iLocalVertex, iNewVertex);

				iCorner = pHalfedgeMeshWrapper->NextAroundVertex(iCorner);

				/*This edge has been cut*/
				if (iCorner == GEO::NO_CORNER)
				{
					break;
				}

				iCurrentVertex = pMesh->facet_corners.vertex(iCorner);

			} while (iCorner != GEO::NO_CORNER && iCurrentVertex == iOldVertex && AttriIsCornerMarked[iCorner]);

			AttriIsCornerMarked[iBeginCorner] = true;
		}

		pMesh->vertices.remove_isolated();
		pMesh->facets.connect();
	}

	void MeshCutAlgorithm::ClearAttribute(In GEO::Mesh * const pMesh)
	{
		GEO::AttributesManager & FacetAttributeManager = pMesh->facets.attributes();
		GEO::AttributesManager & CornerAttributeManager = pMesh->facet_corners.attributes();

		GEO::vector<std::string> AttributeNames;

		FacetAttributeManager.list_attribute_names(AttributeNames);
		for (std::string & Name : AttributeNames)
		{
			if (Name != "region")
			{
				FacetAttributeManager.delete_attribute_store(Name);
			}
		}

		CornerAttributeManager.list_attribute_names(AttributeNames);
		for (std::string & Name : AttributeNames)
		{
			CornerAttributeManager.delete_attribute_store(Name);
		}
	}
}