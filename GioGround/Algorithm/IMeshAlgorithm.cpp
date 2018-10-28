#include "IMeshAlgorithm.hpp"

namespace MeshAlgorithm
{
	IMeshAlgorithm::~IMeshAlgorithm()
	{

	}

	bool IMeshAlgorithm::Execute(InOut GEO::Mesh * pMesh)
	{
		if (pMesh == nullptr)
		{
			return false;
		}

		return CheckAndGetArgs(pMesh);
	}

	bool IMeshAlgorithm::ExecuteOut(In GEO::Mesh * pMesh, Out GEO::Mesh * pOutMesh)
	{
		if (pMesh == nullptr || pOutMesh == nullptr)
		{
			return false;
		}

		pOutMesh->copy(*pMesh);
		return Execute(pOutMesh);
	}

	bool IMeshAlgorithm::Visualize(In GEO::Mesh * const pMesh) const
	{
		return true;
	}

	GEO::index_t IMeshAlgorithm::GetBoundaryNumber(In GEO::Mesh * const pMesh) const
	{
		assert(pMesh != nullptr);

		HalfedgeMeshWrapper Wrapper(pMesh);
		GEO::Attribute<GEO::index_t> AttriCornerBoundaryIndex;

		AttriCornerBoundaryIndex.bind(pMesh->facet_corners.attributes(), "BoundaryIndex");
		AttriCornerBoundaryIndex.fill(GEO::index_t(-1));

		GEO::vector<GEO::index_t> BoundaryCorners;

		for (GEO::index_t iCorner = 0; iCorner < pMesh->facet_corners.nb(); iCorner++)
		{
			if (pMesh->facet_corners.adjacent_facet(iCorner) == GEO::NO_FACET)
			{
				BoundaryCorners.push_back(iCorner);
			}
		}

		GEO::index_t iBoundary = 0;
		for (GEO::index_t i = 0; i < BoundaryCorners.size(); i++)
		{
			GEO::index_t iBeginCorner = BoundaryCorners[i];
			GEO::index_t iCorner = iBeginCorner;

			if (AttriCornerBoundaryIndex[iBeginCorner] != GEO::index_t(-1))
			{
				continue;
			}

			do
			{
				AttriCornerBoundaryIndex[iCorner] = iBoundary;
				iCorner = Wrapper.Next(iCorner);
				do
				{
					iCorner = Wrapper.Corner2Corner[iCorner];
				} while (pMesh->facet_corners.adjacent_facet(iCorner) != GEO::NO_FACET);
			} while (iCorner != iBeginCorner);
			
			iBoundary++;
		}

		AttriCornerBoundaryIndex.unbind();
		pMesh->facet_corners.attributes().delete_attribute_store("BoundaryIndex");

		return iBoundary;
	}

	bool IMeshAlgorithm::IsVertexAdjacent(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper, In GEO::index_t i, In GEO::index_t j) const
	{
		assert(i != j);
		assert(i < pHalfedgeMeshWrapper->pMesh->vertices.nb() && j < pHalfedgeMeshWrapper->pMesh->vertices.nb());

		GEO::index_t iIBeginCorner = pHalfedgeMeshWrapper->Vertex2Corner[i];
		GEO::index_t iICorner = iIBeginCorner;

		do
		{
			if (pHalfedgeMeshWrapper->pMesh->facet_corners.vertex(pHalfedgeMeshWrapper->Prev(iICorner)) == j)
			{
				return true;
			}

			if (pHalfedgeMeshWrapper->pMesh->facet_corners.adjacent_facet(iICorner) == GEO::NO_FACET)
			{
				if (pHalfedgeMeshWrapper->Dest(iICorner) == j)
				{
					return true;
				}
			}

			iICorner = pHalfedgeMeshWrapper->Corner2Corner[iICorner];
		} while (iICorner != iIBeginCorner);
		return false;
	}

	GEO::vector<GEO::index_t> IMeshAlgorithm::GetAdjacentVertices(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper, In GEO::index_t i) const
	{
		assert(i < pHalfedgeMeshWrapper->pMesh->vertices.nb());

		GEO::vector<GEO::index_t> Adj;

		GEO::index_t iBeginCorner = pHalfedgeMeshWrapper->Vertex2Corner[i];
		GEO::index_t iCorner = iBeginCorner;

		do
		{
			Adj.push_back(pHalfedgeMeshWrapper->pMesh->facet_corners.vertex(pHalfedgeMeshWrapper->Prev(iCorner)));
			if (pHalfedgeMeshWrapper->pMesh->facet_corners.adjacent_facet(iCorner) == GEO::NO_FACET)
			{
				Adj.push_back(pHalfedgeMeshWrapper->Dest(iCorner));
			}
			iCorner = pHalfedgeMeshWrapper->Corner2Corner[iCorner];
		} while (iCorner != iBeginCorner);

		return Adj;
	}

	bool IMeshAlgorithm::IsBoundaryVertex(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper, In GEO::index_t i) const
	{
		assert(i < pHalfedgeMeshWrapper->pMesh->vertices.nb());

		GEO::index_t iBeginCorner = pHalfedgeMeshWrapper->Vertex2Corner[i];
		GEO::index_t iCorner = iBeginCorner;

		do
		{
			if (pHalfedgeMeshWrapper->pMesh->facet_corners.adjacent_facet(iCorner) == GEO::NO_FACET)
			{
				return true;
			}
			iCorner = pHalfedgeMeshWrapper->Corner2Corner[iCorner];
		} while (iCorner != iBeginCorner);

		return false;
	}

	bool IMeshAlgorithm::Reset()
	{
		return true;
	}

	bool IMeshAlgorithm::CheckAndGetArgs(In GEO::Mesh * pMesh)
	{
		pMesh->assert_is_valid();
		return true;
	}
}