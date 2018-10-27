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
		GEO::Mesh CopyMesh;
		CopyMesh.copy(*pMesh);

		HalfedgeMeshWrapper Wrapper(&CopyMesh);
		GEO::Attribute<GEO::index_t> AttriCornerBoundaryIndex;

		AttriCornerBoundaryIndex.bind(CopyMesh.facet_corners.attributes(), "BoundaryIndex");
		AttriCornerBoundaryIndex.fill(GEO::index_t(-1));

		GEO::vector<GEO::index_t> BoundaryCorners;

		for (GEO::index_t iCorner = 0; iCorner < CopyMesh.facet_corners.nb(); iCorner++)
		{
			if (CopyMesh.facet_corners.adjacent_facet(iCorner) == GEO::NO_FACET)
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
				} while (CopyMesh.facet_corners.adjacent_facet(iCorner) != GEO::NO_FACET);
			} while (iCorner != iBeginCorner);
			
			iBoundary++;
		}

		return iBoundary;
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