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