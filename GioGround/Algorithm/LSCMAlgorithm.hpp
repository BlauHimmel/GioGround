#pragma once

#include "IMeshAlgorithm.hpp"

namespace MeshAlgorithm
{
	class LSCMAlgorithm : public IMeshAlgorithm
	{
	public:
		virtual bool Execute(InOut GEO::Mesh * pMesh) override;
		virtual bool Visualize(In GEO::Mesh * const pMesh) const override;

	protected:
		virtual bool CheckAndGetArgs(In GEO::Mesh * const pMesh) override;
		virtual bool Reset() override;

	private:

		// Computes the coordinates of the vertices of a triangle
	    // in a local 2D orthonormal basis of the triangle's plane.
		void ProjectTriangle(
			In const GEO::vec3 & P0,
			In const GEO::vec3 & P1,
			In const GEO::vec3 & P2,
			Out GEO::vec2 & Z0,
			Out GEO::vec2 & Z1,
			Out GEO::vec2 & Z2
		) const;

		void SolveLeastSquareEquation(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper);
	};
}