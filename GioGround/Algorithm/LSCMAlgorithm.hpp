#pragma once

#include "IMeshAlgorithm.hpp"

namespace MeshAlgorithm
{
	class LSCMAlgorithm : public IMeshAlgorithm
	{
	private:
		GEO::vector<double/*Dim = 2 UV*/> m_Vertices;

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

		// Get the bounding box of the mesh
		void GetBoundingBox(
			In GEO::Mesh * const pMesh,
			Out double & XMin,
			Out double & YMin,
			Out double & ZMin,
			Out double & XMax,
			Out double & YMax,
			Out double & ZMax
		) const;

		// Get the dimension whose area is biggest(This dimension is perpendicular 
		// to the shortest axis)
		void GetLargestDimension(
			In double XMin,
			In double YMin,
			In double ZMin,
			In double XMax,
			In double YMax,
			In double ZMax,
			Out GEO::vec3 & V1,
			Out GEO::vec3 & V2
		) const;

		// Compute the local coordinate of vertices and fixed two point
		void FixedPoint(
			In GEO::Mesh * const pMesh,
			In GEO::vec3 V1,
			In GEO::vec3 V2,
			Out GEO::index_t & iFixedPointMin,
			Out GEO::index_t & iFixedPointMax
		) const;

		// Step1 : Build linear system and solve it
		void SolveLeastSquareEquation(In GEO::Mesh * const pMesh);

		// Step2 : Normalize the uv so that they fit within the unit square
		void NormalizeUV(In GEO::Mesh * const pMesh);

		// Step3
		void GenerateOutput(In GEO::Mesh * const pMesh);

		// Step4
		void ClearAttribute(In GEO::Mesh * const pMesh);
	};
}