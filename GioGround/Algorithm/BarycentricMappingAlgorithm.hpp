#pragma once

#include "IMeshAlgorithm.hpp"

//#define USE_MKL
//#define USE_PARDISO
#define USE_GEO_NL

namespace MeshAlgorithm
{
	// Params:
	// CoefficientType -> "std::string" Required, ["Mean Value Coordinates", "Discrete Harmonic Coordinates", "Wachspress Coordinates"]
	// DomainShape -> "std::string" Required, ["Square", "Circle"]
	// BoundaryFixWeight -> "std::string" Required, ["Mean", "Length"]
	class BarycentricMappingAlgorithm : public IMeshAlgorithm
	{
	private:
		std::string m_CoefficientType = PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[0];
		std::string m_DomainShape = PARAMS_VALUE_SUPPORTED_DOMAIN_SHAPE[0];
		std::string m_BoundaryFixWeight = PARAMS_VALUE_SUPPORTED_BOUNDARY_FIX_WEIGHT[0];

		GEO::index_t m_nInteriorVertices = GEO::index_t(-1); // n
		GEO::index_t m_nBoundaryVertices = GEO::index_t(-1); // b
		GEO::vector<GEO::index_t> m_iInteriorVertices;
		GEO::vector<GEO::index_t> m_iBoundaryVertices;

		GEO::vector<double/*Dim = 3*/> m_InteriorVertices;
		GEO::vector<double/*Dim = 3*/> m_BoundaryVertices;

		GEO::Mesh m_ParameterizationMesh;

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
		// Step1
		void FindInteriorBoundary(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper);

		// Step2
		void FixBoundaryVertices(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper);
		void FixSquareBoundaryVertices(In GEO::vector<double> & Weight);
		void FixSphereBoundaryVertices(In GEO::vector<double> & Weight);

		// Step3
		void SoloveLinearEquation(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper);

		double Lambda_ij_BarycentricCoordinates(
			In HalfedgeMeshWrapper * pHalfedgeMeshWrapper,
			In GEO::index_t i, /*i = 0....n-1*/
			In GEO::index_t j
		) const;
		
		double w_ij_WachspressCoordinates(
			In HalfedgeMeshWrapper * pHalfedgeMeshWrapper,
			In GEO::index_t i, /*i = 0....n-1*/
			In GEO::index_t j,
			In GEO::index_t iAlpha_ij_Corner,
			In GEO::index_t iAlpha_ji_Corner
		) const;
		double w_ij_DiscreteHarmonicCoordinates(
			In HalfedgeMeshWrapper * pHalfedgeMeshWrapper,
			In GEO::index_t i, /*i = 0....n-1*/
			In GEO::index_t j,
			In GEO::index_t iAlpha_ij_Corner,
			In GEO::index_t iAlpha_ji_Corner
		) const;
		double w_ij_MeanValueCoordinates(
			In HalfedgeMeshWrapper * pHalfedgeMeshWrapper,
			In GEO::index_t i, /*i = 0....n-1*/
			In GEO::index_t j,
			In GEO::index_t iAlpha_ij_Corner,
			In GEO::index_t iAlpha_ji_Corner
		) const;

		// Step 4
		void GenerateOutput(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper);
	};
}