#pragma once
#pragma warning(disable : 4199)

#include <geogram_gfx\basic\GLSL.h>
#include <geogram_gfx\basic\GL.h>
#include <geogram_gfx\GLUP\GLUP.h>
#include <geogram_gfx\glup_viewer\glup_viewer.h>

#include <unordered_map>
#include <string>
#include <any>

#include "..\Macro.hpp"
#include "..\HalfedgeMeshWrapper.hpp"

namespace MeshAlgorithm
{
	
	// VERY IMPORTANT THING : The order of element of the input mesh must be guarantee, otherwise you will obtain an WRONG ANSWER.
	// That is to say, GEO::mesh_reorder() must be called first. In the GiogroundApp, we dont need to care about this 
	// issue since MeshFacetsAABB structure is build automatically, which implicitly call this function. But in other circumstance,
	// GEO::mesh_reorder() must be called first!!!!
	class IMeshAlgorithm
	{
	private:
		std::unordered_map<std::string, std::any> m_Args;

	public:
		virtual ~IMeshAlgorithm();
		virtual bool Execute(InOut GEO::Mesh * pMesh);
		bool ExecuteOut(In GEO::Mesh * pMesh, Out GEO::Mesh * pOutMesh);

		virtual bool Visualize(In GEO::Mesh * const pMesh) const;

	public:
		GEO::index_t GetBoundaryNumber(In GEO::Mesh * const pMesh) const;
		bool IsVertexAdjacent(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper, In GEO::index_t i, In GEO::index_t j) const;
		GEO::vector<GEO::index_t> GetAdjacentVertices(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper, In GEO::index_t i) const;
		bool IsBoundaryVertex(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper, In GEO::index_t i) const;
		void GetBBox(In GEO::Mesh * const pMesh, Out double * XYZMin, Out double * XYZMax, bool bAnimate) const;
				/**		* The SOD is the simplest classification operator. It assigns a weight		* to every edge in the mesh proportional to the dihedral angle defined by		* the normals of its two adjacent triagles.		* Ref : Multiresolution Feature Extraction from Unstructured Meshes		* - Andreas Hubeli, Markus Gross		*/		double SecondOrderDifference(In HalfedgeMeshWrapper * pHalfedgeMeshWrapper, In GEO::index_t iCorner) const;

	protected:
		virtual bool Reset();
		virtual bool CheckAndGetArgs(In GEO::Mesh * pMesh);

	public:
		template <class T>
		bool PutArg(In const std::string & Key, In T Arg)
		{
			if (m_Args.find(Key) != m_Args.end())
			{
				m_Args[Key] = std::make_any<T>(Arg);
				return false;
			}
			m_Args.insert(std::make_pair(Key, std::make_any<T>(Arg)));
			return true;
		}

		template <class T>
		bool GetArg(In const std::string & Key, Out T * Arg) const
		{
			if (m_Args.find(Key) == m_Args.end())
			{
				return false;
			}

			std::any RawArg = m_Args.at(Key);

			if (RawArg.type() != typeid(T))
			{
				return false;
			}

			*Arg = std::any_cast<T>(RawArg);
			return true;
		}
	};
}