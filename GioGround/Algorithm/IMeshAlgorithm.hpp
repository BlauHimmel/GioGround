#pragma once

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
	class IMeshAlgorithm
	{
	private:
		std::unordered_map<std::string, std::any> m_Args;

	public:
		virtual ~IMeshAlgorithm();
		virtual bool Execute(InOut GEO::Mesh * pMesh);
		bool ExecuteOut(In GEO::Mesh * pMesh, Out GEO::Mesh * pOutMesh);

		virtual bool Visualize(In GEO::Mesh * const pMesh) const;

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