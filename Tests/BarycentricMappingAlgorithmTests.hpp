#pragma once

#include <gtest\gtest.h>
#include <geogram\mesh\mesh_io.h>
#include <geogram\mesh\mesh_reorder.h>
#include <geogram\mesh\mesh_AABB.h>

#include <MeshGenerator.hpp>
#include <HalfedgeMeshWrapper.hpp>
#include <Algorithm\BarycentricMappingAlgorithm.hpp>

#include <iostream>

TEST(BarycentricMappingAlgorithm, Function)
{
	{
		std::string Root = "..\\Mesh\\";
		std::vector<std::string> BenchMarkModels =
		{
			"parameterization_basic_benchmark.obj"
		};

		std::vector<double> CorrectU =
		{
			(0.4538 + 1.0) * 0.5
		};

		for (GEO::index_t i = 0; i < BenchMarkModels.size(); i++)
		{
			GEO::Mesh Mesh;
			GEO::mesh_load(Root + BenchMarkModels[i], Mesh);
			GEO::mesh_reorder(Mesh, GEO::MESH_ORDER_MORTON);

			MeshAlgorithm::BarycentricMappingAlgorithm IAlgorithm;

			std::string CoefficientType = MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[0];
			std::string DomainShape = MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_DOMAIN_SHAPE[0];
			std::string BoundaryFixWeight = MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_BOUNDARY_FIX_WEIGHT[0];

			IAlgorithm.PutArg(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_KEY_COEFFICIENT_TYPE, CoefficientType);
			IAlgorithm.PutArg(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_KEY_DOMAIN_SHAPE, DomainShape);
			IAlgorithm.PutArg(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_KEY_BOUNDARY_FIX_WEIGHT, BoundaryFixWeight);
			IAlgorithm.Execute(&Mesh);

			GEO::Attribute<double> AttriVertexTexCoord(Mesh.vertices.attributes(), "tex_coord");

			ASSERT_NEAR(AttriVertexTexCoord[4 * 2 + 0], CorrectU[i], 1e-4);
			ASSERT_NEAR(AttriVertexTexCoord[4 * 2 + 1], 0.5, 1e-4);
		}
	}

	{
		std::string Root = "..\\Mesh\\";
		std::vector<std::string> BenchMarkModels =
		{
			"parameterization_basic_benchmark.obj"
		};

		std::vector<double> CorrectU =
		{
			(2.1138 + 1.0) * 0.5
		};

		for (GEO::index_t i = 0; i < BenchMarkModels.size(); i++)
		{
			GEO::Mesh Mesh;
			GEO::mesh_load(Root + BenchMarkModels[i], Mesh);
			GEO::mesh_reorder(Mesh, GEO::MESH_ORDER_MORTON);

			MeshAlgorithm::BarycentricMappingAlgorithm IAlgorithm;

			std::string CoefficientType = MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[1];
			std::string DomainShape = MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_DOMAIN_SHAPE[0];
			std::string BoundaryFixWeight = MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_BOUNDARY_FIX_WEIGHT[0];

			IAlgorithm.PutArg(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_KEY_COEFFICIENT_TYPE, CoefficientType);
			IAlgorithm.PutArg(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_KEY_DOMAIN_SHAPE, DomainShape);
			IAlgorithm.PutArg(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_KEY_BOUNDARY_FIX_WEIGHT, BoundaryFixWeight);
			IAlgorithm.Execute(&Mesh);

			GEO::Attribute<double> AttriVertexTexCoord(Mesh.vertices.attributes(), "tex_coord");

			ASSERT_NEAR(AttriVertexTexCoord[4 * 2 + 0], CorrectU[i], 1e-4);
			ASSERT_NEAR(AttriVertexTexCoord[4 * 2 + 1], 0.5, 1e-4);
		}
	}

	{
		std::string Root = "..\\Mesh\\";
		std::vector<std::string> BenchMarkModels =
		{
			"parameterization_basic_benchmark.obj"
		};

		std::vector<double> CorrectU =
		{
			(-35.1369 + 1.0) * 0.5
		};

		for (GEO::index_t i = 0; i < BenchMarkModels.size(); i++)
		{
			GEO::Mesh Mesh;
			GEO::mesh_load(Root + BenchMarkModels[i], Mesh);
			GEO::mesh_reorder(Mesh, GEO::MESH_ORDER_MORTON);

			MeshAlgorithm::BarycentricMappingAlgorithm IAlgorithm;

			std::string CoefficientType = MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[2];
			std::string DomainShape = MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_DOMAIN_SHAPE[0];
			std::string BoundaryFixWeight = MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_BOUNDARY_FIX_WEIGHT[0];

			IAlgorithm.PutArg(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_KEY_COEFFICIENT_TYPE, CoefficientType);
			IAlgorithm.PutArg(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_KEY_DOMAIN_SHAPE, DomainShape);
			IAlgorithm.PutArg(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_KEY_BOUNDARY_FIX_WEIGHT, BoundaryFixWeight);
			IAlgorithm.Execute(&Mesh);

			GEO::Attribute<double> AttriVertexTexCoord(Mesh.vertices.attributes(), "tex_coord");

			ASSERT_NEAR(AttriVertexTexCoord[4 * 2 + 0], CorrectU[i], 1e-4);
			ASSERT_NEAR(AttriVertexTexCoord[4 * 2 + 1], 0.5, 1e-4);
		}
	}
}