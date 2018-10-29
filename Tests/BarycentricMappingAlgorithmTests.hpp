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
	GEO::Mesh Mesh;
	MeshGenerator::MeshGenPyramid(&Mesh);
	MeshAlgorithm::BarycentricMappingAlgorithm IAlgorithm;

	std::string CoefficientType = MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[1];
	std::string DomainShape = MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_DOMAIN_SHAPE[0];
	std::string BoundaryFixWeight = MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_BOUNDARY_FIX_WEIGHT[0];

	IAlgorithm.PutArg(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_KEY_COEFFICIENT_TYPE, CoefficientType);
	IAlgorithm.PutArg(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_KEY_DOMAIN_SHAPE, DomainShape);
	IAlgorithm.PutArg(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_KEY_BOUNDARY_FIX_WEIGHT, BoundaryFixWeight);

	GEO::vec3 v0 = Mesh.vertices.point(0);
	GEO::vec3 v1 = Mesh.vertices.point(1);
	GEO::vec3 v2 = Mesh.vertices.point(2);
	GEO::vec3 v3 = Mesh.vertices.point(3);
	GEO::vec3 v4 = Mesh.vertices.point(4);

	/*TODO : Finished it*/
}