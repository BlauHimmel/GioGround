#pragma once

#include <gtest\gtest.h>
#include <geogram\mesh\mesh_io.h>

#include <Algorithm\IMeshAlgorithm.hpp>

TEST(IMeshAlgorithm, Function_Execute)
{
	GEO::Mesh Mesh, OutMesh;
	MeshAlgorithm::IMeshAlgorithm IAlgorithm;
	ASSERT_TRUE(IAlgorithm.Execute(&Mesh));
	ASSERT_TRUE(IAlgorithm.ExecuteOut(&Mesh, &OutMesh));
}

TEST(IMeshAlgorithm, Function_PutArg_SetArg)
{
	MeshAlgorithm::IMeshAlgorithm IAlgorithm;

	int32_t Int32Arg = 123456;
	uint32_t UInt32Arg = 123456;
	float FloatArg = 123456.7f;
	double DoubleArg = 123456.7;
	GEO::vector<int32_t> VectorInt32Arg(4);
	VectorInt32Arg[0] = 0;
	VectorInt32Arg[1] = 1;
	VectorInt32Arg[2] = 2;
	VectorInt32Arg[3] = 3;
	GEO::vector<uint32_t> VectorUInt32Arg(4);
	VectorUInt32Arg[0] = 0u;
	VectorUInt32Arg[1] = 1u;
	VectorUInt32Arg[2] = 2u;
	VectorUInt32Arg[3] = 3u;

	ASSERT_TRUE(IAlgorithm.PutArg("Int32Arg", Int32Arg));
	ASSERT_FALSE(IAlgorithm.PutArg("Int32Arg", UInt32Arg));
	ASSERT_TRUE(IAlgorithm.PutArg("UInt32Arg", UInt32Arg));

	ASSERT_TRUE(IAlgorithm.PutArg("FloatArg", FloatArg));
	ASSERT_FALSE(IAlgorithm.PutArg("FloatArg", DoubleArg));
	ASSERT_TRUE(IAlgorithm.PutArg("DoubleArg", DoubleArg));

	ASSERT_TRUE(IAlgorithm.PutArg("VectorInt32Arg", VectorInt32Arg));
	ASSERT_FALSE(IAlgorithm.PutArg("VectorInt32Arg", VectorUInt32Arg));
	ASSERT_TRUE(IAlgorithm.PutArg("VectorUInt32Arg", VectorUInt32Arg));

	int32_t Int32ArgOut = 0;
	uint32_t UInt32ArgOut = 0u;
	float FloatArgOut = 0.0f;
	double DoubleArgOut = 0.0;
	GEO::vector<int32_t> VectorInt32ArgOut;
	GEO::vector<uint32_t> VectorUInt32ArgOut;

	ASSERT_TRUE(IAlgorithm.GetArg("Int32Arg", &UInt32ArgOut));
	ASSERT_EQ(UInt32ArgOut, UInt32Arg);
	ASSERT_TRUE(IAlgorithm.GetArg("UInt32Arg", &UInt32ArgOut));
	ASSERT_EQ(UInt32ArgOut, UInt32Arg);

	ASSERT_TRUE(IAlgorithm.GetArg("FloatArg", &DoubleArgOut));
	ASSERT_EQ(DoubleArgOut, DoubleArg);
	ASSERT_TRUE(IAlgorithm.GetArg("DoubleArg", &DoubleArgOut));
	ASSERT_EQ(DoubleArgOut, DoubleArg);

	ASSERT_TRUE(IAlgorithm.GetArg("VectorInt32Arg", &VectorUInt32ArgOut));
	for (GEO::index_t i = 0; i < VectorUInt32ArgOut.size(); i++)
	{
		ASSERT_EQ(VectorUInt32Arg[i], VectorUInt32ArgOut[i]);
	}
	ASSERT_TRUE(IAlgorithm.GetArg("VectorUInt32Arg", &VectorUInt32ArgOut));
	for (GEO::index_t i = 0; i < VectorUInt32ArgOut.size(); i++)
	{
		ASSERT_EQ(VectorUInt32Arg[i], VectorUInt32ArgOut[i]);
	}

	Int32ArgOut = 0;
	UInt32ArgOut = 0u;
	FloatArgOut = 0.0f;
	DoubleArgOut = 0.0;
	VectorInt32ArgOut.clear();
	VectorUInt32ArgOut.clear();

	ASSERT_FALSE(IAlgorithm.GetArg("Int32Arg", &Int32ArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("Int32Arg", &FloatArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("Int32Arg", &DoubleArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("Int32Arg", &VectorInt32ArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("Int32Arg", &VectorUInt32ArgOut));

	Int32ArgOut = 0;
	UInt32ArgOut = 0u;
	FloatArgOut = 0.0f;
	DoubleArgOut = 0.0;
	VectorInt32ArgOut.clear();
	VectorUInt32ArgOut.clear();

	ASSERT_FALSE(IAlgorithm.GetArg("UInt32Arg", &Int32ArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("UInt32Arg", &FloatArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("UInt32Arg", &DoubleArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("UInt32Arg", &VectorInt32ArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("UInt32Arg", &VectorUInt32ArgOut));

	Int32ArgOut = 0;
	UInt32ArgOut = 0u;
	FloatArgOut = 0.0f;
	DoubleArgOut = 0.0;
	VectorInt32ArgOut.clear();
	VectorUInt32ArgOut.clear();

	ASSERT_FALSE(IAlgorithm.GetArg("FloatArg", &Int32ArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("FloatArg", &UInt32ArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("FloatArg", &FloatArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("FloatArg", &VectorInt32ArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("FloatArg", &VectorUInt32ArgOut));

	Int32ArgOut = 0;
	UInt32ArgOut = 0u;
	FloatArgOut = 0.0f;
	DoubleArgOut = 0.0;
	VectorInt32ArgOut.clear();
	VectorUInt32ArgOut.clear();

	ASSERT_FALSE(IAlgorithm.GetArg("DoubleArg", &Int32ArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("DoubleArg", &UInt32ArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("DoubleArg", &FloatArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("DoubleArg", &VectorInt32ArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("DoubleArg", &VectorUInt32ArgOut));

	Int32ArgOut = 0;
	UInt32ArgOut = 0u;
	FloatArgOut = 0.0f;
	DoubleArgOut = 0.0;
	VectorInt32ArgOut.clear();
	VectorUInt32ArgOut.clear();

	ASSERT_FALSE(IAlgorithm.GetArg("VectorInt32", &Int32ArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("VectorInt32", &UInt32ArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("VectorInt32", &FloatArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("VectorInt32", &DoubleArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("VectorInt32", &VectorInt32ArgOut));

	Int32ArgOut = 0;
	UInt32ArgOut = 0u;
	FloatArgOut = 0.0f;
	DoubleArgOut = 0.0;
	VectorInt32ArgOut.clear();
	VectorUInt32ArgOut.clear();

	ASSERT_FALSE(IAlgorithm.GetArg("VectorUInt32", &Int32ArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("VectorUInt32", &UInt32ArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("VectorUInt32", &FloatArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("VectorUInt32", &DoubleArgOut));
	ASSERT_FALSE(IAlgorithm.GetArg("VectorUInt32", &VectorInt32ArgOut));
}

TEST(IMeshAlgorithm, Function_GetBoundaryNumber)
{
	std::string Root = "..\\Mesh\\";
	std::vector<std::string> BenchMarkModels =
	{
		"plane_hole1.obj",
		"plane_hole2.obj",
		"donut_genus1.obj",
		"simple_donut_genus1.obj",
		"sphere_hole1.obj",
		"sphere_hole2.obj",
		"kitten.obj",
		"kitten_holes.obj",
		"king.obj"
	};

	std::vector<GEO::index_t> BoundaryNumber =
	{
		2,
		3,
		0,
		0,
		1,
		2,
		0,
		3,
		0
	};

	for (GEO::index_t i = 0; i < BenchMarkModels.size(); i++)
	{
		GEO::Mesh Mesh;
		GEO::mesh_load(Root + BenchMarkModels[i], Mesh);
		std::unique_ptr<MeshAlgorithm::IMeshAlgorithm> IMeshAlgorithm(new MeshAlgorithm::IMeshAlgorithm());
		ASSERT_EQ(IMeshAlgorithm->GetBoundaryNumber(&Mesh), BoundaryNumber[i]);;
	}
}
