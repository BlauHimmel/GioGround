#include <gtest\gtest.h>

#include "HalfedgeMeshWrapperTests.hpp"
#include "IMeshAlgorithmTests.hpp"
#include "MeshCutAlgorithmTests.hpp"

int main(int argc, char ** argv)
{
	testing::InitGoogleTest(&argc, argv);
	RUN_ALL_TESTS();
	system("PAUSE");
	return 0;
}
