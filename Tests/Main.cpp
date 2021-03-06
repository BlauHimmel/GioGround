#include <gtest\gtest.h>

#include "HalfedgeMeshWrapperTests.hpp"
#include "IMeshAlgorithmTests.hpp"
#include "MeshCutAlgorithmTests.hpp"
#include "BarycentricMappingAlgorithmTests.hpp"

int main(int argc, char ** argv)
{
	GEO::initialize();
	GEO::set_assert_mode(GEO::AssertMode::ASSERT_BREAKPOINT);
	GEO::Logger::instance()->set_quiet(true);
	testing::InitGoogleTest(&argc, argv);
	RUN_ALL_TESTS();
	system("PAUSE");
	return 0;
}
