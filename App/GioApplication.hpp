#pragma once

#include <geogram_gfx\basic\GLSL.h>
#include <geogram_gfx\basic\GL.h>
#include <geogram_gfx\GLUP\GLUP.h>
#include <geogram_gfx\glup_viewer\glup_viewer.h>
#include <geogram\basic\command_line.h>
#include <geogram\basic\command_line_args.h>
#include <geogram\basic\logger.h>
#include <geogram\mesh\mesh.h>
#include <geogram\mesh\mesh_AABB.h>
#include <geogram\mesh\mesh_io.h>
#include <geogram\basic\assert.h>
#include <geogram\basic\stopwatch.h>

#include <memory>
#include <limits>
#include <vector>

#include <MeshGenerator.hpp>
#include <Algorithm\MeshCutAlgorithm.hpp>
#include <Algorithm\BarycentricMappingAlgorithm.hpp>
#include <Algorithm\LSCMAlgorithm.hpp>

extern class GioApplication * g_pApp;

class GioApplication : public GEO::SimpleMeshApplication
{
	friend GLboolean MouseCallbackFunc(float X, float Y, int Button, enum GlupViewerEvent Event);

	struct Algorithm
	{
		std::unique_ptr<MeshAlgorithm::IMeshAlgorithm> MeshAlgorithm = nullptr;
		bool bVisualizeAlgorithm = false;
		bool bRunAlgorithm = false;
		bool bShowDialog = false;
	};

protected:
	float m_BorderWidth = 0.1f;
	GEO::vec4f m_BorderColor = GEO::vec4f(1.0f, 0.0f, 0.0f, 1.0f);

	std::unique_ptr<GEO::MeshFacetsAABB> m_MashFacetsAABB;
	GEO::index_t m_iSelectedFacet = GEO::NO_FACET;
	bool m_bSelectingFacet = false;
	GEO::index_t m_iSelectedVertex = GEO::NO_VERTEX;
	bool m_bSelectingVertex = false;
	bool m_bAniamted = false;

	std::vector<Algorithm> m_Algorithms;
	size_t m_iCurrentAlgorithm = size_t(-1);

	const size_t ALGORITHM_NUMBER = 3;
	const size_t CUT_MESH_ALGORITHM_INDEX = 0;
	const size_t BARYCENTRIC_MAPPING_ALGORITHM_INDEX = 1;
	const size_t LSCM_ALGORITHM_INDEX = 2;

public:
	GioApplication(int argc, char ** argv);

	virtual void draw_scene() override;
	virtual void init_graphics() override;
	virtual void draw_gui() override;
	virtual void draw_object_properties() override;
	virtual bool load(const std::string & Filename) override;
	virtual bool save(const std::string & Filename) override;

protected:
	virtual int PANE_WIDTH() const override;
	virtual void draw_application_menus() override;

private:
	void DrawAlgorithmDialog();
	void CloseAlgorithmDialog();

	void DrawMeshCutAlgorithmDialog();
	void CloseMeshCutAlgorithmDialog();

	void DrawBarycentricMappingAlgorithmDialog();
	void CloseBarycentricMappingAlgorithmDialog();

	void DrawLSCMAlgorithmDialog();
	void CloseLSCMAlgorithmDialog();
};

GLboolean MouseCallbackFunc(float X, float Y, int Button, enum GlupViewerEvent Event);