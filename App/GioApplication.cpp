#include "GioApplication.hpp"

GioApplication::GioApplication(int argc, char ** argv) : GEO::SimpleMeshApplication(argc, argv, "")
{

}

void GioApplication::draw_scene()
{
	if (mesh_gfx_.mesh() == nullptr)
	{
		return;
	}

	if (glup_viewer_is_enabled(GLUP_VIEWER_IDLE_REDRAW))
	{
		mesh_gfx_.set_animate(true);
		anim_time_ = float(sin(double(anim_speed_) * GEO::SystemStopwatch::now()));
		anim_time_ = 0.5f * (anim_time_ + 1.0f);
	}

	mesh_gfx_.set_lighting(lighting_);
	mesh_gfx_.set_time(double(anim_time_));

	if (show_attributes_)
	{
		mesh_gfx_.set_scalar_attribute(
			attribute_subelements_, attribute_name_,
			double(attribute_min_), double(attribute_max_),
			current_colormap_texture_, 1
		);
	}
	else 
	{
		mesh_gfx_.unset_scalar_attribute();
	}

	if (show_vertices_) 
	{
		mesh_gfx_.set_points_color(vertices_color_.x, vertices_color_.y, vertices_color_.z);
		mesh_gfx_.set_points_size(vertices_size_);
		mesh_gfx_.draw_vertices();
	}

	if (show_vertices_selection_)
	{
		mesh_gfx_.set_points_color(1.0, 0.0, 0.0);
		mesh_gfx_.set_points_size(2.0f * vertices_size_);
		mesh_gfx_.set_vertices_selection("selection");
		mesh_gfx_.draw_vertices();
		mesh_gfx_.set_vertices_selection("");
	}

	mesh_gfx_.set_mesh_color(0.0, 0.0, 0.0);

	mesh_gfx_.set_surface_color(surface_color_.x, surface_color_.y, surface_color_.z);
	if (show_surface_sides_) 
	{
		mesh_gfx_.set_backface_surface_color(surface_color_2_.x, surface_color_2_.y, surface_color_2_.z);
	}

	mesh_gfx_.set_show_mesh(show_mesh_);
	mesh_gfx_.set_mesh_color(mesh_color_.x, mesh_color_.y, mesh_color_.z);
	mesh_gfx_.set_mesh_width(GEO::index_t(mesh_width_ * 10.0f));

	if (show_surface_)
	{
		float specular_backup = glupGetSpecular();
		glupSetSpecular(0.4f);
		mesh_gfx_.draw_surface(); 
		glupSetSpecular(specular_backup);
	}

	if (show_mesh_)
	{
		mesh_gfx_.draw_edges();
	}

	mesh_gfx_.set_mesh_color(m_BorderColor.x, m_BorderColor.y, m_BorderColor.z);
	mesh_gfx_.set_mesh_border_width(GEO::index_t(m_BorderWidth * 10.0f));
	if (show_surface_borders_)
	{
		mesh_gfx_.draw_surface_borders(); 
	}

	if (show_volume_)
	{
		if (glupIsEnabled(GLUP_CLIPPING) && glupGetClipMode() == GLUP_CLIP_SLICE_CELLS)
		{
			mesh_gfx_.set_lighting(false);
		}

		mesh_gfx_.set_shrink(double(cells_shrink_));
		mesh_gfx_.set_draw_cells(GEO::MESH_HEX, show_hexes_);
		mesh_gfx_.set_draw_cells(GEO::MESH_CONNECTOR, show_connectors_);

		if (show_colored_cells_)
		{
			mesh_gfx_.set_cells_colors_by_type();
		}
		else
		{
			mesh_gfx_.set_cells_color(volume_color_.x, volume_color_.y, volume_color_.z);
		}
		mesh_gfx_.draw_volume();
		mesh_gfx_.set_lighting(lighting_);
	}

	if (m_bSelectingFacet)
	{
		if (!m_bDisplayAllSelectedFacets && m_iSelectedFacet != GEO::NO_FACET)
		{
			GEO::index_t iV1 = mesh_.facets.vertex(m_iSelectedFacet, 0);
			GEO::index_t iV2 = mesh_.facets.vertex(m_iSelectedFacet, 1);
			GEO::index_t iV3 = mesh_.facets.vertex(m_iSelectedFacet, 2);

			GEO::index_t EdgeWidth = 3;
			GEO::vec4f TriangleColor(0.0f, 1.0f, 0.0f, 1.0f);

			glupSetMeshWidth(EdgeWidth);
			glupSetColor4fv(GLUP_FRONT_AND_BACK_COLOR, TriangleColor.data());

			double * pV1 = mesh_.vertices.point_ptr(iV1);
			double * pV2 = mesh_.vertices.point_ptr(iV2);
			double * pV3 = mesh_.vertices.point_ptr(iV3);

			glupBegin(GLUP_LINES);
			glupVertex3d(pV1[0], pV1[1], pV1[2]);
			glupVertex3d(pV2[0], pV2[1], pV2[2]);
			glupVertex3d(pV2[0], pV2[1], pV2[2]);
			glupVertex3d(pV3[0], pV3[1], pV3[2]);
			glupVertex3d(pV3[0], pV3[1], pV3[2]);
			glupVertex3d(pV1[0], pV1[1], pV1[2]);
			glupEnd();
		}
		else
		{
			for (auto Iter = m_iSelectedFacets.begin(); Iter != m_iSelectedFacets.end(); ++Iter)
			{
				GEO::index_t iV1 = mesh_.facets.vertex(*Iter, 0);
				GEO::index_t iV2 = mesh_.facets.vertex(*Iter, 1);
				GEO::index_t iV3 = mesh_.facets.vertex(*Iter, 2);

				GEO::index_t EdgeWidth = 3;
				GEO::vec4f TriangleColor(0.0f, 1.0f, 0.0f, 1.0f);

				glupSetMeshWidth(EdgeWidth);
				glupSetColor4fv(GLUP_FRONT_AND_BACK_COLOR, TriangleColor.data());

				double * pV1 = mesh_.vertices.point_ptr(iV1);
				double * pV2 = mesh_.vertices.point_ptr(iV2);
				double * pV3 = mesh_.vertices.point_ptr(iV3);

				glupBegin(GLUP_LINES);
				glupVertex3d(pV1[0], pV1[1], pV1[2]);
				glupVertex3d(pV2[0], pV2[1], pV2[2]);
				glupVertex3d(pV2[0], pV2[1], pV2[2]);
				glupVertex3d(pV3[0], pV3[1], pV3[2]);
				glupVertex3d(pV3[0], pV3[1], pV3[2]);
				glupVertex3d(pV1[0], pV1[1], pV1[2]);
				glupEnd();
			}
		}
	}

	if (m_bSelectingVertex)
	{
		if (!m_bDisplayAllSelectedVertices && m_iSelectedVertex != GEO::NO_VERTEX)
		{
			double * pVertex = mesh_.vertices.point_ptr(m_iSelectedVertex);

			GLUPfloat PointSize = 8.0;
			GEO::vec4f PointColor(0.0f, 1.0f, 0.0f, 1.0f);

			glupSetPointSize(PointSize);
			glupSetColor4fv(GLUP_FRONT_AND_BACK_COLOR, PointColor.data());

			glupBegin(GLUP_POINTS);
			glupVertex3d(pVertex[0], pVertex[1], pVertex[2]);
			glupEnd();
		}
		else
		{
			for (auto Iter = m_iSelectedVertices.begin(); Iter != m_iSelectedVertices.end(); ++Iter)
			{
				double * pVertex = mesh_.vertices.point_ptr(*Iter);

				GLUPfloat PointSize = 8.0;
				GEO::vec4f PointColor(0.0f, 1.0f, 0.0f, 1.0f);

				glupSetPointSize(PointSize);
				glupSetColor4fv(GLUP_FRONT_AND_BACK_COLOR, PointColor.data());

				glupBegin(GLUP_POINTS);
				glupVertex3d(pVertex[0], pVertex[1], pVertex[2]);
				glupEnd();
			}
		}
	}
	
	if (m_iCurrentAlgorithm != size_t(-1) && m_Algorithms[m_iCurrentAlgorithm].bVisualizeAlgorithm && m_Algorithms[m_iCurrentAlgorithm].MeshAlgorithm != nullptr)
	{
		m_Algorithms[m_iCurrentAlgorithm].MeshAlgorithm->Visualize(&mesh_);
	}
}

void GioApplication::init_graphics()
{
	GEO::SimpleMeshApplication::init_graphics();
	GEO::set_assert_mode(GEO::ASSERT_BREAKPOINT);
	glup_viewer_set_mouse_func(MouseCallbackFunc);
	retina_mode_ = false;
	scaling_ = 1.0f;

	m_Algorithms.resize(ALGORITHM_NUMBER);
	m_Algorithms[CUT_MESH_ALGORITHM_INDEX].MeshAlgorithm = std::make_unique<MeshAlgorithm::MeshCutAlgorithm>();
	m_Algorithms[BARYCENTRIC_MAPPING_ALGORITHM_INDEX].MeshAlgorithm = std::make_unique<MeshAlgorithm::BarycentricMappingAlgorithm>();
	m_Algorithms[LSCM_ALGORITHM_INDEX].MeshAlgorithm = std::make_unique<MeshAlgorithm::LSCMAlgorithm>();
}

void GioApplication::draw_gui()
{
	GEO::SimpleMeshApplication::draw_gui();
}

void GioApplication::draw_object_properties()
{
	ImGui::Checkbox("Attributes", &show_attributes_);
	if (show_attributes_) 
	{
		if (attribute_min_ == 0.0f && attribute_max_ == 0.0f)
		{
			autorange();
		}

		if (ImGui::Button((attribute_ + "##Attribute").c_str(), ImVec2(-1, 0)))
		{
			ImGui::OpenPopup("##Attributes");
		}

		if (ImGui::BeginPopup("##Attributes"))
		{
			std::vector<std::string> attributes;
			GEO::String::split_string(attribute_names(), ';', attributes);
			for (GEO::index_t i = 0; i < attributes.size(); ++i) 
			{
				if (ImGui::Button(attributes[i].c_str())) 
				{
					set_attribute(attributes[i]);
					ImGui::CloseCurrentPopup();
				}
			}
			ImGui::EndPopup();
		}

		ImGui::InputFloat("Min", &attribute_min_);
		ImGui::InputFloat("Max", &attribute_max_);
		if (ImGui::Button("Autorange", ImVec2(-1, 0))) 
		{
			autorange();
		}

		if (ImGui::ImageButton(convert_to_ImTextureID(current_colormap_texture_),ImVec2(115.0f * scaling(), 8.0f * scaling()))) 
		{
			ImGui::OpenPopup("##Colormap");
		}

		if (ImGui::BeginPopup("##Colormap"))
		{
			for (GEO::index_t i = 0; i < colormaps_.size(); ++i)
			{
				if (ImGui::ImageButton(convert_to_ImTextureID(colormaps_[i].texture), ImVec2(100.0f * scaling(), 8.0f * scaling())))
				{
					current_colormap_texture_ = colormaps_[i].texture;
					ImGui::CloseCurrentPopup();
				}
			}
			ImGui::EndPopup();
		}
	}

	if (mesh_.vertices.dimension() >= 6) 
	{
		ImGui::Separator();
		ImGui::Checkbox("Animate [a]", (bool*)glup_viewer_is_enabled_ptr(GLUP_VIEWER_IDLE_REDRAW));
		ImGui::SliderFloat("Speed", &anim_speed_, 1.0f, 10.0f, "%.1f");
		ImGui::SliderFloat("Time", &anim_time_, 0.0f, 1.0f, "%.2f");
	}

	ImGui::Separator();
	ImGui::Checkbox("##VertOnOff", &show_vertices_);
	ImGui::SameLine();
	ImGui::ColorEdit3WithPalette("Vertex[p]", vertices_color_.data());

	if (show_vertices_)
	{
		ImGui::Checkbox("Selection", &show_vertices_selection_);
		ImGui::SliderFloat("Size", &vertices_size_, 0.1f, 5.0f, "%.1f");
	}

	if (mesh_.facets.nb() != 0)
	{
		ImGui::Separator();
		ImGui::Checkbox("##SurfOnOff", &show_surface_);
		ImGui::SameLine();
		ImGui::ColorEdit3WithPalette("Surface [S]", surface_color_.data());
		if (show_surface_) 
		{
			ImGui::Checkbox("##SidesOnOff", &show_surface_sides_);
			ImGui::SameLine();
			ImGui::ColorEdit3WithPalette("2-Sided [c]", surface_color_2_.data());

			ImGui::Checkbox("##MeshOnOff", &show_mesh_);
			ImGui::SameLine();
			ImGui::ColorEdit3WithPalette("Mesh [m]", mesh_color_.data());

			if (show_mesh_)
			{
				ImGui::SliderFloat("Width##1", &mesh_width_, 0.1f, 1.0f, "%.1f");
			}

			ImGui::Checkbox("##BordersOnOff", &show_surface_borders_);
			ImGui::SameLine();
			ImGui::ColorEdit3WithPalette("Borders [B]", m_BorderColor.data());

			if (show_surface_borders_)
			{
				ImGui::SliderFloat("Width##2", &m_BorderWidth, 0.1f, 1.0f, "%.1f");
			}
		}
	}

	if (mesh_.cells.nb() != 0) 
	{
		ImGui::Separator();
		ImGui::Checkbox("##VolumeOnOff", &show_volume_);
		ImGui::SameLine();
		ImGui::ColorEdit3WithPalette("Volume [V]", volume_color_.data());
		if (show_volume_)
		{
			ImGui::SliderFloat("Shrink", &cells_shrink_, 0.0f, 1.0f, "%.2f");
			if (!mesh_.cells.are_simplices())
			{
				ImGui::Checkbox("Colored Cells [C]", &show_colored_cells_);
				ImGui::Checkbox("Hexes [j]", &show_hexes_);
			}
		}
	}
}

bool GioApplication::load(const std::string & Filename)
{
	CloseAlgorithmDialog();
	bool bResult = GEO::SimpleMeshApplication::load(Filename);
	mesh_.vertices.set_double_precision();
	m_MashFacetsAABB.reset(new GEO::MeshFacetsAABB(mesh_));
	return bResult;
}

bool GioApplication::save(const std::string & Filename)
{
	bool bResult = GEO::SimpleMeshApplication::save(Filename);
	return bResult;
}

int GioApplication::PANE_WIDTH() const
{
	return int(180 * scaling());;
}

void GioApplication::draw_application_menus()
{
	if (ImGui::BeginMenu("Algorithm"))
	{
		CloseAlgorithmDialog();
		
		ImGui::MenuItem("Cut Mesh", nullptr, &m_Algorithms[CUT_MESH_ALGORITHM_INDEX].bShowDialog);
		ImGui::Separator();
		ImGui::MenuItem("Barycentric Mapping", nullptr, &m_Algorithms[BARYCENTRIC_MAPPING_ALGORITHM_INDEX].bShowDialog);
		ImGui::MenuItem("LSCM", nullptr, &m_Algorithms[LSCM_ALGORITHM_INDEX].bShowDialog);
		ImGui::EndMenu();
	}

	DrawAlgorithmDialog();
}

void GioApplication::DrawAlgorithmDialog()
{
	DrawMeshCutAlgorithmDialog();
	DrawBarycentricMappingAlgorithmDialog();
	DrawLSCMAlgorithmDialog();
}

void GioApplication::CloseAlgorithmDialog()
{
	m_iCurrentAlgorithm = size_t(-1);
	CloseMeshCutAlgorithmDialog();
	CloseBarycentricMappingAlgorithmDialog();
	CloseLSCMAlgorithmDialog();
}

void GioApplication::DrawMeshCutAlgorithmDialog()
{
	if (m_Algorithms[CUT_MESH_ALGORITHM_INDEX].bShowDialog)
	{
		m_iCurrentAlgorithm = CUT_MESH_ALGORITHM_INDEX;

		ImGui::Begin("Cut Mesh", &m_Algorithms[CUT_MESH_ALGORITHM_INDEX].bShowDialog, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoCollapse);

		ImGui::Text("This algorithm cut an input mesh so that the mesh has only one boundary.");
		ImGui::TextColored(ImVec4(1.0f, 0.0f, 0.0f, 1.0f), "Note: Dose no effect on the mesh that has 0 boundary and 0 genus.\n");

		ImGui::Text("Select a start facet:");
		ImGui::Checkbox("Select Start Facet", &m_bSelectingFacet);
		ImGui::SameLine();
		ImGui::Text("Selected: %d", m_iSelectedFacet == GEO::NO_FACET ? 0 : m_iSelectedFacet);

		if (!m_bSelectingFacet)
		{
			ReleaseSelectingFacet(CUT_MESH_ALGORITHM_INDEX);
		}
		else
		{
			RequestSelectingFacet(CUT_MESH_ALGORITHM_INDEX, true);
		}

		if (ImGui::Button("Run Algorithm"))
		{
			GEO::index_t StartFacet = (m_iSelectedFacet == GEO::NO_FACET ? 0 : m_iSelectedFacet);
			m_Algorithms[CUT_MESH_ALGORITHM_INDEX].MeshAlgorithm->PutArg(MeshAlgorithm::MeshCutAlgorithm::PARAMS_KEY_START_FACET, StartFacet);
			m_Algorithms[CUT_MESH_ALGORITHM_INDEX].MeshAlgorithm->Execute(&mesh_);
			m_Algorithms[CUT_MESH_ALGORITHM_INDEX].bRunAlgorithm = true;
		}

		if (m_Algorithms[CUT_MESH_ALGORITHM_INDEX].bRunAlgorithm)
		{
			ImGui::SameLine();
			ImGui::Checkbox("Visualize", &m_Algorithms[CUT_MESH_ALGORITHM_INDEX].bVisualizeAlgorithm);
		}

		ImGui::End();
	}
	else
	{
		m_Algorithms[CUT_MESH_ALGORITHM_INDEX].bRunAlgorithm = false;
		m_Algorithms[CUT_MESH_ALGORITHM_INDEX].bVisualizeAlgorithm = false;
		ReleaseSelectingFacet(CUT_MESH_ALGORITHM_INDEX);
	}
}

void GioApplication::CloseMeshCutAlgorithmDialog()
{
	if (m_Algorithms[CUT_MESH_ALGORITHM_INDEX].bShowDialog)
	{
		m_Algorithms[CUT_MESH_ALGORITHM_INDEX].bRunAlgorithm = false;
		m_Algorithms[CUT_MESH_ALGORITHM_INDEX].bVisualizeAlgorithm = false;
		m_Algorithms[CUT_MESH_ALGORITHM_INDEX].bShowDialog = false;
		ReleaseSelectingFacet(CUT_MESH_ALGORITHM_INDEX);
	}
}

void GioApplication::DrawBarycentricMappingAlgorithmDialog()
{
	if (m_Algorithms[BARYCENTRIC_MAPPING_ALGORITHM_INDEX].bShowDialog)
	{
		m_iCurrentAlgorithm = BARYCENTRIC_MAPPING_ALGORITHM_INDEX;

		ImGui::Begin("Barycentric Mapping", &m_Algorithms[BARYCENTRIC_MAPPING_ALGORITHM_INDEX].bShowDialog, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoCollapse);
		std::vector<const char*> CoefficientTypes;
		for (size_t i = 0; i < IM_ARRAYSIZE(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE); i++)
		{
			CoefficientTypes.push_back(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_COEFFICIENT_TYPE[i].c_str());
		}

		static int iCoefficientTypesIdx = 0;
		ImGui::Combo("Coefficient Type", &iCoefficientTypesIdx, CoefficientTypes.data(), int(CoefficientTypes.size()));

		std::vector<const char*> DomainShapes;
		for (size_t i = 0; i < IM_ARRAYSIZE(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_DOMAIN_SHAPE); i++)
		{
			DomainShapes.push_back(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_DOMAIN_SHAPE[i].c_str());
		}

		static int iDomainShapesIdx = 0;
		ImGui::Combo("Domain Shape", &iDomainShapesIdx, DomainShapes.data(), int(DomainShapes.size()));

		std::vector<const char*> BoundaryFixWeights;
		for (size_t i = 0; i < IM_ARRAYSIZE(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_BOUNDARY_FIX_WEIGHT); i++)
		{
			BoundaryFixWeights.push_back(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_VALUE_SUPPORTED_BOUNDARY_FIX_WEIGHT[i].c_str());
		}

		static int iBoundaryFixWeight = 0;
		ImGui::Combo("Boundary Fix Weight", &iBoundaryFixWeight, BoundaryFixWeights.data(), int(BoundaryFixWeights.size()));

		static bool bFailure = false;
		if (ImGui::Button("Run Algorithm"))
		{
			std::string CoefficientType = CoefficientTypes[iCoefficientTypesIdx];
			std::string DomainShape = DomainShapes[iDomainShapesIdx];
			std::string BoundaryFixWeight = BoundaryFixWeights[iBoundaryFixWeight];
			m_Algorithms[BARYCENTRIC_MAPPING_ALGORITHM_INDEX].MeshAlgorithm->PutArg(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_KEY_COEFFICIENT_TYPE, CoefficientType);
			m_Algorithms[BARYCENTRIC_MAPPING_ALGORITHM_INDEX].MeshAlgorithm->PutArg(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_KEY_DOMAIN_SHAPE, DomainShape);
			m_Algorithms[BARYCENTRIC_MAPPING_ALGORITHM_INDEX].MeshAlgorithm->PutArg(MeshAlgorithm::BarycentricMappingAlgorithm::PARAMS_KEY_BOUNDARY_FIX_WEIGHT, BoundaryFixWeight);
			if (m_Algorithms[BARYCENTRIC_MAPPING_ALGORITHM_INDEX].MeshAlgorithm->Execute(&mesh_))
			{
				m_Algorithms[BARYCENTRIC_MAPPING_ALGORITHM_INDEX].bRunAlgorithm = true;
				bFailure = false;
			}
			else
			{
				ImGui::OpenPopup("Error##Barycentric Mapping");
				bFailure = true;
			}
		}

		if (bFailure)
		{
			if (ImGui::BeginPopupModal("Error##Barycentric Mapping"), ImGuiWindowFlags_AlwaysAutoResize)
			{
				ImGui::Text("Some errors occurred. Input mesh must have 1 boundary and 0 genus!");
				if (ImGui::Button("OK", ImVec2(120, 0)))
				{
					ImGui::CloseCurrentPopup();
					bFailure = false;
				}
				ImGui::EndPopup();
			}
		}

		if (m_Algorithms[BARYCENTRIC_MAPPING_ALGORITHM_INDEX].bRunAlgorithm)
		{
			ImGui::SameLine();
			ImGui::TextColored(ImVec4(1.0f, 0.0f, 0.0f, 1.0f), "Open Animation to visualize the algorithm");
		}

		ImGui::End();
	}
	else
	{
		m_Algorithms[BARYCENTRIC_MAPPING_ALGORITHM_INDEX].bVisualizeAlgorithm = false;
		m_Algorithms[BARYCENTRIC_MAPPING_ALGORITHM_INDEX].bRunAlgorithm = false;
	}
}

void GioApplication::CloseBarycentricMappingAlgorithmDialog()
{
	if (m_Algorithms[BARYCENTRIC_MAPPING_ALGORITHM_INDEX].bShowDialog)
	{
		m_Algorithms[BARYCENTRIC_MAPPING_ALGORITHM_INDEX].bRunAlgorithm = false;
		m_Algorithms[BARYCENTRIC_MAPPING_ALGORITHM_INDEX].bVisualizeAlgorithm = false;
		m_Algorithms[BARYCENTRIC_MAPPING_ALGORITHM_INDEX].bShowDialog = false;
	}
}

void GioApplication::DrawLSCMAlgorithmDialog()
{
	if (m_Algorithms[LSCM_ALGORITHM_INDEX].bShowDialog)
	{
		m_iCurrentAlgorithm = LSCM_ALGORITHM_INDEX;

		ImGui::Begin("LSCM", &m_Algorithms[LSCM_ALGORITHM_INDEX].bShowDialog, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoCollapse);

		static bool bFailure = false;
		if (ImGui::Button("Run Algorithm"))
		{
			if (m_Algorithms[LSCM_ALGORITHM_INDEX].MeshAlgorithm->Execute(&mesh_))
			{
				m_Algorithms[LSCM_ALGORITHM_INDEX].bRunAlgorithm = true;
				bFailure = false;
			}
			else
			{
				ImGui::OpenPopup("Error##LSCM");
				bFailure = true;
			}
		}

		if (bFailure)
		{
			if (ImGui::BeginPopupModal("Error##LSCM"), ImGuiWindowFlags_AlwaysAutoResize)
			{
				ImGui::Text("Some errors occurred. Input mesh must have 1 boundary and 0 genus!");
				if (ImGui::Button("OK", ImVec2(120, 0)))
				{
					ImGui::CloseCurrentPopup();
					bFailure = false;
				}
				ImGui::EndPopup();
			}
		}

		if (m_Algorithms[LSCM_ALGORITHM_INDEX].bRunAlgorithm)
		{
			ImGui::SameLine();
			ImGui::TextColored(ImVec4(1.0f, 0.0f, 0.0f, 1.0f), "Open Animation to visualize the algorithm");
		}

		ImGui::End();
	}
	else
	{
		m_Algorithms[LSCM_ALGORITHM_INDEX].bVisualizeAlgorithm = false;
		m_Algorithms[LSCM_ALGORITHM_INDEX].bRunAlgorithm = false;
	}
}

void GioApplication::CloseLSCMAlgorithmDialog()
{
	if (m_Algorithms[LSCM_ALGORITHM_INDEX].bShowDialog)
	{
		m_Algorithms[LSCM_ALGORITHM_INDEX].bRunAlgorithm = false;
		m_Algorithms[LSCM_ALGORITHM_INDEX].bVisualizeAlgorithm = false;
		m_Algorithms[LSCM_ALGORITHM_INDEX].bShowDialog = false;
	}
}

void GioApplication::RequestSelectingFacet(size_t iAlgorithm, bool bSingleSelect)
{
	assert(m_iFacetSelectingAlgorithm == size_t(-1) || m_iFacetSelectingAlgorithm == iAlgorithm);
	m_iFacetSelectingAlgorithm = iAlgorithm;
	m_bDisplayAllSelectedFacets = !bSingleSelect;
	m_bSelectingFacet = true;
}

void GioApplication::RequestSelectingVertex(size_t iAlgorithm, bool bSingleSelect)
{
	assert(m_iVertexSelectingAlgorithm == size_t(-1) || m_iVertexSelectingAlgorithm == iAlgorithm);
	m_iVertexSelectingAlgorithm = iAlgorithm;
	m_bDisplayAllSelectedVertices = !bSingleSelect;
	m_bSelectingVertex = true;
}

void GioApplication::ReleaseSelectingFacet(size_t iAlgorithm)
{
	if (m_iFacetSelectingAlgorithm == iAlgorithm)
	{
		m_iSelectedFacet = GEO::NO_FACET;
		m_iSelectedFacets.clear();
		m_iFacetSelectingAlgorithm = size_t(-1);
		m_bSelectingFacet = false;
		m_bDisplayAllSelectedFacets = false;
	}
}

void GioApplication::ReleaseSelectingVertex(size_t iAlgorithm)
{
	if (m_iFacetSelectingAlgorithm == iAlgorithm)
	{
		m_iSelectedVertex = GEO::NO_VERTEX;
		m_iSelectedVertices.clear();
		m_iVertexSelectingAlgorithm = size_t(-1);
		m_bSelectingVertex = false;
		m_bDisplayAllSelectedVertices = false;
	}
}

GLboolean MouseCallbackFunc(float X, float Y, int Button, enum GlupViewerEvent Event)
{
	double HitPoint[3];
	GLboolean bHitBackground;
	GEO::vec3 NearestPoint;
	double SqDistance;
	GEO::index_t iFacet;

	switch (Event)
	{
	case GLUP_VIEWER_DOWN:
		if (g_pApp->m_bSelectingFacet)
		{
			glup_viewer_get_picked_point(HitPoint, &bHitBackground);
			iFacet = g_pApp->m_MashFacetsAABB->nearest_facet(GEO::vec3(HitPoint[0], HitPoint[1], HitPoint[2]), NearestPoint, SqDistance);

			if (SqDistance < 0.05 && Button == 0)
			{
				g_pApp->m_iSelectedFacets.insert(iFacet);
				g_pApp->m_iSelectedFacet = iFacet;
			}

			if (SqDistance < 0.05 && Button == 1)
			{
				if (iFacet == g_pApp->m_iSelectedFacet)
				{
					g_pApp->m_iSelectedFacet = GEO::NO_FACET;
				}
				g_pApp->m_iSelectedFacets.erase(iFacet);
			}
		}

		if (g_pApp->m_bSelectingVertex)
		{
			glup_viewer_get_picked_point(HitPoint, &bHitBackground);
			iFacet = g_pApp->m_MashFacetsAABB->nearest_facet(GEO::vec3(HitPoint[0], HitPoint[1], HitPoint[2]), NearestPoint, SqDistance);

			double MinSqDistance = std::numeric_limits<double>::max();
			GEO::index_t iMinDistanceVertex = GEO::NO_VERTEX;

			for (GEO::index_t i = 0; i < g_pApp->mesh_.facets.nb_vertices(iFacet); i++)
			{
				GEO::index_t iVertex = g_pApp->mesh_.facets.vertex(iFacet, i);
				GEO::vec3 Vertex = g_pApp->mesh_.vertices.point(iVertex);
				double SqDistance = std::pow((Vertex.x - NearestPoint.x), 2.0) + std::pow((Vertex.y - NearestPoint.y), 2.0) + std::pow((Vertex.z - NearestPoint.z), 2.0);
				if (SqDistance < MinSqDistance)
				{
					MinSqDistance = SqDistance;
					iMinDistanceVertex = iVertex;
				}
			}

			if (MinSqDistance < 1e-5 && Button == 0)
			{
				g_pApp->m_iSelectedVertices.insert(iMinDistanceVertex);
				g_pApp->m_iSelectedVertex = iMinDistanceVertex;
			}

			if (MinSqDistance < 1e-5 && Button == 1)
			{
				if (iMinDistanceVertex == g_pApp->m_iSelectedVertex)
				{
					g_pApp->m_iSelectedVertex = GEO::NO_VERTEX;
				}
				g_pApp->m_iSelectedVertices.erase(iMinDistanceVertex);
			}
		}
	case GLUP_VIEWER_MOVE:
		break;
	case GLUP_VIEWER_UP:
		break;
	}
	return false;
}
