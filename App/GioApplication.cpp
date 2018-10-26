#include "GioApplication.hpp"

static MeshAlgorithm::MeshCutAlgorithm MeshCut;

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

	if (m_bSelectingFacet && m_iSelectedFacet != GEO::NO_FACET)
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

	if (m_bVisualizeAlgorithm && m_MeshAlgorithm != nullptr)
	{
		m_MeshAlgorithm->Visualize(&mesh_);
	}
}

void GioApplication::init_graphics()
{
	GEO::SimpleMeshApplication::init_graphics();
	GEO::set_assert_mode(GEO::ASSERT_BREAKPOINT);
	glup_viewer_set_mouse_func(MouseCallbackFunc);
	retina_mode_ = false;
	scaling_ = 1.0f;
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

int GioApplication::PANE_WIDTH() const
{
	return int(180 * scaling());;
}

void GioApplication::draw_application_menus()
{
	if (ImGui::BeginMenu("Algorithm"))
	{
		ImGui::MenuItem("Slice Mesh", nullptr, &m_bShowMeshCutAlgorithmDialog);
		ImGui::EndMenu();
	}

	DrawAlgorithmDialog();
}

void GioApplication::DrawAlgorithmDialog()
{
	DrawMeshCutAlgorithmDialog();
}

void GioApplication::CloseAlgorithmDialog()
{
	CloseMeshCutAlgorithmDialog();
}

void GioApplication::DrawMeshCutAlgorithmDialog()
{
	if (m_bShowMeshCutAlgorithmDialog)
	{
		if (m_MeshAlgorithm == nullptr)
		{
			m_MeshAlgorithm.reset(new MeshAlgorithm::MeshCutAlgorithm());
		}

		ImGui::Begin("Cut Mesh", &m_bShowMeshCutAlgorithmDialog, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoCollapse);

		ImGui::Text("This algorithm cut an input mesh so that the mesh has only one boundary.");
		ImGui::TextColored(ImVec4(1.0f, 0.0f, 0.0f, 1.0f), "Note: Dose no effect on the mesh that has 0 boundary and 0 genus.\n");

		ImGui::Text("Select a start facet:");
		ImGui::Checkbox("Select Start Facet", &m_bSelectingFacet);
		ImGui::SameLine();
		ImGui::Text("Selected: %d", m_iSelectedFacet == GEO::NO_FACET ? 0 : m_iSelectedFacet);

		if (!m_bSelectingFacet)
		{
			m_iSelectedFacet = GEO::NO_FACET;
		}

		if (ImGui::Button("Run Algorithm"))
		{
			GEO::index_t StartFacet = (m_iSelectedFacet == GEO::NO_FACET ? 0 : m_iSelectedFacet);
			m_MeshAlgorithm->PutArg(MeshAlgorithm::MeshCutAlgorithm::PARAMS_KEY_START_FACET, StartFacet);
			m_MeshAlgorithm->Execute(&mesh_);
			m_bRunAlgorithm = true;
		}

		if (m_bRunAlgorithm)
		{
			ImGui::SameLine();
			ImGui::Checkbox("Visualize", &m_bVisualizeAlgorithm);
		}

		ImGui::End();
	}
	else
	{
		m_MeshAlgorithm.release();
		m_bRunAlgorithm = false;
		m_bVisualizeAlgorithm = false;
		m_bSelectingFacet = false;
		m_iSelectedFacet = GEO::NO_FACET;
	}
}

void GioApplication::CloseMeshCutAlgorithmDialog()
{
	if (m_bShowMeshCutAlgorithmDialog)
	{
		m_MeshAlgorithm.release();
		m_bRunAlgorithm = false;
		m_bVisualizeAlgorithm = false;
		m_bSelectingFacet = false;
		m_iSelectedFacet = GEO::NO_FACET;

		m_bShowMeshCutAlgorithmDialog = false;
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
				g_pApp->m_iSelectedFacet = iFacet;
			}
			
			if (SqDistance < 0.05 && Button == 1)
			{
				if (iFacet == g_pApp->m_iSelectedFacet)
				{
					g_pApp->m_iSelectedFacet = GEO::NO_FACET;
				}
			}

		}
	case GLUP_VIEWER_MOVE:
		break;
	case GLUP_VIEWER_UP:
		break;
	}
	return false;
}
