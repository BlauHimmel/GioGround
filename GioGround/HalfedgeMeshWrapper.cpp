#include "HalfedgeMeshWrapper.hpp"

HalfedgeMeshWrapper::HalfedgeMeshWrapper(Ref GEO::Mesh * pMesh) : pMesh(pMesh)
{
	Reset();
}

void HalfedgeMeshWrapper::Reset()
{
	Corner2Corner.resize(pMesh->facet_corners.nb(), GEO::max_index_t());
	Corner2Facet.resize(pMesh->facet_corners.nb(), GEO::max_index_t());
	Vertex2Corner.resize(pMesh->vertices.nb(), GEO::max_index_t());

	for (GEO::index_t iFacet = 0; iFacet < pMesh->facets.nb(); iFacet++)
	{
		for (GEO::index_t iCorner = pMesh->facets.corners_begin(iFacet); iCorner != pMesh->facets.corners_end(iFacet); iCorner++)
		{
			Corner2Corner[iCorner] = iCorner;
			Corner2Facet[iCorner] = iFacet;
			Vertex2Corner[pMesh->facet_corners.vertex(iCorner)] = iCorner;
		}
	}

	for (GEO::index_t iFacet = 0; iFacet < pMesh->facets.nb(); iFacet++)
	{
		for (GEO::index_t iCorner = pMesh->facets.corners_begin(iFacet); iCorner != pMesh->facets.corners_end(iFacet); iCorner++)
		{
			Corner2Corner[iCorner] = Vertex2Corner[pMesh->facet_corners.vertex(iCorner)];
			Vertex2Corner[pMesh->facet_corners.vertex(iCorner)] = iCorner;
		}
	}
}

GEO::index_t HalfedgeMeshWrapper::Next(In GEO::index_t iCorner)
{
	return pMesh->facets.next_corner_around_facet(Corner2Facet[iCorner], iCorner);
}

GEO::index_t HalfedgeMeshWrapper::Prev(In GEO::index_t iCorner)
{
	return pMesh->facets.prev_corner_around_facet(Corner2Facet[iCorner], iCorner);
}

GEO::index_t HalfedgeMeshWrapper::Opposite(In GEO::index_t iCorner)
{
	GEO::index_t AdjFacet = pMesh->facet_corners.adjacent_facet(iCorner);

	if (AdjFacet == GEO::NO_FACET)
	{
		return GEO::NO_CORNER;
	}

	GEO::index_t V = pMesh->facet_corners.vertex(iCorner);

	for (GEO::index_t iCorner = pMesh->facets.corners_begin(AdjFacet); iCorner != pMesh->facets.corners_end(AdjFacet); iCorner++)
	{
		if (pMesh->facet_corners.vertex(Next(iCorner)) == V)
		{
			return iCorner;
		}
	}
	
	return GEO::NO_CORNER;
}

GEO::index_t HalfedgeMeshWrapper::Facet(In GEO::index_t iCorner)
{
	return Corner2Facet[iCorner];
}

GEO::index_t HalfedgeMeshWrapper::NextAroundVertex(In GEO::index_t iCorner)
{
	return Opposite(Prev(iCorner));
}

GEO::index_t HalfedgeMeshWrapper::Origin(In GEO::index_t iCorner)
{
	return pMesh->facet_corners.vertex(iCorner);
}

GEO::index_t HalfedgeMeshWrapper::Dest(In GEO::index_t iCorner)
{
	return pMesh->facet_corners.vertex(Next(iCorner));
}
