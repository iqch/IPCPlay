#include "pxr/pxr.h"
#include "pxr/usdImaging/usdImaging/api.h"
#include "pxr/usdImaging/usdImaging/primAdapter.h"
#include "pxr/usdImaging/usdImaging/gprimAdapter.h"

//#include <iostream>

//#include "pxr/usdImaging/usdImaging/tetraMeshAdapter.h"

#include "pxr/usdImaging/usdImaging/debugCodes.h"
#include "pxr/usdImaging/usdImaging/delegate.h"
#include "pxr/usdImaging/usdImaging/indexProxy.h"
#include "pxr/usdImaging/usdImaging/tokens.h"

#include "pxr/imaging/hd/mesh.h"
#include "pxr/imaging/hd/geomSubset.h"
#include "pxr/imaging/hd/perfLog.h"

#include "pxr/imaging/pxOsd/meshTopology.h"
#include "pxr/imaging/pxOsd/tokens.h"

#include "pxr/usd/usdGeom/mesh.h"
#include "pxr/usd/usdGeom/primvarsAPI.h"
#include "pxr/usd/usdGeom/subset.h"
#include "pxr/usd/usdGeom/xformCache.h"

#include "pxr/base/tf/type.h"

PXR_NAMESPACE_OPEN_SCOPE



/// \class UsdImagingTetraMeshAdapter
///
/// Delegate support for UsdGeomMesh. // SAME SIMPLIFIED MESH-ADAPTOR WITH CUSTOMIZED TOPOLOGY
///
class UsdImagingTetraMeshAdapter : public UsdImagingGprimAdapter
{
public:
    using BaseAdapter = UsdImagingGprimAdapter;

    UsdImagingTetraMeshAdapter() : UsdImagingGprimAdapter() {};
    USDIMAGING_API ~UsdImagingTetraMeshAdapter() override {};

    USDIMAGING_API
        SdfPath Populate(
            UsdPrim const& prim,
            UsdImagingIndexProxy* index,
            UsdImagingInstancerContext const* instancerContext = nullptr) override;

    USDIMAGING_API
        bool IsSupported(UsdImagingIndexProxy const* index) const override;

    /// Thread Safe.
    USDIMAGING_API
        void TrackVariability(UsdPrim const& prim,
            SdfPath const& cachePath,
            HdDirtyBits* timeVaryingBits,
            UsdImagingInstancerContext const*
            instancerContext = nullptr) const override;


    /// Thread Safe.
    USDIMAGING_API // OBSOLETE
    void UpdateForTime(UsdPrim const& prim,
                       SdfPath const& cachePath, 
                       UsdTimeCode time,
                       HdDirtyBits requestedBits,
                       UsdImagingInstancerContext const* 
                           instancerContext = nullptr) const override;


 

    // ---------------------------------------------------------------------- //
    /// \name Change Processing
    // ---------------------------------------------------------------------- //

    USDIMAGING_API
        HdDirtyBits ProcessPropertyChange(UsdPrim const& prim,
            SdfPath const& cachePath,
            TfToken const& propertyName) override;

    // ---------------------------------------------------------------------- //
    /// \name Data access
    // ---------------------------------------------------------------------- //

    USDIMAGING_API
        //VtValue GetTopology(UsdPrim const& prim,
        void GetTopology(UsdPrim const& prim,
            //SdfPath const& cachePath,
            VtValue* topo,
            UsdTimeCode time) const; // NEW override;

    //USDIMAGING_API - OBSOLETE 
    //    VtValue Get(UsdPrim const& prim,
    //        SdfPath const& cachePath,
    //        TfToken const& key,
    //        UsdTimeCode time,
    //        VtIntArray* outIndices) const; // NEW override;

protected:
};


TF_REGISTRY_FUNCTION(TfType) // +
{
    typedef UsdImagingTetraMeshAdapter Adapter;
    TfType t = TfType::Define<Adapter, TfType::Bases<Adapter::BaseAdapter> >();
    t.SetFactory< UsdImagingPrimAdapterFactory<Adapter> >();
};

//UsdImagingTetraMeshAdapter::~UsdImagingTetraMeshAdapter() {};

// +
bool UsdImagingTetraMeshAdapter::IsSupported(UsdImagingIndexProxy const* index) const
{
    return index->IsRprimTypeSupported(HdPrimTypeTokens->mesh);
}

SdfPath // +
UsdImagingTetraMeshAdapter::Populate(UsdPrim const& prim,
    UsdImagingIndexProxy* index,
    UsdImagingInstancerContext const* instancerContext)
{
    SdfPath cachePath = _AddRprim(HdPrimTypeTokens->mesh, prim, index, GetMaterialUsdPath(prim), instancerContext);
    return cachePath;
}

//TfToken tetraVertexIndices("tetraVertexIndices");

void // +-
UsdImagingTetraMeshAdapter::TrackVariability(UsdPrim const& prim,
    SdfPath const& cachePath,
    HdDirtyBits* timeVaryingBits,
    UsdImagingInstancerContext const*
    instancerContext) const
{
    BaseAdapter::TrackVariability(prim, cachePath, timeVaryingBits, instancerContext);

    // WARNING: This method is executed from multiple threads, the value cache
    // has been carefully pre-populated to avoid mutating the underlying
    // container during update.

    // Discover time-varying points.
    _IsVarying(prim,
        UsdGeomTokens->points,
        HdChangeTracker::DirtyPoints,
        UsdImagingTokens->usdVaryingPrimvar,
        timeVaryingBits,
        /*isInherited*/false);

     _IsVarying(prim,
        /*UsdGeomTokens->*/TfToken("tetraVertexIndices"),
        HdChangeTracker::DirtyTopology,
        UsdImagingTokens->usdVaryingTopology,
        timeVaryingBits,
        /*isInherited*/false);

    _IsVarying(prim,
        UsdGeomTokens->orientation,
        HdChangeTracker::DirtyTopology,
        UsdImagingTokens->usdVaryingTopology,
        timeVaryingBits,
        /*isInherited*/false);
}

HdDirtyBits //+
UsdImagingTetraMeshAdapter::ProcessPropertyChange(UsdPrim const& prim,
    SdfPath const& cachePath,
    TfToken const& propertyName)
{
    if (propertyName == UsdGeomTokens->points)
        return HdChangeTracker::DirtyPoints;

    if (propertyName == /*UsdGeomTokens->*/TfToken("tetraVertexIndices") ||
        propertyName == UsdGeomTokens->orientation)
    {
        return HdChangeTracker::DirtyTopology;
    }

    // Allow base class to handle change processing.
    return BaseAdapter::ProcessPropertyChange(prim, cachePath, propertyName);
}

void
UsdImagingTetraMeshAdapter::UpdateForTime(UsdPrim const& prim,
    SdfPath const& cachePath,
    UsdTimeCode time,
    HdDirtyBits requestedBits,
    UsdImagingInstancerContext const*
    instancerContext) const
{
    //TF_DEBUG(USDIMAGING_CHANGES).Msg("[UpdateForTime] Mesh path: <%s>\n",
    //    prim.GetPath().GetText());

    BaseAdapter::UpdateForTime(
        prim, cachePath, time, requestedBits, instancerContext);

    UsdImagingValueCache* valueCache = _GetValueCache();
    HdPrimvarDescriptorVector& primvars = valueCache->GetPrimvars(cachePath);

    if (requestedBits & HdChangeTracker::DirtyTopology) {
        VtValue& topology = valueCache->GetTopology(cachePath);
        GetTopology(prim, &topology, time);
    }
}

/*virtual*/
/*VtValue*/ void
UsdImagingTetraMeshAdapter::GetTopology(UsdPrim const& prim,
    //SdfPath const& cachePath,
    VtValue* topo,
    UsdTimeCode time) const
{
    TRACE_FUNCTION();
    HF_MALLOC_TAG_FUNCTION();

    //TfToken schemeToken = PxOsdOpenSubdivTokens->none;
    //_GetPtr(prim, UsdGeomTokens->subdivisionScheme, time, &schemeToken);

    typedef VtArray<GfVec4i> VtInt4Array; // NEWLY IMPLEMENTED CODE FOR BUILDIN A BUNCH OF TETRAHEDRON
    VtInt4Array TOPO;
    //= _Get<VtInt4Array>(prim, /*UsdGeomTokens->*/tetraVertexIndices, time);

    prim.GetAttribute(TfToken("tetraVertexIndices")).Get<VtInt4Array>(&TOPO, time);

    VtIntArray T(TOPO.size() * 12);

    // ...POSSIBLE PLAY WITH TBB LATER

    int index = 0;
    for (GfVec4i t : TOPO)
    {
        T[index] = t[0]; index++;
        T[index] = t[1]; index++;
        T[index] = t[2]; index++;

        T[index] = t[1]; index++;
        T[index] = t[0]; index++;
        T[index] = t[3]; index++;

        T[index] = t[2]; index++;
        T[index] = t[1]; index++;
        T[index] = t[3]; index++;

        T[index] = t[0]; index++;
        T[index] = t[2]; index++;
        T[index] = t[3]; index++;
    };

    VtIntArray F(TOPO.size() * 4, 3);

    HdMeshTopology meshTopo(
        PxOsdOpenSubdivTokens->none,
        _Get<TfToken>(prim, UsdGeomTokens->orientation, time),
        F, T, VtIntArray());

    topo->Swap(meshTopo);
    // return VtValue(meshTopo); // return that mesh!
    // NOTE : I DON'T CHECK ADJACENT FACES DUPLICATION - I SEEM, IT IS BETTER TO RENDER AND CULL THEM LTER WITH RENDERER
}   // YOU MAY IMPLEMENT IT IN DIFFERENT WAY

/*virtual*/
//VtValue // +
//UsdImagingTetraMeshAdapter::Get(UsdPrim const& prim,
//    SdfPath const& cachePath,
//    TfToken const& key,
//    UsdTimeCode time,
//    VtIntArray* outIndices) const
//{
//    TRACE_FUNCTION();
//    HF_MALLOC_TAG_FUNCTION();
//
//    return BaseAdapter::Get(prim, cachePath, key, time, outIndices);
//}

PXR_NAMESPACE_CLOSE_SCOPE
