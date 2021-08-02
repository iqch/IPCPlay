//
// USD LIB
#include "pxr/pxr.h"
#include "pxr/usd/usdGeom/api.h"
#include "pxr/usd/usdGeom/pointBased.h"
#include "pxr/usd/usd/prim.h"
#include "pxr/usd/usd/stage.h"
#include "pxr/usd/usdGeom/tokens.h"

#include "pxr/usd/usd/timeCode.h" 

#include "pxr/base/vt/value.h"

#include "pxr/base/gf/vec3d.h"
#include "pxr/base/gf/vec3f.h"
#include "pxr/base/gf/vec4i.h"
#include "pxr/base/gf/matrix4d.h"

#include "pxr/base/tf/token.h"
#include "pxr/base/tf/type.h"

// 
//#include "pxr/usd/usdGeom/tetramesh.h"
#include "pxr/usd/usd/schemaRegistry.h"
#include "pxr/usd/usd/typed.h"

#include "pxr/usd/sdf/types.h"
#include "pxr/usd/sdf/assetPath.h"

PXR_NAMESPACE_OPEN_SCOPE

class SdfAssetPath;

class UsdGeomTetraMesh : public UsdGeomPointBased
{
public:
    static const UsdSchemaType schemaType = UsdSchemaType::ConcreteTyped; // 20.08
    //static const UsdSchemaKind schemaKind = UsdSchemaKind::ConcreteTyped; // 21.08
    explicit UsdGeomTetraMesh(const UsdPrim& prim = UsdPrim())
        : UsdGeomPointBased(prim) {};

    explicit UsdGeomTetraMesh(const UsdSchemaBase& schemaObj)
        : UsdGeomPointBased(schemaObj) {};

    //USDGEOM_API virtual ~UsdGeomTetraMesh() {};

    USDGEOM_API
        static const TfTokenVector&
        GetSchemaAttributeNames(bool includeInherited = true);

    USDGEOM_API
        static UsdGeomTetraMesh
        Get(const UsdStagePtr& stage, const SdfPath& path);

    USDGEOM_API
        static UsdGeomTetraMesh
        Define(const UsdStagePtr& stage, const SdfPath& path);

protected:
    USDGEOM_API UsdSchemaType _GetSchemaType() const override; // 20.08
    //USDGEOM_API UsdSchemaKind _GetSchemaKind() const override; //21.08

private:
    friend class UsdSchemaRegistry;
    USDGEOM_API static const TfType& _GetStaticTfType();

    static bool _IsTypedSchema();

    // override SchemaBase virtuals.
    USDGEOM_API const TfType& _GetTfType() const override;

public:
    // --------------------------------------------------------------------- //
    // TETRAEVERTEXINDICES 
    // --------------------------------------------------------------------- //
    ///
    /// | ||
    /// | -- | -- |
    /// | Declaration | `int4[] tetraVertexIndices` |
    /// | C++ Type | VtArray<int4> |
    /// | \ref Usd_Datatypes "Usd Type" | SdfValueTypeNames->Int4Array |
    USDGEOM_API
        UsdAttribute GetTetraVertexIndicesAttr() const;

    USDGEOM_API
        UsdAttribute CreateTetraVertexIndicesAttr(
            VtValue const& defaultValue, bool writeSparsely) const;

    USDGEOM_API
        size_t GetTetraCount(UsdTimeCode timeCode = UsdTimeCode::Default()) const;

};

// Register the schema with the TfType system.
TF_REGISTRY_FUNCTION(TfType)
{
    TfType::Define<UsdGeomTetraMesh, TfType::Bases< UsdGeomPointBased > >();

    // Register the usd prim typename as an alias under UsdSchemaBase. This
    // enables one to call
    // TfType::Find<UsdSchemaBase>().FindDerivedByName("Mesh")
    // to find TfType<UsdGeomMesh>, which is how IsA queries are
    // answered.
    TfType::AddAlias<UsdSchemaBase, UsdGeomTetraMesh>("TetraMesh");
};

/* static */
UsdGeomTetraMesh
UsdGeomTetraMesh::Get(const UsdStagePtr& stage, const SdfPath& path)
{
    if (!stage) {
        TF_CODING_ERROR("Invalid stage");
        return UsdGeomTetraMesh();
    }
    return UsdGeomTetraMesh(stage->GetPrimAtPath(path));
};

/* static */
UsdGeomTetraMesh
UsdGeomTetraMesh::Define(
    const UsdStagePtr& stage, const SdfPath& path)
{
    static TfToken usdPrimTypeName("TetraMesh");
    if (!stage) {
        TF_CODING_ERROR("Invalid stage");
        return UsdGeomTetraMesh();
    };
    return UsdGeomTetraMesh(stage->DefinePrim(path, usdPrimTypeName));
};

/* virtual */
UsdSchemaType UsdGeomTetraMesh::_GetSchemaType() const
{
    return UsdGeomTetraMesh::schemaType;
};

/* static */
const TfType&
UsdGeomTetraMesh::_GetStaticTfType()
{
    static TfType tfType = TfType::Find<UsdGeomTetraMesh>();
    return tfType;
};

/* static */
bool
UsdGeomTetraMesh::_IsTypedSchema()
{
    static bool isTyped = _GetStaticTfType().IsA<UsdTyped>();
    return isTyped;
};

/* virtual */
const TfType&
UsdGeomTetraMesh::_GetTfType() const
{
    return _GetStaticTfType();
};

/// ATTRIBUTES

static TfToken tetraVertexIndices("tetraVertexIndices");

UsdAttribute
UsdGeomTetraMesh::GetTetraVertexIndicesAttr() const
{
    return GetPrim().GetAttribute(tetraVertexIndices);
};

UsdAttribute
UsdGeomTetraMesh::CreateTetraVertexIndicesAttr(VtValue const& defaultValue, bool writeSparsely) const
{
    return UsdSchemaBase::_CreateAttr(
        tetraVertexIndices, SdfValueTypeNames->Int4Array,
        /* custom = */ false,
        SdfVariabilityUniform,
        defaultValue, writeSparsely);
};


namespace {
    static inline TfTokenVector
        _ConcatenateAttributeNames(const TfTokenVector& left, const TfTokenVector& right)
    {
        TfTokenVector result;
        result.reserve(left.size() + right.size());
        result.insert(result.end(), left.begin(), left.end());
        result.insert(result.end(), right.begin(), right.end());
        return result;
    }
}

/*static*/
const TfTokenVector&
UsdGeomTetraMesh::GetSchemaAttributeNames(bool includeInherited)
{
    static TfTokenVector localNames = {
        tetraVertexIndices,
    };
    static TfTokenVector allNames =
        _ConcatenateAttributeNames(
            UsdGeomPointBased::GetSchemaAttributeNames(true),
            localNames);

    if (includeInherited)
        return allNames;
    else
        return localNames;
}

//PXR_NAMESPACE_CLOSE_SCOPE

// ===================================================================== //
// Feel free to add custom code below this line. It will be preserved by
// the code generator.
//
// Just remember to wrap code in the appropriate delimiters:
// 'PXR_NAMESPACE_OPEN_SCOPE', 'PXR_NAMESPACE_CLOSE_SCOPE'.
// ===================================================================== //
// --(BEGIN CUSTOM CODE)--

//#include "pxr/base/arch/hints.h"
//#include "pxr/base/tf/stringUtils.h"
//
//#include <numeric>

size_t
UsdGeomTetraMesh::GetTetraCount(UsdTimeCode timeCode) const
{
    UsdAttribute tetCountsAttr = GetTetraVertexIndicesAttr();
    VtVec4iArray tetCounts;
    tetCountsAttr.Get(&tetCounts, timeCode);
    return tetCounts.size();
}


/// PY
#include "pxr/pxr.h"
#include "pxr/base/tf/registryManager.h"
#include "pxr/base/tf/scriptModuleLoader.h"
#include "pxr/base/tf/token.h"

#include <vector>

//#include "pxr/base/tf/pySafePython.h"
#include "pxr/pxr.h"
//#include "pxr/base/tf/pyModule.h"

#include "pxr/usd/usd/schemaBase.h"

#include "pxr/usd/sdf/primSpec.h"

//#include "pxr/usd/usd/pyConversions.h"
//#include "pxr/base/tf/pyContainerConversions.h"
//#include "pxr/base/tf/pyResultConversions.h"
//#include "pxr/base/tf/pyUtils.h"
//#include "pxr/base/tf/wrapTypeHelpers.h"

//#include <hboost/python.hpp>

#include <string>


//TF_REGISTRY_FUNCTION(TfScriptModuleLoader) 
//{
//    // List of direct dependencies for this library.
//    const std::vector<TfToken> reqs = {
//        TfToken("sdf"),
//        TfToken("tf"),
//        TfToken("usd"),
//        TfToken("vt"),
//        TfToken("usdGeom"),
//    };
//    TfScriptModuleLoader::GetInstance().
//        RegisterLibrary(TfToken("usdMore"), TfToken("pxr.UsdMore"), reqs);
//}

PXR_NAMESPACE_CLOSE_SCOPE


//PXR_NAMESPACE_USING_DIRECTIVE

//TF_WRAP_MODULE
//{
//    TF_WRAP(UsdSchemaExamplesSimple);
//    TF_WRAP(UsdSchemaExamplesComplex);
//    TF_WRAP(UsdSchemaExamplesParamsAPI);
//    TF_WRAP(UsdSchemaExamplesTokens);
//}


//using namespace hboost::python;

PXR_NAMESPACE_USING_DIRECTIVE

namespace {
//
//#define WRAP_CUSTOM                                                     \
//    template <class Cls> static void _CustomWrapCode(Cls &_class)
//
//    // fwd decl.
//    WRAP_CUSTOM;
//
//
//    static UsdAttribute
//        _CreateTetraVertexIndicesAttr(UsdGeomTetraMesh& self,
//            object defaultVal, bool writeSparsely)
//    {
//        return self.CreateTetraVertexIndicesAttr(
//            UsdPythonToSdfType(defaultVal, SdfValueTypeNames->Int4Array), writeSparsely);
//    }
//
//    static std::string
//        _Repr(const UsdGeomTetraMesh& self)
//    {
//        std::string primRepr = TfPyRepr(self.GetPrim());
//        return TfStringPrintf(
//            "UsdGeom.TetraMesh(%s)",
//            primRepr.c_str());
//    }

} // anonymous namespace

//void wrapUsdGeomTetraMesh()
//{
//    typedef UsdGeomTetraMesh This;
//
//    class_<This, bases<UsdGeomPointBased> >
//        cls("TetraMesh");
//
//    cls
//        .def(init<UsdPrim>(arg("prim")))
//        .def(init<UsdSchemaBase const&>(arg("schemaObj")))
//        .def(TfTypePythonClass())
//
//        .def("Get", &This::Get, (arg("stage"), arg("path")))
//        .staticmethod("Get")
//
//        .def("Define", &This::Define, (arg("stage"), arg("path")))
//        .staticmethod("Define")
//
//        .def("GetSchemaAttributeNames",
//            &This::GetSchemaAttributeNames,
//            arg("includeInherited") = true,
//            return_value_policy<TfPySequenceToList>())
//        .staticmethod("GetSchemaAttributeNames")
//
//        .def("_GetStaticTfType", (TfType const& (*)()) TfType::Find<This>,
//            return_value_policy<return_by_value>())
//        .staticmethod("_GetStaticTfType")
//
//        .def(!self)
//
//
//        .def("GetTetraVertexIndicesAttr",
//            &This::GetTetraVertexIndicesAttr)
//        .def("CreateTetraVertexIndicesAttr",
//            &_CreateTetraVertexIndicesAttr,
//            (arg("defaultValue") = object(),
//                arg("writeSparsely") = false))
//
//        .def("__repr__", ::_Repr)
//        ;
//
//    _CustomWrapCode(cls);
//}


namespace {
    //WRAP_CUSTOM{
    //    typedef UsdGeomTetraMesh This;

    //    _class
    //        .def("GetTetraCount", &UsdGeomTetraMesh::GetTetraCount,
    //            arg("timeCode") = UsdTimeCode::Default());
    //    ;
    //}

} // anonymous namespace 
