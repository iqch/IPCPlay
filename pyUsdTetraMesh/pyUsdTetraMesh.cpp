#include "pxr/base/tf/pySafePython.h"
#include "pxr/pxr.h"
#include "pxr/base/tf/pyModule.h"

PXR_NAMESPACE_USING_DIRECTIVE

TF_WRAP_MODULE
{
    TF_WRAP(UsdSchemaExamplesSimple);
    TF_WRAP(UsdSchemaExamplesComplex);
    TF_WRAP(UsdSchemaExamplesParamsAPI);
    TF_WRAP(UsdSchemaExamplesTokens);
}

#include "pxr/pxr.h"
#include "pxr/base/tf/registryManager.h"
#include "pxr/base/tf/scriptModuleLoader.h"
#include "pxr/base/tf/token.h"

#include <vector>

PXR_NAMESPACE_OPEN_SCOPE

TF_REGISTRY_FUNCTION(TfScriptModuleLoader) {
    // List of direct dependencies for this library.
    const std::vector<TfToken> reqs = {
        TfToken("sdf"),
        TfToken("tf"),
        TfToken("usd"),
        TfToken("vt"),
        TfToken("usdGeom"),
    };
    TfScriptModuleLoader::GetInstance().
        RegisterLibrary(TfToken("usdMore"), TfToken("pxr.UsdMore"), reqs);
}

///// INCLUDE


PXR_NAMESPACE_CLOSE_SCOPE

//#include "pxr/usd/usdGeom/tetramesh.h"
#include "pxr/usd/usd/schemaBase.h"

#include "pxr/usd/sdf/primSpec.h"

#include "pxr/usd/usd/pyConversions.h"
#include "pxr/base/tf/pyContainerConversions.h"
#include "pxr/base/tf/pyResultConversions.h"
#include "pxr/base/tf/pyUtils.h"
#include "pxr/base/tf/wrapTypeHelpers.h"

#include <hboost/python.hpp>

#include <string>

using namespace hboost::python;

PXR_NAMESPACE_USING_DIRECTIVE

namespace {

#define WRAP_CUSTOM                                                     \
    template <class Cls> static void _CustomWrapCode(Cls &_class)

// fwd decl.
WRAP_CUSTOM;

        
static UsdAttribute
_CreateTetraVertexIndicesAttr(UsdGeomTetraMesh &self,
                                      object defaultVal, bool writeSparsely)
{
    return self.CreateTetraVertexIndicesAttr(
        UsdPythonToSdfType(defaultVal, SdfValueTypeNames->Int4Array), writeSparsely);
}
        
static std::string
_Repr(const UsdGeomTetraMesh &self)
{
    std::string primRepr = TfPyRepr(self.GetPrim());
    return TfStringPrintf(
        "UsdGeom.TetraMesh(%s)",
        primRepr.c_str());
}

} // anonymous namespace

void wrapUsdGeomTetraMesh()
{
    typedef UsdGeomTetraMesh This;

    class_<This, bases<UsdGeomPointBased> >
        cls("TetraMesh");

    cls
        .def(init<UsdPrim>(arg("prim")))
        .def(init<UsdSchemaBase const&>(arg("schemaObj")))
        .def(TfTypePythonClass())

        .def("Get", &This::Get, (arg("stage"), arg("path")))
        .staticmethod("Get")

        .def("Define", &This::Define, (arg("stage"), arg("path")))
        .staticmethod("Define")

        .def("GetSchemaAttributeNames",
             &This::GetSchemaAttributeNames,
             arg("includeInherited")=true,
             return_value_policy<TfPySequenceToList>())
        .staticmethod("GetSchemaAttributeNames")

        .def("_GetStaticTfType", (TfType const &(*)()) TfType::Find<This>,
             return_value_policy<return_by_value>())
        .staticmethod("_GetStaticTfType")

        .def(!self)

        
        .def("GetTetraVertexIndicesAttr",
             &This::GetTetraVertexIndicesAttr)
        .def("CreateTetraVertexIndicesAttr",
             &_CreateTetraVertexIndicesAttr,
             (arg("defaultValue")=object(),
              arg("writeSparsely")=false))
        
        //.def("GetFaceVertexCountsAttr",
        //     &This::GetFaceVertexCountsAttr)
        //.def("CreateFaceVertexCountsAttr",
        //     &_CreateFaceVertexCountsAttr,
        //     (arg("defaultValue")=object(),
        //      arg("writeSparsely")=false))
        
        //.def("GetSubdivisionSchemeAttr",
        //     &This::GetSubdivisionSchemeAttr)
        //.def("CreateSubdivisionSchemeAttr",
        //     &_CreateSubdivisionSchemeAttr,
        //     (arg("defaultValue")=object(),
        //      arg("writeSparsely")=false))
        
        //.def("GetInterpolateBoundaryAttr",
        //     &This::GetInterpolateBoundaryAttr)
        //.def("CreateInterpolateBoundaryAttr",
        //     &_CreateInterpolateBoundaryAttr,
        //     (arg("defaultValue")=object(),
        //      arg("writeSparsely")=false))
        
        //.def("GetFaceVaryingLinearInterpolationAttr",
        //     &This::GetFaceVaryingLinearInterpolationAttr)
        //.def("CreateFaceVaryingLinearInterpolationAttr",
        //     &_CreateFaceVaryingLinearInterpolationAttr,
        //     (arg("defaultValue")=object(),
        //      arg("writeSparsely")=false))
        
        //.def("GetTriangleSubdivisionRuleAttr",
        //     &This::GetTriangleSubdivisionRuleAttr)
        //.def("CreateTriangleSubdivisionRuleAttr",
        //     &_CreateTriangleSubdivisionRuleAttr,
        //     (arg("defaultValue")=object(),
        //      arg("writeSparsely")=false))
        
        //.def("GetHoleIndicesAttr",
        //     &This::GetHoleIndicesAttr)
        //.def("CreateHoleIndicesAttr",
        //     &_CreateHoleIndicesAttr,
        //     (arg("defaultValue")=object(),
        //      arg("writeSparsely")=false))
        
        //.def("GetCornerIndicesAttr",
        //     &This::GetCornerIndicesAttr)
        //.def("CreateCornerIndicesAttr",
        //     &_CreateCornerIndicesAttr,
        //     (arg("defaultValue")=object(),
        //      arg("writeSparsely")=false))
        
        //.def("GetCornerSharpnessesAttr",
        //     &This::GetCornerSharpnessesAttr)
        //.def("CreateCornerSharpnessesAttr",
        //     &_CreateCornerSharpnessesAttr,
        //     (arg("defaultValue")=object(),
        //      arg("writeSparsely")=false))
        
        //.def("GetCreaseIndicesAttr",
        //     &This::GetCreaseIndicesAttr)
        //.def("CreateCreaseIndicesAttr",
        //     &_CreateCreaseIndicesAttr,
        //     (arg("defaultValue")=object(),
        //      arg("writeSparsely")=false))
        
        //.def("GetCreaseLengthsAttr",
        //     &This::GetCreaseLengthsAttr)
        //.def("CreateCreaseLengthsAttr",
        //     &_CreateCreaseLengthsAttr,
        //     (arg("defaultValue")=object(),
        //      arg("writeSparsely")=false))
        
        //.def("GetCreaseSharpnessesAttr",
        //     &This::GetCreaseSharpnessesAttr)
        //.def("CreateCreaseSharpnessesAttr",
        //     &_CreateCreaseSharpnessesAttr,
        //     (arg("defaultValue")=object(),
        //      arg("writeSparsely")=false))

        .def("__repr__", ::_Repr)
    ;

    _CustomWrapCode(cls);
}

// ===================================================================== //
// Feel free to add custom code below this line, it will be preserved by 
// the code generator.  The entry point for your custom code should look
// minimally like the following:
//
// WRAP_CUSTOM {
//     _class
//         .def("MyCustomMethod", ...)
//     ;
// }
//
// Of course any other ancillary or support code may be provided.
// 
// Just remember to wrap code in the appropriate delimiters:
// 'namespace {', '}'.
//
// ===================================================================== //
// --(BEGIN CUSTOM CODE)--

namespace {

//
//tuple
//_ValidateTopology(const VtIntArray& faceVertexIndices,
//                  const VtIntArray& faceVertexCounts,
//                  size_t numPoints)
//{
//    std::string reason;
//    bool valid = UsdGeomMesh::ValidateTopology(faceVertexIndices,
//                                               faceVertexCounts,
//                                               numPoints, &reason);
//    return boost::python::make_tuple(valid, reason);
//}


WRAP_CUSTOM {
    typedef UsdGeomTetraMesh This;

    _class
        //.def("ValidateTopology", &_ValidateTopology,
        //     (arg("faceVertexIndices"),
        //      arg("faceVertexCounts"),
        //      arg("numPoints")))
        .def("GetTetraCount", &UsdGeomTetraMesh::GetTetraCount,
            arg("timeCode")=UsdTimeCode::Default());
    ;
}

} // anonymous namespace 
