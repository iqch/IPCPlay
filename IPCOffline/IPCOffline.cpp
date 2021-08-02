//#include <sys/stat.h> // for mkdir
//#include <stdio.h>

#include <fstream>
#include <string>
#include <ctime>

#include <ghc/fs_std.hpp> // filesystem
#include <spdlog/spdlog.h>

//#include <igl/readOBJ.h>
//#ifdef USE_OPENGL
//#include <igl/opengl/glfw/Viewer.h>
//#endif

#include "F:\\Projects\\Repos\\IPC-master\\build\\api.h"

//#include <CLI/CLI.hpp>

//#include "Types.hpp"
//#include "IglUtils.hpp"
//#include "Config.hpp"
#include "Solver.hpp"
//#include "NeoHookeanEnergy.hpp"
//#include "FixedCoRotEnergy.hpp"
//#include "Timer.hpp"
//#include "getRSS.hpp"
//#include "CCDUtils.hpp"

//////
#define BOOST_ALL_NO_LIB

// USD
#include <pxr/base/tf/token.h>
#include <pxr/usd/usd/common.h>
////#include <pxr/usd/usd/treeIterator.h>
//
#include <pxr/base/gf/vec2f.h>
#include <pxr/base/gf/vec3f.h>

#include <pxr/usd/ar/resolver.h>
//
#include <pxr/usd/sdf/layer.h>
#include <pxr/usd/sdf/path.h>
//
#include <pxr/usd/usd/attribute.h>
#include <pxr/usd/usd/timeCode.h>
#include <pxr/usd/usd/stage.h>
#include <pxr/usd/usdGeom/metrics.h>
#include <pxr/usd/usdGeom/gprim.h>
#include <pxr/usd/usdGeom/scope.h>
//
#include <pxr/usd/usdGeom/scope.h>
#include <pxr/usd/usdGeom/pointBased.h>
#include <pxr/usd/usdGeom/xform.h>

#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/subset.h>

#include <pxr/usd/usdGeom/tetramesh.h>
//
using namespace pxr;


void printProgressBar(
    int cur_iter,
    int max_iter,
    int bar_width = 70,
    bool clear_end = false)
{
    float progress = std::max(0.0f, std::min(1.0f, cur_iter / float(max_iter)));
    int pos = bar_width * progress;
    int len = fmt::format("{:d}", max_iter).length();
    fmt::print(
        "{3:2d}%|{0:█>{1}}{0: >{2}}| {4: >{6}d}/{5:d}\r", //
        "", pos, std::max(bar_width - pos, 0), int(progress * 100.0), cur_iter,
        max_iter, len);
    if (progress >= 1) {
        if (clear_end) {
            fmt::print("{0: {1}}", "", bar_width * 2);
        }
        else {
            std::cout << std::endl;
        }
    }
    std::cout.flush();
}



int main(int argc, char* argv[])
{
#if defined(USE_TBB) && defined(TBB_NUM_THREADS)
    tbb::task_scheduler_init init(TBB_NUM_THREADS);
    //tbb::task_scheduler_init init(16);
    //std::cout << " TBB TH : " << TBB_NUM_THREADS << std::endl;
#endif


    spdlog::set_level(static_cast<spdlog::level::level_enum>(0/*args.logLevel*/));
    //if (args.logLevel == 6) 
    //{
    //    std::cout.setstate(std::ios_base::failbit);
    //}

    std::string  input = argv[1];

    IPC::Solver<DIM> solver(input);

    if (!solver.valid)
    {
        spdlog::error("failed to load script/mesh!");
        return -1;
    };   

    // LOAD SUCCEED
    //assert(DIM == 3);
    // Load mesh
//    Eigen::MatrixXd V, UV, N;
//    Eigen::MatrixXi F, FUV, FN, E;
//
//    std::vector<std::vector<int>> borderVerts_primitive;
//    std::vector<int> componentNodeRange, componentSFRange, componentCERange, componentCoDim;
//    std::vector<std::pair<Eigen::Vector3i, Eigen::Vector3d>> componentMaterial, componentLVels, componentAVels;
//    std::vector<std::pair<Eigen::Vector3i, std::array<Eigen::Vector3d, 2>>> componentInitVels;
//    std::vector<IPC::DirichletBC> DirichletBCs;
//    std::vector<IPC::NeumannBC> NeumannBCs;
//    std::vector<std::pair<int, std::string>> meshSeqFolderPath;
//
//    Eigen::MatrixXd newV;
//    Eigen::MatrixXi newF, newSF, newE;
//
//    V.resize(0, 3);
//    F.resize(0, 4);
//
//    IPC::SF.resize(0, 3);
//
//    E.resize(0, 2);
//    componentNodeRange.emplace_back(0);
//    componentSFRange.emplace_back(0);
//    componentCERange.emplace_back(0);
//
//    int DBCI = 0, NBCI = 0;
//    for (int i = 0; i < (int)config.inputShapePaths.size(); ++i) {
//        componentCoDim.emplace_back(3);
//
//        const auto& inputShapePathStr = config.inputShapePaths[i];
//        const fs::path inputShapePath(inputShapePathStr);
//        const auto& translate = config.inputShapeTranslates[i];
//        const auto& rotate = config.inputShapeRotates[i];
//        const auto& scale = config.inputShapeScales[i];
//        if (!inputShapePath.has_extension()) {
//            IPC::IglUtils::readNodeEle(inputShapePathStr, newV, newF, newSF);
//        }
//        else {
//            const std::string meshFileSuffix = inputShapePath.extension().string();
//            const std::string inputShapePathNoSuffix = (inputShapePath.parent_path() / inputShapePath.stem()).string();
//            if (meshFileSuffix == ".msh") {
//                if (!IPC::IglUtils::readTetMesh(inputShapePathStr, newV, newF, newSF)) {
//                    spdlog::error("Unable to read input msh file: {:s}", inputShapePathStr);
//                    exit(1);
//                }
//                newE.resize(0, 2);
//
//                if (config.inputShapeMaterials[i][0] > 0.0 && config.inputShapeMaterials[i][1] > 0.0 && config.inputShapeMaterials[i][2] > 0.0) {
//                    Eigen::Vector3i startToEnd;
//                    startToEnd[0] = i;
//                    startToEnd[1] = F.rows();
//                    startToEnd[2] = F.rows() + newF.rows();
//                    componentMaterial.emplace_back(startToEnd, config.inputShapeMaterials[i]);
//                }
//            }
//            else if (meshFileSuffix == ".ele") {
//                IPC::IglUtils::readNodeEle(inputShapePathNoSuffix, newV, newF, newSF);
//                newE.resize(0, 2);
//
//                if (config.inputShapeMaterials[i][0] > 0.0 && config.inputShapeMaterials[i][1] > 0.0 && config.inputShapeMaterials[i][2] > 0.0) {
//                    Eigen::Vector3i startToEnd;
//                    startToEnd[0] = i;
//                    startToEnd[1] = F.rows();
//                    startToEnd[2] = F.rows() + newF.rows();
//                    componentMaterial.emplace_back(startToEnd, config.inputShapeMaterials[i]);
//                }
//            }
//            else if (meshFileSuffix == ".obj") {
//                // for kinematic object
//                componentCoDim.back() = 2;
//                if (!igl::readOBJ(inputShapePathStr, newV, newSF)) {
//                    spdlog::error("Unable to read input obj file: {:s}", inputShapePathStr);
//                    exit(-1);
//                }
//                newF.resize(0, 4);
//                newE.resize(0, 2);
//            }
//            else if (meshFileSuffix == ".seg") {
//                // for kinematic object
//                componentCoDim.back() = 1;
//                if (!IPC::IglUtils::readSEG(inputShapePathStr, newV, newE)) {
//                    Eigen::MatrixXi tempF;
//                    if (!igl::readOBJ(inputShapePathNoSuffix + ".obj", newV, tempF)) {
//                        spdlog::error("Unable to read input seg or obj file: {:s}", inputShapePathStr);
//                        exit(1);
//                    }
//
//                    std::set<std::pair<int, int>> edgesSet;
//                    for (int sfI = 0; sfI < tempF.rows(); ++sfI) {
//                        auto finder = edgesSet.find(std::pair<int, int>(tempF(sfI, 1), tempF(sfI, 0)));
//                        if (finder == edgesSet.end()) {
//                            edgesSet.insert(std::pair<int, int>(tempF(sfI, 0), tempF(sfI, 1)));
//                        }
//
//                        finder = edgesSet.find(std::pair<int, int>(tempF(sfI, 2), tempF(sfI, 1)));
//                        if (finder == edgesSet.end()) {
//                            edgesSet.insert(std::pair<int, int>(tempF(sfI, 1), tempF(sfI, 2)));
//                        }
//
//                        finder = edgesSet.find(std::pair<int, int>(tempF(sfI, 0), tempF(sfI, 2)));
//                        if (finder == edgesSet.end()) {
//                            edgesSet.insert(std::pair<int, int>(tempF(sfI, 2), tempF(sfI, 0)));
//                        }
//                    }
//
//                    newE.resize(edgesSet.size(), 2);
//                    int newEI = 0;
//                    for (const auto& eI : edgesSet) {
//                        newE(newEI, 0) = eI.first;
//                        newE(newEI, 1) = eI.second;
//                        ++newEI;
//                    }
//                }
//                newF.resize(0, 4);
//                newSF.resize(0, 3);
//            }
//            else if (meshFileSuffix == ".pt") {
//                // for kinematic object
//                componentCoDim.back() = 0;
//                Eigen::MatrixXi temp;
//                if (!igl::readOBJ(inputShapePathStr, newV, temp)) {
//                    if (!igl::readOBJ(inputShapePathNoSuffix + ".obj", newV, temp)) {
//                        spdlog::error("Unable to read input pt or obj file: {:s}", inputShapePathStr);
//                        exit(1);
//                    }
//                }
//                newF.resize(0, 4);
//                newSF.resize(0, 3);
//                newE.resize(0, 2);
//            }
//            else {
//                spdlog::error("Unsupported tet mesh file format: {:s}", meshFileSuffix);
//                exit(1);
//            }
//
//            if (!std::isnan(config.inputShapeLVels[i][0]) && !std::isnan(config.inputShapeLVels[i][1]) && !std::isnan(config.inputShapeLVels[i][2])) {
//                Eigen::Vector3i startToEnd;
//                startToEnd[0] = i;
//                startToEnd[1] = V.rows();
//                startToEnd[2] = V.rows() + newV.rows();
//                componentLVels.emplace_back(startToEnd, config.inputShapeLVels[i]);
//            }
//            if (!std::isnan(config.inputShapeAVels[i][0]) && !std::isnan(config.inputShapeAVels[i][1]) && !std::isnan(config.inputShapeAVels[i][2])) {
//                Eigen::Vector3i startToEnd;
//                startToEnd[0] = i;
//                startToEnd[1] = V.rows();
//                startToEnd[2] = V.rows() + newV.rows();
//                componentAVels.emplace_back(startToEnd, config.inputShapeAVels[i]);
//            }
//            if (!std::isnan(config.inputShapeInitVels[i][0][0]) && !std::isnan(config.inputShapeInitVels[i][0][1]) && !std::isnan(config.inputShapeInitVels[i][0][2]) && !std::isnan(config.inputShapeInitVels[i][1][0]) && !std::isnan(config.inputShapeInitVels[i][1][1]) && !std::isnan(config.inputShapeInitVels[i][1][2])) {
//                Eigen::Vector3i startToEnd;
//                startToEnd[0] = i;
//                startToEnd[1] = V.rows();
//                startToEnd[2] = V.rows() + newV.rows();
//                componentInitVels.emplace_back(startToEnd, config.inputShapeInitVels[i]);
//            }
//
//            while (DBCI < config.inputShapeDBC.size() && config.inputShapeDBC[DBCI].first == i) {
//                // vertex selection
//                std::vector<int> selectedVerts;
//                const auto& inputDBC = config.inputShapeDBC[DBCI].second;
//                IPC::IglUtils::Init_Dirichlet(newV, inputDBC.minBBox, inputDBC.maxBBox, selectedVerts);
//                for (auto& i : selectedVerts) {
//                    i += V.rows();
//                }
//                if (selectedVerts.size()) {
//                    DirichletBCs.emplace_back(selectedVerts, inputDBC.linearVelocity, inputDBC.angularVelocity, inputDBC.timeRange);
//                }
//                ++DBCI;
//            }
//            while (NBCI < config.inputShapeNBC.size() && config.inputShapeNBC[NBCI].first == i) {
//                // vertex selection
//                std::vector<int> selectedVerts;
//                const auto& inputNBC = config.inputShapeNBC[NBCI].second;
//                IPC::IglUtils::Init_Dirichlet(newV, inputNBC.minBBox, inputNBC.maxBBox, selectedVerts);
//                for (auto& i : selectedVerts) {
//                    i += V.rows();
//                }
//                if (selectedVerts.size()) {
//                    NeumannBCs.emplace_back(selectedVerts, inputNBC.force, inputNBC.timeRange);
//                }
//                ++NBCI;
//            }
//        }
//
//        int existVrt = (int)V.rows();
//        for (int i = 0; i < (int)newV.rows(); ++i) {
//            Eigen::Vector3d p = newV.row(i).transpose();
//            newV.row(i) = (rotate * p.cwiseProduct(scale) + translate).transpose();
//        }
//
//        if (newF.rows()) {
//            newF.array() += existVrt;
//        }
//        if (newSF.rows()) {
//            newSF.array() += existVrt;
//        }
//        if (newE.rows()) {
//            newE.array() += existVrt;
//        }
//
//        auto append = [](auto& dst, const auto& src) {
//            int dstRows = dst.rows();
//            int srcRows = src.rows();
//            int srcCols = src.cols();
//            dst.conservativeResize(dstRows + srcRows, Eigen::NoChange);
//            dst.block(dstRows, 0, srcRows, srcCols) = src;
//        };
//        append(V, newV);
//        if (newF.rows()) {
//            append(F, newF);
//        }
//        if (newSF.rows()) {
//            append(IPC::SF, newSF);
//        }
//        if (newE.rows()) {
//            append(E, newE);
//        }
//
//        componentNodeRange.emplace_back(V.rows());
//        componentSFRange.emplace_back(IPC::SF.rows());
//        componentCERange.emplace_back(E.rows());
//
//        IPC::compVAccSize.emplace_back(V.rows());
//        IPC::compFAccSize.emplace_back(F.rows());
//    }
//
//    UV = V.leftCols(DIM);
//    if (config.rotDeg != 0.0) {
//        const Eigen::Matrix3d rotMtr = Eigen::AngleAxis<double>(config.rotDeg / 180.0 * M_PI,
//            config.rotAxis)
//            .toRotationMatrix();
//        Eigen::Matrix<double, 1, DIM> center = (UV.colwise().maxCoeff() - UV.colwise().minCoeff()) / 2.0;
//#ifdef USE_TBB
//        tbb::parallel_for(0, (int)UV.rows(), 1, [&](int vI)
//#else
//        for (int vI = 0; vI < UV.rows(); ++vI)
//#endif
//        {
//            if constexpr (DIM == 3) {
//                UV.row(vI) = (rotMtr * (UV.row(vI) - center).transpose()).transpose() + center;
//            }
//            else {
//                Eigen::Vector3d pos;
//                pos.head(2) = (UV.row(vI) - center).transpose();
//                pos[2] = 0.0;
//                UV.row(vI) = (rotMtr * pos).transpose().head(2) + center;
//            }
//        }
//#ifdef USE_TBB
//        );
//#endif
//    }
//    if (config.size > 0.0) {
//        V *= config.size / (UV.colwise().maxCoeff() - UV.colwise().minCoeff()).maxCoeff();
//        UV *= config.size / (UV.colwise().maxCoeff() - UV.colwise().minCoeff()).maxCoeff();
//        V.leftCols(DIM).rowwise() -= UV.colwise().minCoeff();
//        UV.rowwise() -= UV.colwise().minCoeff();
//    }
//
//#ifdef FIRST_TIME_STEP
//
//    Eigen::RowVector3d maxCoords = V.colwise().maxCoeff();
//    Eigen::RowVector3d minCoords = V.colwise().minCoeff();
//    Eigen::RowVector3d coordsRange = maxCoords - minCoords;
//#ifdef USE_TBB
//    tbb::parallel_for(0, (int)V.rows(), 1, [&](int vI)
//#else
//    for (int vI = 0; vI < V.rows(); ++vI)
//#endif
//    {
//        // nonuniform scale
//        UV(vI, 0) *= 1.0 + 0.1 * (UV(vI, 0) - minCoords[0]) / coordsRange[0];
//
//        double yRatio = (UV(vI, 1) - minCoords[1]) / coordsRange[1];
//        UV(vI, 1) += yRatio * 0.1 * UV(vI, 2);
//
//        UV(vI, 2) *= 1.0 - 0.1 * (UV(vI, 2) - minCoords[2]) / coordsRange[2];
//    }
//#ifdef USE_TBB
//    );
//#endif
//
//#endif // FIRST_TIME_STEP
//
//    // nonuniform scale
//    //                UV.col(0) *= 1.1;
//    //                UV.col(1) *= 1.2;
//    //                UV.col(2) *= 1.3;
//    // shear
//    //                UV.col(0) += 0.1 * UV.col(1);
//
//    IPC::IglUtils::findBorderVerts(V, borderVerts_primitive, config.handleRatio);
//
//    IPC::IglUtils::buildSTri2Tet(F, IPC::SF, IPC::sTri2Tet);
//
//
//    //int  vertAmt_input = V.rows();
//
//    // construct mesh data structure
//    IPC::Mesh<DIM> Source(V, F, IPC::SF, E, UV,
//        componentNodeRange, componentSFRange, componentCERange, componentCoDim,
//        componentMaterial, componentLVels, componentAVels, componentInitVels, DirichletBCs, NeumannBCs,
//        config.inputShapeMeshSeqFolderPath,
//        config.YM, config.PR, config.rho);
//    
//    // primitive test cases
//    Source.borderVerts_primitive = borderVerts_primitive;

//    if (config.ccdMethod == ccd::CCDMethod::TIGHT_INCLUSION)
//    {
//        // Compute a conservative error for the tight inclusion CCD
//        computeTightInclusionError(Source, config.meshCollisionObjects);
//    }
//    else if (config.ccdMethod == ccd::CCDMethod::FLOATING_POINT_ROOT_PARITY)
//    {
//#ifdef USE_FPRP_CCD
//        // shift entire mesh so the CCD will be exact in doubles
//        IPC::invShift = shiftVertices(*temp, config.meshCollisionObjects);
//#else
//        spdlog::error("FPRP CCD is disabled in CMake (IPC_WITH_FPRP=OFF)!");
//        exit(1);
//#endif
//    }

    //{
    //    // for output surface mesh
    //    IPC::isSurfNode.resize(0);
    //    IPC::isSurfNode.resize(Source.V.rows(), false);
    //    for (int tI = 0; tI < IPC::SF.rows(); ++tI) {
    //        IPC::isSurfNode[IPC::SF(tI, 0)] = true;
    //        IPC::isSurfNode[IPC::SF(tI, 1)] = true;
    //        IPC::isSurfNode[IPC::SF(tI, 2)] = true;
    //    }

    //    IPC::tetIndToSurf.resize(0);
    //    IPC::tetIndToSurf.resize(Source.V.rows(), -1);
    //    IPC::surfIndToTet.resize(0);
    //    IPC::surfIndToTet.resize(Source.V.rows(), -1);
    //    int sVI = 0;
    //    for (int vI = 0; vI < IPC::isSurfNode.size(); ++vI)
    //    {
    //        if (IPC::isSurfNode[vI]) 
    //        {
    //            IPC::tetIndToSurf[vI] = sVI;
    //            IPC::surfIndToTet[sVI] = vI;
    //            ++sVI;
    //        }
    //    }

    //    IPC::V_surf.resize(sVI, 3);
    //    IPC::F_surf.resize(IPC::SF.rows(), 3);
    //    for (int tI = 0; tI < IPC::SF.rows(); ++tI)
    //    {
    //        IPC::F_surf(tI, 0) = IPC::tetIndToSurf[IPC::SF(tI, 0)];
    //        IPC::F_surf(tI, 1) = IPC::tetIndToSurf[IPC::SF(tI, 1)];
    //        IPC::F_surf(tI, 2) = IPC::tetIndToSurf[IPC::SF(tI, 2)];
    //    }
    //}

    //{
    //    // for output codimensional segment mesh
    //    IPC::isCENode.resize(0);
    //    IPC::isCENode.resize(Source.V.rows(), false);
    //    for (int ceI = 0; ceI < Source.CE.rows(); ++ceI)
    //    {
    //        IPC::isCENode[Source.CE(ceI, 0)] = true;
    //        IPC::isCENode[Source.CE(ceI, 1)] = true;
    //    }

    //    IPC::tetIndToCE.resize(0);
    //    IPC::tetIndToCE.resize(Source.V.rows(), -1);
    //    IPC::CEIndToTet.resize(0);
    //    IPC::CEIndToTet.resize(Source.V.rows(), -1);
    //    int ceVI = 0;
    //    for (int vI = 0; vI < IPC::isCENode.size(); ++vI) 
    //    {
    //        if (IPC::isCENode[vI]) {
    //            IPC::tetIndToCE[vI] = ceVI;
    //            IPC::CEIndToTet[ceVI] = vI;
    //            ++ceVI;
    //        }
    //    }

    //    IPC::V_CE.resize(ceVI, 3);
    //    IPC::F_CE.resize(Source.CE.rows(), 2);
    //    for (int ceI = 0; ceI < Source.CE.rows(); ++ceI)
    //    {
    //        IPC::F_CE(ceI, 0) = IPC::tetIndToCE[Source.CE(ceI, 0)];
    //        IPC::F_CE(ceI, 1) = IPC::tetIndToCE[Source.CE(ceI, 1)];
    //    };
    //}

    fs::path ip(input);
    auto pp = ip.parent_path();

    auto fn = ip.filename();

    time_t rawTime = std::time(NULL);
    char buf[BUFSIZ];
    std::strftime(buf, sizeof(buf), "_%Y%m%d%H%M%S", std::localtime(&rawTime));
   
    std::string clz(buf);

    std::string stagePath = input + clz + ".usd";

    spdlog::critical("STAGE : {:s}", stagePath);

    const char* sp = stagePath.c_str();
    UsdStageRefPtr stage = UsdStage::CreateInMemory();

    TfToken upAxis(std::string("Y"));
    UsdGeomSetStageUpAxis(stage, upAxis);

    stage->SetStartTimeCode(0);

    for (int coI = 0; coI < solver.config.collisionObjects.size(); ++coI)
    {
        IPC::CollisionObject<DIM>& C = *solver.config.collisionObjects[coI];

        std::string name("/CollisionObject_");
        name += C.name;
        //name += std::to_string(coI);

        SdfPath meshpath(name.c_str());

        UsdGeomMesh mesh = UsdGeomMesh::Define(stage, meshpath);

        UsdAttribute mesh_p_attr(mesh.GetPointsAttr());
        UsdAttribute mesh_fc_attr(mesh.GetFaceVertexCountsAttr());
        UsdAttribute mesh_fi_attr(mesh.GetFaceVertexIndicesAttr());

        UsdAttribute mesh_ext_attr(mesh.GetExtentAttr());

        Eigen::MatrixXd& V_surf = C.V;
        Eigen::MatrixXi& F_surf = C.F;

        int VCNT = V_surf.rows();
        VtVec3fArray mesh_P(VCNT);

        for (int i = 0; i < VCNT; i++)
        {
            mesh_P[i][0] = V_surf.row(i)(0);
            mesh_P[i][1] = V_surf.row(i)(1);
            mesh_P[i][2] = V_surf.row(i)(2);
        };

        int FCNT = F_surf.rows();
        VtIntArray mesh_FI(3 * FCNT);
        VtIntArray mesh_FC(FCNT);

        for (int i = 0; i < FCNT; i++)
        {
            mesh_FC[i] = 3;

            mesh_FI[3 * i + 0] = F_surf.row(i)(0);
            mesh_FI[3 * i + 1] = F_surf.row(i)(1);
            mesh_FI[3 * i + 2] = F_surf.row(i)(2);
        };

        mesh_p_attr.Set(mesh_P, 0);
        mesh_fc_attr.Set(mesh_FC, 0);
        mesh_fi_attr.Set(mesh_FI, 0);

        VtVec3fArray extent(2);
        mesh.ComputeExtent(mesh_P, &extent);
        mesh_ext_attr.Set(extent, 0);
    };

    for(int mcoI = 0; mcoI < solver.config.meshCollisionObjects.size(); ++mcoI)
    {
        IPC::CollisionObject<DIM>& C = *solver.config.meshCollisionObjects[mcoI];

        std::string name("/CollisionMesh_");
        name += C.name;

        //name += std::to_string(mcoI);

        SdfPath meshpath(name.c_str());

        UsdGeomMesh mesh = UsdGeomMesh::Define(stage, meshpath);

        UsdAttribute mesh_p_attr(mesh.GetPointsAttr());
        UsdAttribute mesh_fc_attr(mesh.GetFaceVertexCountsAttr());
        UsdAttribute mesh_fi_attr(mesh.GetFaceVertexIndicesAttr());

        UsdAttribute mesh_ext_attr(mesh.GetExtentAttr());

        Eigen::MatrixXd& V_surf = C.V;
        Eigen::MatrixXi& F_surf = C.F;

        int VCNT = V_surf.rows();
        VtVec3fArray mesh_P(VCNT);

        for (int i = 0; i < VCNT; i++)
        {
            mesh_P[i][0] = V_surf.row(i)(0);
            mesh_P[i][1] = V_surf.row(i)(1);
            mesh_P[i][2] = V_surf.row(i)(2);
        };

        int FCNT = F_surf.rows();
        VtIntArray mesh_FI(3 * FCNT);
        VtIntArray mesh_FC(FCNT);

        for (int i = 0; i < FCNT; i++)
        {
            mesh_FC[i] = 3;

            mesh_FI[3 * i + 0] = F_surf.row(i)(0);
            mesh_FI[3 * i + 1] = F_surf.row(i)(1);
            mesh_FI[3 * i + 2] = F_surf.row(i)(2);
        };

        mesh_p_attr.Set(mesh_P, 0);
        mesh_fc_attr.Set(mesh_FC, 0);
        mesh_fi_attr.Set(mesh_FI, 0);

        VtVec3fArray extent(2);
        mesh.ComputeExtent(mesh_P, &extent);
        mesh_ext_attr.Set(extent, 0);
    };


 

    //IPC::Solver<DIM> Solver(config); // , energyTerms, energyParams);

    // ..SAVE MESH
    SdfPath meshpath("/Geometry");

    UsdGeomMesh mesh = UsdGeomMesh::Define(stage, meshpath);

    UsdAttribute mesh_p_attr(mesh.GetPointsAttr());
    UsdAttribute mesh_fc_attr(mesh.GetFaceVertexCountsAttr());
    UsdAttribute mesh_fi_attr(mesh.GetFaceVertexIndicesAttr());

    //UsdAttribute mesh_ti_attr(mesh.GetTetraVertexIndicesAttr());
    UsdAttribute mesh_ext_attr(mesh.GetExtentAttr());

    Eigen::MatrixXd V_surf;
    Eigen::MatrixXi F_surf;
    solver.GetSurfTopology(V_surf, F_surf);

    int VCNT = V_surf.rows();
    VtVec3fArray mesh_P(VCNT);

    for (int i = 0; i < VCNT; i++)
    {
        mesh_P[i][0] = V_surf.row(i)(0);
        mesh_P[i][1] = V_surf.row(i)(1);
        mesh_P[i][2] = V_surf.row(i)(2);
    };

    int FCNT = F_surf.rows();
    VtIntArray mesh_FI(3* FCNT);
    VtIntArray mesh_FC(FCNT);

    for (int i = 0; i < FCNT; i++)
    {
        mesh_FC[i] = 3;

        mesh_FI[3 * i + 0] = F_surf.row(i)(0);
        mesh_FI[3 * i + 1] = F_surf.row(i)(1);
        mesh_FI[3 * i + 2] = F_surf.row(i)(2);
    };

    mesh_p_attr.Set(mesh_P, 0);
    mesh_fc_attr.Set(mesh_FC, 0);
    mesh_fi_attr.Set(mesh_FI, 0);

    VtVec3fArray extent(2);
    mesh.ComputeExtent(mesh_P, &extent);
    mesh_ext_attr.Set(extent,0);

    // ...FACESETS
    int fci = 0;
    VtIntArray vi;
    for (int i = 0; i < solver.config.inputShapeNames.size(); i++)
    {
        std::string sname = "/Geometry/" + solver.config.inputShapeNames[i];
        SdfPath setpath(sname);
        UsdGeomSubset set = UsdGeomSubset::Define(stage, setpath);

        UsdAttribute sa = set.CreateIndicesAttr();

        int FCNUM = solver.config.inputShapeFaceSets[i];

        vi.resize(FCNUM);

        for (int j = 0; j < FCNUM; j++)
        {
            vi[j] = fci; fci++;
        };

        sa.Set(vi);

        spdlog::info("SET {:s} :: {:d}", sname, FCNUM);
    };

    stage->GetRootLayer()->Export(sp);

    //return 0;

    ///////////////////

    int previousSaveIter = std::numeric_limits<int>::min() / 2;

    double save_dt = 1e-2;
    int itersBetweenSaves = round(save_dt / solver.config.dt);

    bool optimization_on = false;
    int iterNum = 0;
    int converged = 0;
    bool outerLoopFinished = false;

    double secPast = 0.0;
    time_t lastStart_world;

    time(&lastStart_world);

    do
    {
        converged = solver.solve(1);
        iterNum = solver.getIterNum();

        if(spdlog::get_level() >= spdlog::level::level_enum::warn)
        {
            printProgressBar(iterNum, solver.getFrameAmt());
        };

        if (!(iterNum != previousSaveIter && iterNum - previousSaveIter < itersBetweenSaves)
            || (converged == 1))
        {
            previousSaveIter = iterNum;

            solver.GetSurfPoints(V_surf);

            for (int i = 0; i < VCNT; i++)
            {
                mesh_P[i][0] = V_surf.row(i)(0);
                mesh_P[i][1] = V_surf.row(i)(1);
                mesh_P[i][2] = V_surf.row(i)(2);
            };

            mesh_p_attr.Set(mesh_P, iterNum);

            VtVec3fArray extent(2);
            mesh.ComputeExtent(mesh_P, &extent);
            mesh_ext_attr.Set(extent, iterNum);

            stage->SetEndTimeCode(iterNum);

            stage->GetRootLayer()->Export(sp);

            spdlog::critical("SAVED FRAME : {:d}", iterNum);

            if (converged == 1) break;
        };

    }
    while(converged != 1);

    secPast += difftime(time(NULL), lastStart_world);

    spdlog::info("optimization converged, with {:d} inner iterations in {:g}s.", solver.getInnerIterAmt(), secPast);
    spdlog::critical("simulation finished");

    return 0;
}
