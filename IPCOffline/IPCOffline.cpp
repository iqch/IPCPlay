////////////////////////////////////////
//
// IPC Offline simulator v0.9
// Based on https://ipc-sim.github.io/ sources
// minimal-interface playing utility with USD export
// 
////////////////////////////////////////

// CRT
#include <string>
#include <ctime>

// THIRD PARTIES
#include <ghc/fs_std.hpp>
#include <spdlog/spdlog.h>

#define IPC_API __declspec(dllimport)

#include "Solver.hpp"

#define BOOST_ALL_NO_LIB

// USD
#include <pxr/base/tf/token.h>
#include <pxr/usd/usd/common.h>
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

Eigen::Matrix<double, 3, 1> FORCE;

Eigen::Matrix<double, 3, 1> antigravity(const Eigen::Matrix<double, 3, 1>&)
{
    return FORCE;
};

int main(int argc, char* argv[])
{

    FORCE.Zero();
    //FORCE[1] = -9.80665;

#if defined(USE_TBB) && defined(TBB_NUM_THREADS)
    tbb::task_scheduler_init init(TBB_NUM_THREADS);
#endif

    spdlog::set_level(static_cast<spdlog::level::level_enum>(0));

    std::string  input = argv[1];

    IPC::Solver solver(input, antigravity);

    if (!solver.valid)
    {
        spdlog::error("failed to load scene!");
        return -1;
    };   

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
        IPC::CollisionObject& C = *solver.config.collisionObjects[coI];

        std::string name("/CollisionObject_");
        name += C.name;

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
        IPC::CollisionObject& C = *solver.config.meshCollisionObjects[mcoI];

        std::string name("/CollisionMesh_");
        name += C.name;

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

    SdfPath meshpath("/Geometry");

    UsdGeomMesh mesh = UsdGeomMesh::Define(stage, meshpath);

    UsdAttribute mesh_p_attr(mesh.GetPointsAttr());
    UsdAttribute mesh_fc_attr(mesh.GetFaceVertexCountsAttr());
    UsdAttribute mesh_fi_attr(mesh.GetFaceVertexIndicesAttr());

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

    // FACESETS
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

    // ...TETRAMESHES

    // INITIAL EXPORT
    stage->GetRootLayer()->Export(sp);

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
