#ifndef SOLTRACE_TRACE_H
#define SOLTRACE_TRACE_H

#include <cstdint>

#include <simulation_data.hpp>
#include <simulation_runner.hpp>

#include "mtrand.hpp"
#include "native_runner_types.hpp"


namespace SolTrace::NativeRunner
{

    class GlobalRay_refactored
    {
    public:
        GlobalRay_refactored() // : active(true)
        {
            Num = 0;
            for (int i = 0; i < 3; i++)
                Pos[i] = Cos[i] = 0.0;
        }

        double Pos[3];
        double Cos[3];
        uint_fast64_t Num;
        // bool active;
    };

    // void SpencerandMurtySurfaceClosedForm(TElement *Element,
    //                                       double PosLoc[3],
    //                                       double CosLoc[3],
    //                                       double PosXYZ[3],
    //                                       double DFXYZ[3],
    //                                       double *PathLength,
    //                                       int *ErrorFlag);

    // bool LoadExistingStage0Ray(int index,
    //                            std::vector<std::vector<double>> *raydat,
    //                            double PosRayGlobal[3],
    //                            double CosRayGlobal[3],
    //                            st_uint_t &ElementNum,
    //                            st_uint_t &RayNum);

    // bool LoadExistingStage1Ray(int index,
    //                            std::vector<std::vector<double>> *raydat,
    //                            double PosRayGlobal[3],
    //                            double CosRayGlobal[3],
    //                            int &raynum);

    // bool AperturePlane(TElement *Element);

    // void NewZStartforCubicSplineSurf(double CRadius,
    //                                  double PosLoc[3],
    //                                  double CosLoc[3],
    //                                  char AperShapeIndex,
    //                                  double *NewZStart,
    //                                  double *PLength,
    //                                  int *EFlag);

    // // the 0.0's are values for DeltaX and DeltaY;
    // **[need to look at this further]**
    // void EvalPoly(double ax,
    //               double ay,
    //               std::vector<double> &Coeffs,
    //               int POrder,
    //               double *az);

    // void PolySlope(std::vector<double> &Coeffs,
    //                int POrder,
    //                double ax,
    //                double ay,
    //                double *dzdx,
    //                double *dzdy);

    // void spline(std::vector<double> &x,
    //             std::vector<double> &y,
    //             int n,
    //             double yp1, double ypn,
    //             std::vector<double> &y2);

    // void EvalMono(double ax,
    //               double ay,
    //               HPM2D &B,
    //               int order,
    //               double DeltaX,
    //               double DeltaY,
    //               double *az);

    // void FEInterpolate(double Xray, double Yray, double Delta, double Density,
    //       HPM2D &FEData, int NumFEPoints,
    //       double *z, double *zx, double *zy);

    // void MonoSlope(HPM2D &B,
    //                int order,
    //                double sxp,
    //                double syp,
    //                double *dzdx,
    //                double *dzdy);

    // void VSHOTInterpolateModShepard(double Xray,
    //                                 double Yray, double Density,
    //                                 HPM2D &VSHOTData,
    //                                 int NumVSHOTPoints,
    //                                 double *zx,
    //                                 double *zy,
    //                                 int *ErrorFlag);

    // void FEInterpNew(double Xray,
    //                  double Yray,
    //                  double Density,
    //                  HPM2D &FEData,
    //                  int NumFEPoints,
    //                  double *zr);
    // void FEInterpGM(double Xray, double Yray, GaussMarkov *gm, double *zr);
    // void FEInterpKD(double Xray,
    //                 double Yray,
    //                 FEDataObj *kd,
    //                 double step,
    //                 double *zr,
    //                 double *dzrdx,
    //                 double *dzrdy);

    // bool InitGeometries(TSystem *sys);
    // bool TranslateSurfaceParams( TElement *elm, double params[8]);
    // bool ReadSurfaceFile( const char *file, TElement *elm );

    // bool TranslateSurfaceParams(TSystem *sys, TElement *elm, double params[8]);
    // bool ReadSurfaceFile(const char *file, TElement *elm, TSystem *sys);

    SolTrace::Runner::RunnerStatus trace_native(
        TSystem *System,
        unsigned int seed,
        uint_fast64_t NumberOfRays,
        uint_fast64_t MaxNumberOfRays,
        bool IncludeSunShape,
        bool IncludeErrors,
        bool AsPowerTower);

} // namespace SolTrace::NativeRunner

#endif
