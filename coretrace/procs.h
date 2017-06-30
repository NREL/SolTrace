#ifndef __procs_h
#define __procs_h

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "types.h"
#include "mtrand.h"
#include "stapi.h"

void Intersect( 
			double PosLoc[3], 
			double CosLoc[3],
			TElement *Element,
			double PosXYZ[3], 
			double CosKLM[3],
			double DFXYZ[3],
			double *PathLength,
			int *ErrorFlag );
			
void Surface(
			double PosXYZ[3],
			TElement *Element,
			double *FXYZ,
			double DFXYZ[3],
			int *ErrorFlag );

void QuadricSurfaceClosedForm(
			TElement *Element,
			double PosLoc[3],
			double CosLoc[3],
			double PosXYZ[3],
			double DFXYZ[3],
			double *PathLength,
			int *ErrorFlag);

void TorusClosedForm(
			TElement *Element,
			double PosLoc[3],
			double CosLoc[3],
			double PosXYZ[3],
			double DFXYZ[3],
			double *PathLength,
			int *ErrorFlag);

void SpencerandMurtySurfaceClosedForm(
			TElement *Element,
			double PosLoc[3],
			double CosLoc[3],
			double PosXYZ[3],
			double DFXYZ[3],
			double *PathLength,
			int *ErrorFlag );

void Interaction(
			MTRand &myrng,
			double PosXYZ[3],
			double CosKLM[3],
			double DFXYZ[3],
			int InteractionType,
			TOpticalProperties *Opticl,
			double Wavelength,
			double PosOut[3],
			double CosOut[3],
			int *ErrorFlag );

void GenerateRay(
			MTRand &myrng,
			double PosSunStage[3],
			double Origin[3],
			double RLocToRef[3][3],
			TSun *Sun,
			double PosRayGlobal[3],
			double CosRayGlobal[3],
            double PosRaySun[3]
            );

bool LoadExistingStage0Ray(
            int index,
            std::vector<std::vector< double> > *raydat,
			double PosRayGlobal[3],
            double CosRayGlobal[3],
            st_uint_t &ElementNum,
            st_uint_t &RayNum );

bool LoadExistingStage1Ray(
            int index,
            std::vector<std::vector< double> > *raydat,
			double PosRayGlobal[3],
            double CosRayGlobal[3],
            int &raynum);

bool SunToPrimaryStage(
				TSystem *System,
				TStage *Stage,
				TSun *Sun,
				double PosSunStage[3]);
				
bool AperturePlane(
			TElement *Element);

void Errors(
			MTRand &myrng,
			double CosIn[3],
			int Source,
			TSun *Sun,
			TElement *Element,
			TOpticalProperties *OptProperties,
			double CosOut[3],
			double DFXYZ[3] );


void SurfaceNormalErrors( MTRand &myrng, double CosIn[3],
						 TOpticalProperties *OptProperties,
						 double CosOut[3] )  throw(nanexcept);

void DetermineElementIntersectionNew(
			TElement *Element,
			double PosRayIn[3],
			double CosRayIn[3],
			double PosRayOut[3],
			double CosRayOut[3],
			double DFXYZ[3],
			double *PathLength,
			int *ErrorFlag,
			int *Intercept,
			int *BacksideFlag );

void NewZStartforCubicSplineSurf(
			double CRadius,
			double PosLoc[3],
			double CosLoc[3],
			char AperShapeIndex,
			double *NewZStart,
			double *PLength,
			int *EFlag );

void SurfaceZatXYPair(
			double PosXYZ[3],
			TElement *Element,
			double *FXYZ,
			int *ErrorFlag );
			



void MatrixVectorMult(double M[3][3], double V[3], double MxV[3]);
void MatrixTranspose(double InputMatrix[3][3], int NumRowsCols, double OutputMatrix[3][3]);

double DOT(double A[3], double B[3]);

void TransformToLocal(double PosRef[3], double CosRef[3], double Origin[3], 
		double RRefToLoc[3][3], 
		double PosLoc[3], double CosLoc[3]);

void TransformToReference(double PosLoc[3], double CosLoc[3], double Origin[3], 
		double RLocToRef[3][3], 
		double PosRef[3], double CosRef[3]);

void CalculateTransformMatrices(double Euler[3], double RRefToLoc[3][3], double RLocToRef[3][3]);

void EvalPoly(double ax, double ay, std::vector<double> &Coeffs, int POrder, double *az); //the 0.0's are values for DeltaX and DeltaY; **[need to look at this further]**

void PolySlope( std::vector<double> &Coeffs, int POrder, double ax, double ay, double *dzdx, double *dzdy);

bool splint( std::vector<double> &xa,
			std::vector<double> &ya,
			std::vector<double> &y2a,
			int n,
			double x,
			double *y,
			double *dydx );
                         
void spline( std::vector<double> &x, 
			std::vector<double> &y,
			int n,
			double yp1, double ypn,
			std::vector<double> &y2 );
			
void piksrt(int n, double arr[5] );

void EvalMono(double ax, double ay, HPM2D &B, int order, double DeltaX, double DeltaY, double *az);

void FEInterpolate(double Xray, double Yray, double Delta, double Density, 
			HPM2D &FEData, int NumFEPoints, 
			double *z, double *zx, double *zy);
			
void MonoSlope(HPM2D &B, int order, double sxp, double syp, double *dzdx, double *dzdy);

void VSHOTInterpolateModShepard( double Xray, double Yray, double Density,
			HPM2D &VSHOTData, int NumVSHOTPoints,
			double *zx, double *zy, int *ErrorFlag);

void FEInterpNew(double Xray, double Yray, double Density,
			HPM2D &FEData, int NumFEPoints,
			double *zr);
			
void Root_432(int order, double Coeffs[5][5], double RealRoots[5], double *ImRoot1, double *ImRoot2);

bool InitGeometries(TSystem *sys);
bool TranslateSurfaceParams( TElement *elm, double params[8]);
bool ReadSurfaceFile( const char *file, TElement *elm );

bool TranslateSurfaceParams( TSystem *sys, TElement *elm, double params[8]);
bool ReadSurfaceFile(const char *file, TElement *elm, TSystem *sys);

inline void CopyVec3( double dest[3], const std::vector<double> &src );
inline void CopyVec3( std::vector<double> &dest, double src[3] );
inline void CopyVec3( double dest[3], double src[3] );


int intri(double x1, double y1,
				 double x2, double y2,
				 double x3, double y3,
				 double xt, double yt);

int inquad(double x1, double y1,
				 double x2, double y2,
				 double x3, double y3,
				 double x4, double y4,
				 double xt, double yt);

bool Trace(TSystem *System, unsigned int seed,
		   st_uint_t NumberOfRays, 
		   st_uint_t MaxNumberOfRays,
		   bool IncludeSunShape, 
		   bool IncludeErrors,
           bool AsPowerTower,
		   int (*callback)(st_uint_t ntracedtotal, st_uint_t ntraced, st_uint_t ntotrace, st_uint_t curstage, st_uint_t nstages, void *data),
		   void *cbdata,
           std::vector<std::vector< double > > *stage0data = 0,
           std::vector<std::vector< double > > *stage1in = 0,
           bool save_stage_data = false);

bool DumpSystem(const char *file, TSystem *sys);

//std::vector< std::string > split( const std::string &str, const std::string &delim, bool ret_empty = false, bool ret_delim = false );

#endif
