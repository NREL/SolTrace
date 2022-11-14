
/*******************************************************************************************************
*  Copyright 2018 Alliance for Sustainable Energy, LLC
*
*  NOTICE: This software was developed at least in part by Alliance for Sustainable Energy, LLC
*  ("Alliance") under Contract No. DE-AC36-08GO28308 with the U.S. Department of Energy and the U.S.
*  The Government retains for itself and others acting on its behalf a nonexclusive, paid-up,
*  irrevocable worldwide license in the software to reproduce, prepare derivative works, distribute
*  copies to the public, perform publicly and display publicly, and to permit others to do so.
*
*  Redistribution and use in source and binary forms, with or without modification, are permitted
*  provided that the following conditions are met:
*
*  1. Redistributions of source code must retain the above copyright notice, the above government
*  rights notice, this list of conditions and the following disclaimer.
*
*  2. Redistributions in binary form must reproduce the above copyright notice, the above government
*  rights notice, this list of conditions and the following disclaimer in the documentation and/or
*  other materials provided with the distribution.
*
*  3. The entire corresponding source code of any redistribution, with or without modification, by a
*  research entity, including but not limited to any contracting manager/operator of a United States
*  National Laboratory, any institution of higher learning, and any non-profit organization, must be
*  made publicly available under this license for as long as the redistribution is made available by
*  the research entity.
*
*  4. Redistribution of this software, without modification, must refer to the software by the same
*  designation. Redistribution of a modified version of this software (i) may not refer to the modified
*  version by the same designation, or by any confusingly similar designation, and (ii) must refer to
*  the underlying software originally provided by Alliance as "SolTrace". Except to comply with the 
*  foregoing, the term "SolTrace", or any confusingly similar designation may not be used to refer to 
*  any modified version of this software or any modified version of the underlying software originally 
*  provided by Alliance without the prior written consent of Alliance.
*
*  5. The name of the copyright holder, contributors, the United States Government, the United States
*  Department of Energy, or any of their employees may not be used to endorse or promote products
*  derived from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
*  IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
*  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER,
*  CONTRIBUTORS, UNITED STATES GOVERNMENT OR UNITED STATES DEPARTMENT OF ENERGY, NOR ANY OF THEIR
*  EMPLOYEES, BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
*  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
*  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
*  IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
*  THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*******************************************************************************************************/


#ifndef __types_h
#define __types_h

#include <vector>
#include <string>
#include <exception>

#include "stapi.h"
#include "mtrand.h"
#include "hpvm.h"
#include "interpolate.h"
#include "treemesh.h"

#define ACOSM1O180 0.017453292519943295 // acos(-1)/180.0
#ifndef M_PI
	#define M_PI 3.141592653589793238462643
#endif

class nanexcept : public std::exception
{
	std::string m_text;
public:
	nanexcept( const char *text ) : m_text(text) { }
	virtual ~nanexcept() throw() {  }
	virtual const char *what() const throw() { return m_text.c_str(); }
};

class FEDataObj : public st_hash_tree
{
public:
	MatDoub nodes;
};


class TOpticalProperties
{
public:
	TOpticalProperties();
	TOpticalProperties &operator=(const TOpticalProperties &rhs);
	
	char DistributionType;
	int OpticSurfNumber;
	int ApertureStopOrGratingType;
	int DiffractionOrder;
	double Reflectivity;
	double Transmissivity;
	double RMSSlopeError;
	double RMSSpecError;

	double RefractiveIndex[4];
	double AB12[4];	

	bool UseReflectivityTable;
	struct refldat { double angle; double refl; };
	std::vector<refldat> ReflectivityTable;
	bool UseTransmissivityTable;
	struct transdat { double angle; double trans; };
	std::vector<transdat> TransmissivityTable;
};

class TOpticalPropertySet
{
public:
	std::string Name;
	TOpticalProperties Front;
	TOpticalProperties Back;	
};

struct TElement
{
	TElement();
	
	bool Enabled;
	
	/////////// ORIENTATION PARAMETERS ///////////////
	double Origin[3];
	double AimPoint[3];
	double ZRot;
	double Euler[3]; // calculated
	double RRefToLoc[3][3]; // calculated
	double RLocToRef[3][3]; // calculated
	double PosSunCoords[3]; // calculated -- position in sun plane coordinates - mw
	
	/////////// APERTURE PARAMETERS ///////////////
	char ShapeIndex;
	double ParameterA;
	double ParameterB;
	double ParameterC;
	double ParameterD;
	double ParameterE;
	double ParameterF;
	double ParameterG;
	double ParameterH;
	
	double ApertureArea; // calculated
	double ZAperture; // calculated 
	
	/////////// SURFACE PARAMETERS ///////////////
	char SurfaceIndex;
	int SurfaceType; // calculated
	std::string SurfaceFile;
	
	double Kappa;
	double Alpha[5];
	double VertexCurvX;
	double VertexCurvY;
	double AnnularRadius;
	double CrossSectionRadius;
	double ConeHalfAngle;
	double CurvOfRev;
		
	int FitOrder;
	

	// Zernike (*.mon) monomial coeffs
	// (also used for VSHOT Zernike fits)
	HPM2D BCoefficients;

	// Rotationally symmetric polynomial coeffs
	std::vector< double > PolyCoeffs;
	
	// Rotationally symmetric cubic spline
	std::vector< double > CubicSplineXData;
	std::vector< double > CubicSplineYData; 
	std::vector< double > CubicSplineY2Data;   
	double CubicSplineDYDXbc1;
	double CubicSplineDYDXbcN;
	
	// VSHOT file data
	HPM2D VSHOTData;
	double VSHOTRMSSlope;
	double VSHOTRMSScale;
	double VSHOTRadius;
	double VSHOTFocLen;
	double VSHOTTarDis;
	
	// Finite Element data coeffs
	//HPM2D FEData;	
	//GaussMarkov* FEMeshInterp;
	//GaussMarkov FEData;
	FEDataObj FEData;
	
	/////////// OPTICAL PARAMETERS ///////////////
	int InteractionType;
	TOpticalPropertySet *Optics;
	std::string OpticName;

	std::string Comment;	
    int element_number;     //mjw element number in the stage - unique ID in order of addition to element list
};

struct TSun
{
	TSun();
	void Reset();
	
	char ShapeIndex;	
	double Sigma;
	bool PointSource;
	
	std::vector<double> SunShapeAngle;
	std::vector<double> SunShapeIntensity;
	double MaxAngle;
	double MaxIntensity;
	
	
	double Origin[3];
	
	//calculated
	double Euler[3];
	double RRefToLoc[3][3];
	double RLocToRef[3][3];
	
	double MaxRad;
	double Xcm;
	double Ycm;
	double MinXSun;
	double MaxXSun;
	double MinYSun;
	double MaxYSun;
};

class TRayData
{
public:
	TRayData();
	~TRayData();

	struct ray_t
	{
		double pos[3];
		double cos[3];
		int element;
		int stage;
		unsigned int raynum;
	};

	ray_t *Append( double pos[3],
					 double cos[3],
					 int element,
					 int stage,
					 unsigned int raynum );

	bool Overwrite( unsigned int idx,
					double pos[3],
					double cos[3],
					int element,
					int stage,
					unsigned int raynum);

	bool Query( unsigned int idx,
					double pos[3],
					double cos[3],
					int *element,
					int *stage,
					unsigned int *raynum);

	void Merge( TRayData &dest );

	void Clear();

	void Print();

	st_uint_t Count();

	ray_t *Index(st_uint_t i, bool write_access);

private:
	static const unsigned int block_size = 8192;


	struct block_t
	{
		ray_t data[block_size];
		st_uint_t count;
	};

	std::vector<block_t*> m_blockList;
	st_uint_t m_dataCount;
	st_uint_t m_dataCapacity;
};


struct TStage
{
	TStage();
	~TStage();
		
	bool MultiHitsPerRay;
	bool Virtual;
	bool TraceThrough;
	
	double Origin[3];
	double AimPoint[3];
	double ZRot;

	std::vector<TElement*> ElementList;
	
	// calculated
	double Euler[3];
	double RRefToLoc[3][3];
	double RLocToRef[3][3];
	
	TRayData RayData;
};

struct TSystem
{
	TSystem();
	~TSystem();

	void ClearAll();
	
	TSun Sun;
	std::vector<TOpticalPropertySet*> OpticsList;
	std::vector<TStage*> StageList;


	// system simulation context data
	int sim_raycount;
	int sim_raymax;
	bool sim_errors_sunshape;
	bool sim_errors_optical;

	// simulation outputs
	TRayData AllRayData;
	st_uint_t SunRayCount;

	std::vector<std::string> messages;

	void errlog(const char *fmt, ...);
};

#endif
