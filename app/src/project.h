
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


#ifndef __sysdata_h
#define __sysdata_h

#include <exception>
#include <vector>

#include <wx/string.h>

#include <hpvm.h>
#include <stapi.h>

struct PointF
{
	PointF() :x(0),y(0) {}
	PointF(double _x, double _y) :x(_x),y(_y) {}
	double x, y;
};

class SunShape
{
public:
	SunShape();

	void ResetToDefaults();
	bool Write(FILE *fp);
	bool Read(FILE *fp);

	enum{ GAUSSIAN, PILLBOX, USER_DEFINED };

	bool PointSource;
	int Shape;
	double Sigma;
	double HalfWidth;
	double X,Y,Z;
	bool UseLDHSpec;
	double Latitude, Day, Hour;
	std::vector<PointF> UserShapeData;
};

class SurfaceOptic
{
public:
	SurfaceOptic();
	bool Write(FILE *fp);
	bool Read(FILE *fp, bool oldfmt=false);

	char ErrorDistribution;
	int ApertureStopOrGratingType;
	int OpticalSurfaceNumber;
	int DiffractionOrder;
	double Reflectivity;
	double Transmissivity;
	double RMSSlope;
	double RMSSpecularity;
	double RefractionIndexReal;
	double RefractionIndexImag;
	double GratingCoeffs[4];
	bool UseReflectivityTable;
	bool UseTransmissivityTable;
	std::vector<PointF> ReflectivityTable;
	std::vector<PointF> TransmissivityTable;
};

class Optical
{
public:
	Optical();
	bool Write(FILE *fp);
	bool Read(FILE *fp, const wxString &old_name=wxEmptyString);

	wxString Name;
	SurfaceOptic Front;
	SurfaceOptic Back;
};

class Element
{
public:
	Element();
	bool Write(FILE *fp);
	bool Read(FILE *fp);

	bool Enabled;
	
	double X,Y,Z;
	double AX,AY,AZ;
	double ZRot;

	char ApertureIndex;
	double ApertureParams[8];

	char SurfaceIndex;
	double SurfaceParams[8];
	wxString SurfaceFile;

	wxString OpticName;

	enum { REFRACTION=1, REFLECTION=2 };
	int InteractionType;

	wxString Comment;

	
	// initialized when system is simulated
	void RecomputeTransforms();

	int RayHits;
        int FinalRayHits;
	double Euler[3];
	double RRefToLoc[3][3];
	double RLocToRef[3][3];
};

class Stage
{
public:
	Stage();
	~Stage();

	void ClearElements();
	bool Write(FILE *fp);
	bool Read(FILE *fp);

	double X,Y,Z; // coordinates
	double AX,AY,AZ; // aim point
	double ZRot; // Z rotation

	wxString Name;

	bool Virtual;
	bool MultiHit;
	bool TraceThrough;

	std::vector<Element*> ElementList;

	// initialized when system is simulated
	void RecomputeTransforms();

	int RayHits;
	double Euler[3];
	double RRefToLoc[3][3];
	double RLocToRef[3][3];
};

class Project;

class RayData
{
public:
	RayData();
	~RayData();

	enum { COORD_GLOBAL, COORD_STAGE, COORD_ELEMENT };

	bool AllocMemory(size_t nintersect);
	void FreeMemory();

	bool ReadResultsFromContext(st_context_t spcxt);
	bool ReadResultsFromContextList(const std::vector<st_context_t> &list);
		
	size_t PrepareExport( int coords, int istage = 0 );
	bool GetNextExport( Project &prj, double Pos[3], double Cos[3], 
		int &Elm, int &Stg, int &Ray );	
	bool Transform( Project &prj, int coords, size_t idx, double Pos[3], double Cos[3], 
		int &Elm, int &Stg, int &Ray );
		
	bool WriteDataFile( const wxString &file );
	bool ReadDataFile( const wxString &file );
	
	size_t Length;
	double *Xi, *Yi, *Zi;
	double *Xc, *Yc, *Zc;
	int *ElementMap, *RayNumbers, *StageMap;
	double SunXMin, SunXMax, SunYMin, SunYMax;
	int SunRayCount;

private:
	int m_expStage, m_expCoords;
	size_t m_index;

};

class Project
{
public:
	Project();
	~Project();

	void ClearOptics();
	void ClearStages();
	void RecomputeTransforms();
	Element *GetElement(int nstage, int nelement);
	Stage *GetStage(int nstage);
	bool Write(FILE *fp);
	bool Read(FILE *fp);

	SunShape Sun;
	std::vector<Optical*> OpticsList;
	std::vector<Stage*> StageList;

	RayData Results;

};

class ElementStatistics
{
public:
	ElementStatistics( Project &prj );

	bool Compute(
			int stageIdx, int elementIdx,
			int nbinsx, int nbinsy,
			bool autoscale, bool finalonly,
			double dni,
			double &minx, double &miny,
			double &maxx, double &maxy );


	HPM2D fluxGrid;
	std::vector<double> xValues, yValues;
	double binszx, binszy;
	double PowerPerRay;
	double PeakFlux, PeakFluxUncertainty;
	double AveFlux, AveFluxUncertainty;
	double MinFlux, SigmaFlux, Uniformity;
	double Centroid[3];
	double zScale;
	double Radius;
	double DNI;
	size_t NumberOfRays;

private:
	void ResetData();
	void BinRaysXY( Element *elm,
					int stageIdx,
					int elementIdx,
					bool AbsorbedOnly,
					double Radius,
					double xmin,
					double ymin,
					double &ZVal,
					size_t &RayCount);
	Project &m_prj;
};


#endif
