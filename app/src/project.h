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
	std::vector<PointF> ReflectivityTable;
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
