#ifndef __types_h
#define __types_h

#include <vector>
#include <string>
#include <exception>

#include "stapi.h"
#include "mtrand.h"
#include "hpvm.h"

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
	HPM2D FEData;	
	
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
