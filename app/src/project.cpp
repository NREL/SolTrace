
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


#include <stdio.h>
#include <stdlib.h>
#include <stapi.h>
#include <math.h>

#ifndef M_PI
	#define M_PI 3.141592653589793238462643
#endif

#include <wx/tokenzr.h>
#include <wx/arrstr.h>
#include <wx/wxcrt.h>

#include "project.h"

static void read_line(char *buf, int len, FILE *fp )
{
	fgets(buf, len, fp);
	int nch = strlen(buf);
	if (nch > 0 && buf[nch-1] == '\n')
		buf[nch-1] = 0;
	if (nch-1 > 0 && buf[nch-2] == '\r')
		buf[nch-2] = 0;
}

static void ZeroTransform(double ref[3][3], double loc[3][3], double eu[3])
{
	for (int i=0;i<3;i++)
	{
		eu[i]=0;
		for (int j=0;j<3;j++)
			ref[i][j]=loc[i][j]=0;
	}
}

SunShape::SunShape()
{
	ResetToDefaults();
}

bool SunShape::Write(FILE *fp)
{
	if (!fp) return false;

	char cshape = 'g';
	if (Shape==PILLBOX) cshape = 'p';
	if (Shape==USER_DEFINED) cshape = 'd';
	fprintf(fp, "SUN\tPTSRC\t%d\tSHAPE\t%c\tSIGMA\t%lg\tHALFWIDTH\t%lg\n",
		PointSource?1:0, cshape, Sigma, HalfWidth);
	fprintf(fp, "XYZ\t%lg\t%lg\t%lg\tUSELDH\t%d\tLDH\t%lg\t%lg\t%lg\n",
		X, Y, Z, UseLDHSpec?1:0, Latitude, Day, Hour);
	fprintf(fp, "USER SHAPE DATA\t%d\n", (int)UserShapeData.size());
	for (size_t i=0;i<UserShapeData.size();i++)
		fprintf(fp, "%lg\t%lg\n", UserShapeData[i].x, UserShapeData[i].y);

	return true;
}

void SunShape::ResetToDefaults()
{
	Shape = PILLBOX;
	HalfWidth = 4.65;
	PointSource = false;
	X = 0;
	Y = 0;
	Z = 100;
	UseLDHSpec = false;
	Latitude = 39.5;
	Day = 90;
	Hour = 10;
	UserShapeData.clear();	
}

bool SunShape::Read(FILE *fp)
{
	if (!fp) return false;

	char buf[1024];
	int bi = 0, count = 0;
	char cshape = 'g';

	read_line( buf, 1023, fp );

	sscanf(buf, "SUN\tPTSRC\t%d\tSHAPE\t%c\tSIGMA\t%lg\tHALFWIDTH\t%lg",
		&bi, &cshape, &Sigma, &HalfWidth);
	PointSource = (bi!=0);

	switch(tolower(cshape))
	{
	case 'p': Shape = PILLBOX; break;
	case 'd': Shape = USER_DEFINED; break;
	default: Shape = GAUSSIAN; break;
	}

	read_line( buf, 1023, fp );
	sscanf(buf, "XYZ\t%lg\t%lg\t%lg\tUSELDH\t%d\tLDH\t%lg\t%lg\t%lg",
		&X, &Y, &Z, &bi, &Latitude, &Day, &Hour);
	UseLDHSpec = (bi!=0);

	read_line( buf, 1023, fp );
	sscanf(buf, "USER SHAPE DATA\t%d", &count);
	if (count > 0)
	{
		UserShapeData.resize(count);
		for (int i=0;i<count;i++)
		{
			double x, y;
			read_line( buf, 1023, fp );
			sscanf(buf, "%lg\t%lg", &x, &y);
			UserShapeData[i] = PointF(x,y);
		}
	}

	return true;
}

SurfaceOptic::SurfaceOptic()
{
	ErrorDistribution = 'g';
	ApertureStopOrGratingType = 3;
	OpticalSurfaceNumber = 1;
	DiffractionOrder = 4;
	Reflectivity = 0.96;
	Transmissivity = 1.0;
	RMSSlope = 0.95;
	RMSSpecularity = 0.2;
	RefractionIndexReal = 1.1;
	RefractionIndexImag = 1.2;
	GratingCoeffs[0] = 1.1;
	GratingCoeffs[1] = 1.2;
	GratingCoeffs[2] = 1.3;
	GratingCoeffs[3] = 1.4;
	UseReflectivityTable = false;
	UseTransmissivityTable = false;
}

bool SurfaceOptic::Write(FILE *fp)
{
	if (!fp) return false;
	fprintf(fp, 
		"OPTICAL.v2\t%c\t"
		"%d\t%d\t%d\t"
		"%lg\t%lg\t%lg\t%lg\t"
		"%lg\t%lg\t"
		"%lg\t%lg\t%lg\t%lg\t"
		"%d\t%d\t"
		"%d\t%d\n",

		ErrorDistribution,
		ApertureStopOrGratingType, OpticalSurfaceNumber, DiffractionOrder,
		Reflectivity, Transmissivity, RMSSlope, RMSSpecularity,
		RefractionIndexReal, RefractionIndexImag,
		GratingCoeffs[0], GratingCoeffs[1], GratingCoeffs[2], GratingCoeffs[3],
		UseReflectivityTable ? 1 : 0, (int)ReflectivityTable.size(),
		UseTransmissivityTable ? 1 : 0, (int)TransmissivityTable.size()
		);

	if (UseReflectivityTable)
		for (size_t i=0;i<ReflectivityTable.size();i++)
			fprintf(fp, "%lg %lg\n", ReflectivityTable[i].x, ReflectivityTable[i].y );
	if (UseTransmissivityTable)
		for (size_t i = 0; i < TransmissivityTable.size(); i++)
			fprintf(fp, "%lg %lg\n", TransmissivityTable[i].x, TransmissivityTable[i].y);

	return true;
}

bool SurfaceOptic::Read(FILE *fp, bool oldfmt)
{
	if (!fp) return false;
	char buf[1024];

	if (oldfmt)
	{
		// read in old format
		read_line( buf, 1023, fp );
		read_line( buf, 1023, fp );
		read_line( buf, 1023, fp ); OpticalSurfaceNumber = atoi(buf);
		read_line( buf, 1023, fp );
		read_line( buf, 1023, fp ); RefractionIndexReal = atof(buf);
		read_line( buf, 1023, fp ); RefractionIndexImag = atof(buf);
		read_line( buf, 1023, fp );
		read_line( buf, 1023, fp ); ApertureStopOrGratingType = atoi(buf);
		read_line( buf, 1023, fp );
		read_line( buf, 1023, fp ); DiffractionOrder = atoi(buf);
		read_line( buf, 1023, fp );
		read_line( buf, 1023, fp ); GratingCoeffs[0] = atof(buf);
		read_line( buf, 1023, fp ); GratingCoeffs[1] = atof(buf);
		read_line( buf, 1023, fp ); GratingCoeffs[2] = atof(buf);
		read_line( buf, 1023, fp ); GratingCoeffs[3] = atof(buf);
		read_line( buf, 1023, fp );
		read_line( buf, 1023, fp ); Reflectivity = atof(buf);
		read_line( buf, 1023, fp );
		read_line( buf, 1023, fp ); Transmissivity = atof(buf);
		read_line( buf, 1023, fp );
		read_line( buf, 1023, fp ); RMSSlope = atof(buf);
		read_line( buf, 1023, fp );
		read_line( buf, 1023, fp ); RMSSpecularity = atof(buf);
		read_line( buf, 1023, fp );
		read_line( buf, 1023, fp ); if (strlen(buf) > 0) ErrorDistribution = buf[0];
		UseReflectivityTable = false;
		ReflectivityTable.clear();
		/*
		>> not sure this is needed for old format?? mjw
		UseTransmissivityTable= false;
		TransmissivityTable.clear();*/
	}
	else
	{
		read_line(buf, 1023, fp);
		
		// TODO: make sure wxStringTokenize keeps empty parts
		wxArrayString parts = wxStringTokenize( buf, '\t', wxTOKEN_RET_EMPTY_ALL );

		if (parts.size() < 15)
		{
			return false;
		}

		if (parts[1].length() > 0)
			ErrorDistribution = parts[1][0];

		ApertureStopOrGratingType = wxAtoi( parts[2]);
		OpticalSurfaceNumber = wxAtoi( parts[3]);
		DiffractionOrder = wxAtoi( parts[4]);
		Reflectivity = wxAtof( parts[5]);
		Transmissivity = wxAtof( parts[6]);
		RMSSlope = wxAtof( parts[7]);
		RMSSpecularity = wxAtof( parts[8]);
		RefractionIndexReal = wxAtof( parts[9]);
		RefractionIndexImag = wxAtof( parts[10]);
		GratingCoeffs[0] = wxAtof( parts[11]);
		GratingCoeffs[1] = wxAtof( parts[12]);
		GratingCoeffs[2] = wxAtof( parts[13]);
		GratingCoeffs[3] = wxAtof( parts[14]);

		UseReflectivityTable = false;
		ReflectivityTable.clear();
		UseTransmissivityTable = false;
		TransmissivityTable.clear();
		
		int refl_count = 0;
		int trans_count = 0;
		if (parts.size() > 15)
		{
			UseReflectivityTable = (wxAtoi( parts[15] ) > 0);
			refl_count = wxAtoi( parts[16] );
			if (parts.size() > 17)
			{
				UseTransmissivityTable = (wxAtoi(parts[17]) > 0);
				trans_count = wxAtoi(parts[18]);
			}
		}

		if (UseReflectivityTable)
		{
			ReflectivityTable.clear();
			for (int i=0;i< refl_count;i++)
			{
				read_line(buf,1023,fp);
				double x = 0.0, y = 0.0;
				sscanf(buf, "%lg %lg", &x, &y);
				ReflectivityTable.push_back( PointF(x, y) );
			}
		}
		if (UseTransmissivityTable)
		{
			TransmissivityTable.clear();
			for (int i = 0; i < trans_count; i++)
			{
				read_line(buf, 1023, fp);
				double x = 0.0, y = 0.0;
				sscanf(buf, "%lg %lg", &x, &y);
				TransmissivityTable.push_back(PointF(x, y));
			}
		}

	}

	return true;
}

Optical::Optical()
{
	Name = "New Optic";
}

bool Optical::Write(FILE *fp)
{
	if (!fp) return false;

	fprintf(fp, "OPTICAL PAIR\t%s\n", (const char*)Name.c_str());
	Front.Write(fp);
	Back.Write(fp);
	return true;
}

bool Optical::Read(FILE *fp, const wxString &old_name)
{
	if (!fp) return false;
	char buf[1024];
	read_line( buf, 1023, fp );

	if (strncmp( buf, "OPTICAL PAIR", 12) == 0)
	{
		Name = wxString( (const char*)(buf+13) );
		Front.Read(fp);
		Back.Read(fp);
	}
	else
	{
		Name = old_name;
		Front.Read(fp, true);
		read_line( buf, 1023, fp ); // skip next old format line (OPTICAL PROPERTIES)
		Back.Read(fp, true);
	}

	return true;
}


Element::Element()
{
	Enabled = true;
	X=Y=Z=0;
	AX=AY=0; AZ=1;
	ZRot = 0;
	ApertureIndex='r';
	SurfaceIndex='f';
	
	for (int i=0;i<8;i++)
		ApertureParams[i] = SurfaceParams[i] = 0;

	InteractionType = REFLECTION;

	RayHits = 0;
	FinalRayHits = 0;
	ZeroTransform( RRefToLoc, RLocToRef, Euler );
}

void Element::RecomputeTransforms()
{
	double origin[3], aimpoint[3];
	origin[0] = X;
	origin[1] = Y;
	origin[2] = Z;
	aimpoint[0] = AX;
	aimpoint[1] = AY;
	aimpoint[2] = AZ;

	::st_calc_euler_angles( origin, aimpoint, ZRot, Euler);
	::st_calc_transform_matrices( Euler, RRefToLoc, RLocToRef );
}

bool Element::Write(FILE *fp)
{
	if (!fp) return false;

	fprintf(fp, 
		"%d\t"
		"%lg\t%lg\t%lg\t"
		"%lg\t%lg\t%lg\t"
		"%lg\t"
		"%c\t"
		"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t"
		"%c\t"
		"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t"
		"%s\t"
		"%s\t%d\t"
		"%s\n",

		Enabled?1:0,
		X, Y, Z,
		AX, AY, AZ,
		ZRot,
		ApertureIndex,
			ApertureParams[0], ApertureParams[1], ApertureParams[2], ApertureParams[3],
			ApertureParams[4], ApertureParams[5], ApertureParams[6], ApertureParams[7],
		SurfaceIndex,
			SurfaceParams[0], SurfaceParams[1], SurfaceParams[2], SurfaceParams[3],
			SurfaceParams[4], SurfaceParams[5], SurfaceParams[6], SurfaceParams[7],
		(const char*)SurfaceFile.c_str(),
		(const char*)OpticName.c_str(), InteractionType,
		(const char*)Comment.c_str());


	return true;
}

bool Element::Read(FILE *fp)
{
	if (!fp) return false;

	char buf[1024];
	read_line(buf, 1023, fp);

	wxArrayString tok = wxStringTokenize( buf, '\t', wxTOKEN_RET_EMPTY_ALL );
	if (tok.size() < 30)
		return false;
	

	Enabled = wxAtoi(tok[0]) ? true : false ;
	X = wxAtof( tok[1] );
	Y = wxAtof( tok[2] );
	Z = wxAtof( tok[3] );
	AX = wxAtof( tok[4] );
	AY = wxAtof( tok[5] );
	AZ = wxAtof( tok[6] );
	ZRot = wxAtof( tok[7] );
	if (tok[8].length() > 0)
		ApertureIndex = tok[8][0];

	ApertureParams[0] = wxAtof( tok[9] );
	ApertureParams[1] = wxAtof( tok[10] );
	ApertureParams[2] = wxAtof( tok[11] );
	ApertureParams[3] = wxAtof( tok[12] );
	ApertureParams[4] = wxAtof( tok[13] );
	ApertureParams[5] = wxAtof( tok[14] );
	ApertureParams[6] = wxAtof( tok[15] );
	ApertureParams[7] = wxAtof( tok[16] );

	if (tok[17].length() > 0)
		SurfaceIndex = tok[17][0];

	SurfaceParams[0] = wxAtof( tok[18] );
	SurfaceParams[1] = wxAtof( tok[19] );
	SurfaceParams[2] = wxAtof( tok[20] );
	SurfaceParams[3] = wxAtof( tok[21] );
	SurfaceParams[4] = wxAtof( tok[22] );
	SurfaceParams[5] = wxAtof( tok[23] );
	SurfaceParams[6] = wxAtof( tok[24] );
	SurfaceParams[7] = wxAtof( tok[25] );

	SurfaceFile = tok[26];
	
	OpticName = tok[27];

	InteractionType = wxAtoi( tok[28] );
	Comment = tok[29];

	return true;
}

Stage::Stage()
{
	Name = "Untitled";
	X=Y=Z=0;
	AX=AY=0; AZ=1;
	ZRot = 0;
	Virtual = false;
	MultiHit = true;
	TraceThrough = false;

	RayHits = 0;
	ZeroTransform( RRefToLoc, RLocToRef, Euler );
}

Stage::~Stage()
{
	ClearElements();
}

void Stage::RecomputeTransforms()
{
	double origin[3], aimpoint[3];
	origin[0] = X;
	origin[1] = Y;
	origin[2] = Z;
	aimpoint[0] = AX;
	aimpoint[1] = AY;
	aimpoint[2] = AZ;

	::st_calc_euler_angles( origin, aimpoint, ZRot, Euler);
	::st_calc_transform_matrices( Euler, RRefToLoc, RLocToRef );

	for ( size_t i=0;i<ElementList.size();i++)
		ElementList[i]->RecomputeTransforms();
}

void Stage::ClearElements()
{
	for (size_t i=0;i<ElementList.size();i++)
		delete ElementList[i];
	ElementList.clear();
}

bool Stage::Write(FILE *fp)
{
	if (!fp) return false;

	fprintf(fp, "STAGE\tXYZ\t%lg\t%lg\t%lg\tAIM\t%lg\t%lg\t%lg\tZROT\t%lg\tVIRTUAL\t%d\tMULTIHIT\t%d\tELEMENTS\t%d\tTRACETHROUGH\t%d\n",
		X,Y,Z,
		AX,AY,AZ,
		ZRot,
		Virtual?1:0,
		MultiHit?1:0,
		(int)ElementList.size(),
		TraceThrough?1:0);

	fprintf(fp, "%s\n", (const char*)Name.c_str());

	for (size_t i=0;i<ElementList.size();i++)
		ElementList[i]->Write(fp);

	return true;
}

bool Stage::Read(FILE *fp)
{
	if (!fp) return false;

	char buf[1024];
	read_line( buf, 1023, fp );

	int virt=0,multi=1,count=0,tr=0;

	sscanf(buf, "STAGE\tXYZ\t%lg\t%lg\t%lg\tAIM\t%lg\t%lg\t%lg\tZROT\t%lg\tVIRTUAL\t%d\tMULTIHIT\t%d\tELEMENTS\t%d\tTRACETHROUGH\t%d",
		&X, &Y, &Z,
		&AX, &AY, &AZ,
		&ZRot,
		&virt,
		&multi,
		&count,
		&tr );

	read_line( buf, 1023, fp );
	Name = wxString(buf);
	
	Virtual = (virt!=0);
	MultiHit = (multi!=0);
	TraceThrough = (tr!=0);

	bool ok = true;
	ClearElements();
	for (int i=0;i<count;i++)
	{
		Element *elm = new Element;
		if (!elm->Read(fp))
		{
			ok = false;
			delete elm;
		}
		else
			ElementList.push_back(elm);
	}

	return ok;
}

Project::Project()
{
	/* nothing to do */
}

Project::~Project()
{
	ClearOptics();
	ClearStages();
}

void Project::RecomputeTransforms()
{
	for (size_t i=0;i<StageList.size();i++)
		StageList[i]->RecomputeTransforms();
}

Element *Project::GetElement(int nstage, int nelement)
{
	Stage *stage = GetStage(nstage);
	if (stage && nelement >= 0 && nelement < (int)stage->ElementList.size())
		return stage->ElementList[nelement];
	else
		return NULL;
}

Stage *Project::GetStage(int nstage)
{
	if (nstage >= 0 && nstage < (int)StageList.size())
		return StageList[nstage];
	else
		return NULL;
}

void Project::ClearOptics()
{
	for (size_t i=0;i<OpticsList.size();i++)
		delete OpticsList[i];
	OpticsList.clear();
}

void Project::ClearStages()
{
	for (size_t i=0;i<StageList.size();i++)
		delete StageList[i];
	StageList.clear();
}

bool Project::Write(FILE *fp)
{
	if (!fp) return false;

	extern int version_major, version_minor, version_micro;
	fprintf(fp, "# SOLTRACE VERSION %d.%d.%d INPUT FILE\n", version_major, version_minor, version_micro );

	Sun.Write(fp);

	fprintf(fp, "OPTICS LIST COUNT\t%d\n", (int)OpticsList.size());
	for (size_t i=0;i<OpticsList.size();i++)
		OpticsList[i]->Write( fp );

	fprintf(fp, "STAGE LIST COUNT\t%d\n", (int)StageList.size());
	for (size_t i=0;i<StageList.size();i++)
		StageList[i]->Write( fp );

	return true;
}

bool Project::Read(FILE *fp)
{
	if (!fp) return false;

	char buf[1024];

	char c = fgetc(fp);
	if ( c == '#' )
	{
		/* 
		Do a version compatibility check. In the old paradigm, the version was year.month.day of release. 
		Open source reset to version 3.0.0 with symantic versions. Version 3.0.0 is compatible with SolTrace
		major versions 2012-2016 prior.
		*/
		int vmaj = 0, vmin = 0, vmic = 0;
		read_line( buf, 1023, fp ); sscanf( buf, " SOLTRACE VERSION %d.%d.%d INPUT FILE", &vmaj, &vmin, &vmic);

		extern int version_major, version_minor, version_micro;
		unsigned int cur_version = version_major*10000 + version_minor*100 + version_micro;
		unsigned int file_version = vmaj*10000 + vmin*100 + vmic;

		if (vmaj < 2010) //current version is compatible with year versions
		{
			//if file is later than current, or if major version has changed
			if (file_version > cur_version || version_major != vmaj)
				return false;
		}
	}
	else
	{
		ungetc( c, fp );
	}

	if (!Sun.Read(fp)) return false;
	
	bool ok = true;
	int count = 0;

	count = 0;
	read_line( buf, 1023, fp ); sscanf(buf, "OPTICS LIST COUNT\t%d", &count);
	ClearOptics();
	for (int i=0;i<count;i++)
	{
		Optical *opt = new Optical;
		if (!opt->Read( fp ))
		{
			ok = false;
			delete opt;
		}
		else
			OpticsList.push_back( opt );
	}

	count = 0;
	read_line( buf, 1023, fp ); sscanf(buf, "STAGE LIST COUNT\t%d", &count);
	ClearStages();
	for (int i=0;i<count;i++)
	{
		Stage *stage = new Stage;
		if (!stage->Read(fp))
		{
			ok = false;
			delete stage;
		}
		else
			StageList.push_back(stage);
	}

	return ok;
}



RayData::RayData()
{
	Xi = Yi = Zi = 0;
	Xc = Yc = Zc = 0;
	ElementMap = StageMap = RayNumbers = 0;
	Length = 0;
	SunXMin=SunXMax=SunYMin=SunYMax=0.0;
	SunRayCount = 0;

	m_index = 0;
	m_expStage = m_expCoords = 0;
}

RayData::~RayData()
{
	FreeMemory();
}

bool RayData::AllocMemory(size_t npoints)
{
	if (npoints <= 0) return false;

	if (npoints < Length)
	{
		Length = npoints;
		return true;
	}
	try
	{
		FreeMemory();

		Xi = new double[npoints];
		Yi = new double[npoints];
		Zi = new double[npoints];

		Xc = new double[npoints];
		Yc = new double[npoints];
		Zc = new double[npoints];

		ElementMap = new int[npoints];
		StageMap = new int[npoints];
		RayNumbers = new int[npoints];

		Length = npoints;
		return true;
	}
	catch( const std::exception & )
	{
		return false;
	}
}


void RayData::FreeMemory()
{
	if (Xi) delete [] Xi;
	if (Yi) delete [] Yi;
	if (Zi) delete [] Zi;
	if (Xc) delete [] Xc;
	if (Yc) delete [] Yc;
	if (Zc) delete [] Zc;
	if (ElementMap) delete [] ElementMap;
	if (StageMap) delete [] StageMap;
	if (RayNumbers) delete [] RayNumbers;

	Xi=Yi=Zi = 0;
	Xc=Yc=Zc = 0;
	ElementMap = 0;
	StageMap = 0;
	RayNumbers = 0;

	Length = 0;

	SunXMin=SunXMax=SunYMin=SunYMax=0.0;
	SunRayCount = 0;

}

bool RayData::ReadResultsFromContextList(const std::vector<st_context_t> &list)
{
	if ( list.size() == 0 ) return false;

	size_t i, j, npoints = 0;
	for (i=0;i<list.size();i++)
		npoints += ::st_num_intersections( list[i] );

	if (npoints < 0)
		return false;

	if (!AllocMemory( npoints ))
		return false;

	double *xi, *yi, *zi, *xc, *yc, *zc;
	int *em, *sm, *rn;

	xi = Xi;
	yi = Yi;
	zi = Zi;
	xc = Xc;
	yc = Yc;
	zc = Zc;
	em = ElementMap;
	sm = StageMap;
	rn = RayNumbers;

	int max_raynum = 0;

	// initialize sun extents and ray count
	::st_sun_stats( list[0], &SunXMin, &SunXMax, &SunYMin, &SunYMax, &SunRayCount );
	for (i=0;i<list.size();i++)
	{
		if (i>=1)
		{
			// accumulate sun ray counts and extend extents.
			double sxmin, sxmax, symin, symax;
			int srays;
			::st_sun_stats( list[i], &sxmin, &sxmax, &symin, &symax, &srays );
			SunRayCount += srays;

			if (sxmin < SunXMin) SunXMin = sxmin;
			if (sxmax > SunXMax) SunXMax = sxmax;
			if (symin < SunYMin) SunYMin = symin;
			if (symax > SunYMax) SunYMax = symax;
		}
	
		// extract data
		::st_locations( list[i], xi, yi, zi );
		::st_cosines( list[i], xc, yc, zc );
		::st_elementmap( list[i], em );
		::st_stagemap( list[i], sm );
		::st_raynumbers( list[i], rn );

		size_t segment_len = ::st_num_intersections( list[i] );

		// adjust ray number based on current max from previous segment
		for (j=0;j<segment_len;j++)
			rn[j] += max_raynum;

		// find the max in this segment for the next iteration
		for (j=0;j<segment_len;j++)
			if ( rn[j] > max_raynum )
				max_raynum = rn[j];

		// increment storage array pointers
		xi = xi + segment_len;
		yi = yi + segment_len;
		zi = zi + segment_len;

		xc = xc + segment_len;
		yc = yc + segment_len;
		zc = zc + segment_len;

		em = em + segment_len;
		sm = sm + segment_len;
		rn = rn + segment_len;
	}

	return true;
}

bool RayData::WriteDataFile(const wxString &file)
{
	FILE *fp = fopen( file.c_str(), "wb" );
	if (!fp) return false;

	double buf[9];
	buf[0] = SunXMin;
	buf[1] = SunXMax;
	buf[2] = SunYMin;
	buf[3] = SunYMax;
	buf[4] = SunRayCount;
	buf[5] = Length;

	fwrite( (void*)buf, sizeof(double), 6,fp );

	for (size_t i=0;i<Length;i++)
	{
		buf[0] = Xi[i];
		buf[1] = Yi[i];
		buf[2] = Zi[i];
		buf[3] = Xc[i];
		buf[4] = Yc[i];
		buf[5] = Zc[i];
		buf[6] = StageMap[i];
		buf[7] = ElementMap[i];
		buf[8] = RayNumbers[i];

		fwrite( (void*)buf, sizeof(double), 9, fp );
	}

	fclose(fp);

	return true;
}

bool RayData::ReadDataFile(const wxString &file)
{
	FILE *fp = fopen( file.c_str(), "rb" );
	if (!fp) return false;

	double buf[9];
	FreeMemory();

	fread(buf, sizeof(double), 6, fp);
	SunXMin = buf[0];
	SunXMax = buf[1];
	SunYMin = buf[2];
	SunYMax = buf[3];
	SunRayCount = (int) buf[4];
	int NPoints = (int) buf[5];

	if (!AllocMemory(NPoints))
	{
		fclose(fp);
		return false;
	}

	for (int i=0;i<NPoints;i++)
	{
		fread(buf,sizeof(double),9,fp);
		Xi[i] = buf[0];
		Yi[i] = buf[1];
		Zi[i] = buf[2];
		Xc[i] = buf[3];
		Yc[i] = buf[4];
		Zc[i] = buf[5];
		StageMap[i] = (int)buf[6];
		ElementMap[i] = (int)buf[7];
		RayNumbers[i] = (int)buf[8];
	}

	fclose(fp);

	return true;

}

bool RayData::ReadResultsFromContext(st_context_t spcxt)
{
	int npoints = ::st_num_intersections(spcxt);
	if (npoints < 0)
		return false;

	if (!AllocMemory(npoints))
		return false;
	
	::st_locations( spcxt, Xi, Yi, Zi );
	::st_cosines( spcxt, Xc, Yc, Zc );
	::st_elementmap( spcxt, ElementMap );
	::st_stagemap( spcxt, StageMap );
	::st_raynumbers( spcxt, RayNumbers );

	::st_sun_stats( spcxt, &SunXMin, &SunXMax, &SunYMin, &SunYMax, &SunRayCount );
	return true;
}

size_t RayData::PrepareExport( int coords, int istage )
{
	m_index = 0;
	m_expCoords = coords;
	m_expStage = istage;

	size_t nwrite = 0;
	if ( istage == 0 ) nwrite = Length;
	else
	{
		for( size_t i=0;i<Length;i++ )
			if ( istage == StageMap[i] )
				nwrite++;
	}

	return nwrite;
}

bool RayData::GetNextExport( Project &prj, double Pos[3], double Cos[3], 
	int &Elm, int &Stg, int &Ray )
{
	if ( m_index >= Length )
		return false;

	if ( m_expStage > 0 )
	{
		while( m_index < Length && m_expStage != StageMap[m_index] )
			m_index++;

		if ( m_index >= Length )
			return false;
	}

	bool ok = Transform( prj, m_expCoords, m_index, Pos, Cos, Elm, Stg, Ray );
	
	m_index++;

	return ok;
}

bool RayData::Transform( Project &prj, int coord, size_t idx, double Pos[3], double Cos[3], 
	int &Elm, int &Stg, int &Ray )
{

	if ( idx >= Length ) return false;

	double origin[3], posin[3], cosin[3];

	Pos[0] = posin[0] = Xi[idx];
	Pos[1] = posin[1] = Yi[idx];
	Pos[2] = posin[2] = Zi[idx];

	Cos[0] = cosin[0] = Xc[idx];
	Cos[1] = cosin[1] = Yc[idx];
	Cos[2] = cosin[2] = Zc[idx];
		
	if ( coord == RayData::COORD_GLOBAL)
	{
		// global
		if( Stage *stage = prj.GetStage( StageMap[idx]-1 ) )
		{
			origin[0] = stage->X;
			origin[1] = stage->Y;
			origin[2] = stage->Z;

			::st_transform_to_reference( posin, cosin, 
				origin, stage->RLocToRef, 
				Pos, Cos );
		}
		else
		{
			Pos[0] = Pos[1] = Pos[2] = -999;
			Cos[0] = Cos[1] = Cos[2] = -999;
		}


	}
	else if (coord == RayData::COORD_ELEMENT)
	{
		// element

		if ( Element *element = prj.GetElement( StageMap[idx]-1, 
			abs( ElementMap[idx] )-1 ) )
		{
			origin[0] = element->X;
			origin[1] = element->Y;
			origin[2] = element->Z;

			::st_transform_to_local( posin, cosin, 
				origin, element->RRefToLoc,
				Pos, Cos);
		}
		else
		{
			Pos[0] = Pos[1] = Pos[2] = -999;
			Cos[0] = Cos[1] = Cos[2] = -999;
		}
	}

	Elm = ElementMap[idx];
	Stg = StageMap[idx];
	Ray = RayNumbers[idx];

	return true;
}

ElementStatistics::ElementStatistics( Project &prj )
	: m_prj( prj )
{
	ResetData();

}
void ElementStatistics::ResetData()
{
	xValues.clear();
	yValues.clear();
	fluxGrid.clear();

	binszx = binszy = 0.0;
	PowerPerRay = 0.0;
	PeakFlux = PeakFluxUncertainty = 0.0;
	AveFlux = AveFluxUncertainty = 0.0;
	MinFlux = SigmaFlux = Uniformity = 0.0;
	Centroid[0] = Centroid[1] = Centroid[2] = 0.0;
	zScale = Radius = DNI = 0.0;
	NumberOfRays = 0;
}

bool ElementStatistics::Compute(
		int stageIdx, int elementIdx,
		int nbinsx, int nbinsy,
		bool autoscale, bool finalonly,
		double dni,
		double &minx, double &miny,
		double &maxx, double &maxy )
{
	ResetData();

	if (stageIdx < 0 || stageIdx >= m_prj.StageList.size()
		|| elementIdx < 0 || elementIdx >= m_prj.StageList[stageIdx]->ElementList.size())
	{
		return false;
	}

	Element *elm = m_prj.StageList[stageIdx]->ElementList[elementIdx];
	char surf = tolower(elm->SurfaceIndex);
	char aper = tolower(elm->ApertureIndex);
	if (surf != 't' && surf != 'f')
	{
		if (aper !='l' || surf !='s')
			return false;
	}

	RayData *rd = &m_prj.Results;
	if ( autoscale )
	{
		minx = miny = 1e199;
		maxx = maxy = -1e199;
		double Origin[3], PosStage[3], CosStage[3], PosElement[3], CosElement[3];

		// automatically size the min/max x
		for (size_t i=0;i<rd->Length;i++)
		{
			if (abs(rd->ElementMap[i]) == elementIdx+1
				&& rd->StageMap[i] == stageIdx+1)
			{
				// convert from stage to element coordinate system
				PosStage[0] = rd->Xi[i];
				PosStage[1] = rd->Yi[i];
				PosStage[2] = rd->Zi[i];
				CosStage[0] = rd->Xc[i];
				CosStage[1] = rd->Yc[i];
				CosStage[2] = rd->Zc[i];
				Origin[0] = elm->X;
				Origin[1] = elm->Y;
				Origin[2] = elm->Z;

				::st_transform_to_local( PosStage, CosStage,
					Origin, elm->RRefToLoc,
					PosElement, CosElement );

				double x = PosElement[0];
				double y = PosElement[1];

				if (x < minx) minx = x;
				if (x > maxx) maxx = x;
				if (y < miny) miny = y;
				if (y > maxy) maxy = y;
			}
		}
	}

	if (nbinsx <= 1 || nbinsy <= 1
		|| maxx <= minx || maxy <= miny )
	{
		return false;
	}

	double gridszx = maxx - minx;
	double gridszy = maxy - miny;

	if (tolower(elm->SurfaceIndex)=='t' || tolower(elm->SurfaceIndex)=='s')
	{
		Radius = 1.0/elm->SurfaceParams[0]; // CurvOfRev
		if (tolower(elm->SurfaceIndex) == 't') // Full cylinder
		{
			minx = -M_PI * Radius;
			maxx = M_PI * Radius;
			gridszx = 2.0 * M_PI * Radius;
		}
		else  // Partial cylinder
		{
			double angle_min, angle_max;
			angle_min = asin(elm->ApertureParams[0] / Radius);
			angle_max = asin(elm->ApertureParams[1] / Radius);
			minx = Radius * angle_min;
			maxx = Radius * angle_max;
			gridszx = (angle_max - angle_min) * Radius;
		}
	}

	binszx = gridszx / nbinsx;
	binszy = gridszy / nbinsy;

	xValues.resize(nbinsx);
	yValues.resize(nbinsy);
	fluxGrid.resize(nbinsx, nbinsy);

	double zval = 0;
	BinRaysXY(elm,
			  stageIdx,
			  elementIdx,
		finalonly,  // final rays only (i.e. only absorbed ones)
		Radius,
		minx, miny,
		zval, NumberOfRays); // 'out' variables

	double SumFlux, SumFlux2;
	PeakFlux=SumFlux=SumFlux2=0;
	MinFlux = 1e99;

	int NRaysInMinFluxBin=-1, NRaysInPeakFluxBin=-1;
	this->DNI = dni;
	PowerPerRay = dni
		* (rd->SunXMax - rd->SunXMin)
		* (rd->SunYMax - rd->SunYMin)
		/ rd->SunRayCount;


	zScale = PowerPerRay/(binszx*binszy);

	for (size_t r=0;r<fluxGrid.nrows();r++)
	{
		for (size_t c=0;c<fluxGrid.ncols();c++)
		{
			double z = fluxGrid.at(r,c)*zScale;
			SumFlux += z;
			SumFlux2 += z*z;

			if (z>PeakFlux)
			{
				PeakFlux = z;
				NRaysInPeakFluxBin = (int)fluxGrid.at(r,c);
			}

			if (z<MinFlux)
			{
				MinFlux = z;
				NRaysInMinFluxBin = (int)fluxGrid.at(r,c);
			}
		}
	}

	AveFlux = SumFlux/(nbinsx*nbinsy);
	SigmaFlux = sqrt( (nbinsx*nbinsy*SumFlux2-SumFlux*SumFlux)/(nbinsx*nbinsy*nbinsx*nbinsy) );
	Uniformity = SigmaFlux/AveFlux;
	PeakFluxUncertainty = 100/sqrt((double)NRaysInPeakFluxBin);
	AveFluxUncertainty = 100/sqrt((double)rd->Length);

	return true;
}

void ElementStatistics::BinRaysXY( Element *elm,
								  int stageIdx,
								  int elementIdx,
							  bool AbsorbedOnly,
							  double Radius,
							  double xmin,
							  double ymin,
							  double &ZVal,
							  size_t &RayCount)
{
	 /*
  Input -   SurfaceType = 1 Plane
						  2 Cylinder
			Radius = Radius of cylinder if SurfaceType = 2
			NumberOfRays = Total Number of Rays traced from sun
			iNumBinsX( or Y) = Number of bins in X and Y directions
			BinSizeX( or Y)  = Size of rectangular bin in X and Y directions
			XMin and YMin = minimum values for data grid
  Output -  fluxGrid = Data grid containing tallied raycounts within each bin
			ZVal = Z value of plane for planar surfaces or Radius of cylinder for cylindrical surfaces
			raycountf = total number of rays which intersected the plane or the cylinder}

*/
	RayData *rd = &m_prj.Results; // untranslated in stage coordinates
	int i;
	double x, y, z;
	int GridIncrementX, GridIncrementY;
	std::vector<int> ElementIndexArray;

	//PosStage, CosStage, PosElement, CosElement: T1DArray;
	RayCount = 0;
	size_t NotBinned = 0;

	size_t npoints = 0;
	Centroid[0] = Centroid[1] = Centroid[2] = 0.0;

	fluxGrid.fill( 0.0 );

	// pre collect element indexes of interest
	for (i=0;i<(int)rd->Length;i++)
		if (abs(rd->ElementMap[i]) == elementIdx+1 && rd->StageMap[i] == stageIdx+1)
			ElementIndexArray.push_back( i );

	double Origin[3], PosStage[3], CosStage[3], PosElement[3], CosElement[3];

	Origin[0] = elm->X;
	Origin[1] = elm->Y;
	Origin[2] = elm->Z;

	bool iscylinder = (tolower(elm->SurfaceIndex)=='t');
	bool ispartialcylinder = (tolower(elm->SurfaceIndex) == 's' && tolower(elm->ApertureIndex) == 'l');

	for (size_t j=0;j<ElementIndexArray.size();j++)
	{

		if (j==ElementIndexArray.size()-1
			|| rd->RayNumbers[ ElementIndexArray[j] ] != rd->RayNumbers[ ElementIndexArray[j+1] ])
		{
			i = ElementIndexArray[j];

			if (AbsorbedOnly && rd->ElementMap[i] >= 0)
				continue; //if accounting for absorbed rays only, (not all final ray intersections, only absorbed)

			PosStage[0] = rd->Xi[i]; //convert to element coordinate system
			PosStage[1] = rd->Yi[i];
			PosStage[2] = rd->Zi[i];
			CosStage[0] = rd->Xc[i];
			CosStage[1] = rd->Yc[i];
			CosStage[2] = rd->Zc[i];

			::st_transform_to_local( PosStage, CosStage,
									 Origin, elm->RRefToLoc,
									 PosElement, CosElement );

			x = PosElement[0];
			y = PosElement[1];
			z = PosElement[2];

			Centroid[0] += x;
			Centroid[1] += y;
			Centroid[2] += z;
			npoints++;

			if (iscylinder)
			{
				ZVal = Radius;
				if (z<=Radius)
					x = Radius*asin(x/Radius);
				else if (z > Radius)
				{
					if (x < 0) x = -(M_PI*Radius/2.0 + Radius*acos(fabs(x)/Radius));
					if (x >= 0) x = M_PI*Radius/2.0 + Radius*acos(x/Radius);
				}
			}
			else if (ispartialcylinder)
			{
				ZVal = Radius;
				x = Radius * asin(x / Radius);
			}
			else
			{
				ZVal = z; //offload z position of flux surface
			}

			GridIncrementX = -1; //initialize grid increment counters
			GridIncrementY = -1;

			//determine which bin the ray falls into
			while ((xmin + (GridIncrementX+1)*binszx) < x)
				GridIncrementX++;

			while ((ymin + (GridIncrementY+1)*binszy) < y)
				GridIncrementY++;

			if (GridIncrementX >= 0 && GridIncrementX < (int)fluxGrid.nrows()
				&& GridIncrementY >= 0 && GridIncrementY < (int)fluxGrid.ncols() )
			{
				fluxGrid.at(GridIncrementX,GridIncrementY) += 1;//if ray falls inside a bin, increment count for that bin
				RayCount++;  //increment ray intersection counter
			}
			else
			{
				NotBinned++;
			//	qDebug("Not binned: [%d %d],  x=%lg, y=%lg", GridIncrementX, GridIncrementY, x, y);
			}
		}
	}

	//qDebug("BinRaysXY: RayCount=%u NotBinned=%u", RayCount, NotBinned);

	//calculate midpoints of bins
	for (i=0;i<xValues.size();i++)
		xValues[i] = xmin + binszx/2.0 + i*binszx;

	for (i=0;i<yValues.size();i++)
		yValues[yValues.size()-i-1] = ymin + binszy/2.0 + i*binszy;

	if (npoints > 0)
	{
		Centroid[0] /= npoints;
		Centroid[1] /= npoints;
		Centroid[2] /= npoints;
	}
}

