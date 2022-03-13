
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
#include <ctype.h>
#include <string.h>
#include <stdarg.h>

#include "types.h"
#include "procs.h"

TOpticalProperties::TOpticalProperties()
{
	for (int i=0;i<4;i++) RefractiveIndex[i] = AB12[i] = 0;
	
	OpticSurfNumber = 1;
	ApertureStopOrGratingType = 0;
	DiffractionOrder = 0;
	Reflectivity = 0;
	Transmissivity = 0;
	RMSSlopeError = 0;
	RMSSpecError = 0;
	DistributionType = 'g';
	UseReflectivityTable = false;
	UseTransmissivityTable = false;
}

TOpticalProperties &TOpticalProperties::operator=(const TOpticalProperties &rhs)
{
	DistributionType = rhs.DistributionType;
	OpticSurfNumber = rhs.OpticSurfNumber;
	ApertureStopOrGratingType = rhs.ApertureStopOrGratingType;
	DiffractionOrder = rhs.DiffractionOrder;
	Reflectivity = rhs.Reflectivity;
	Transmissivity = rhs.Transmissivity;
	RMSSlopeError = rhs.RMSSlopeError;
	RMSSpecError = rhs.RMSSpecError;
	UseReflectivityTable = rhs.UseReflectivityTable;
	UseTransmissivityTable = rhs.UseTransmissivityTable;
	
	
	for (int i=0;i<4;i++)
	{
		RefractiveIndex[i] = rhs.RefractiveIndex[i];
		AB12[i] = rhs.AB12[i];
	}

	ReflectivityTable.resize( rhs.ReflectivityTable.size() );
	for (size_t i=0;i< rhs.ReflectivityTable.size();i++)
	{
		ReflectivityTable[i].angle = rhs.ReflectivityTable[i].angle;
		ReflectivityTable[i].refl = rhs.ReflectivityTable[i].refl;
	}
	TransmissivityTable.resize(rhs.TransmissivityTable.size());
	for (size_t i = 0; i < rhs.TransmissivityTable.size(); i++)
	{
		TransmissivityTable[i].angle = rhs.TransmissivityTable[i].angle;
		TransmissivityTable[i].trans = rhs.TransmissivityTable[i].trans;
	}

	return *this;
}

TElement::TElement()
{
	int i, j;
	for (i=0;i<3;i++) Origin[i] = AimPoint[i] = Euler[i] = PosSunCoords[i] = 0;
	for (i=0;i<3;i++) for (j=0;j<3;j++) RRefToLoc[i][j]=RLocToRef[i][j]=0;
	for (i=0;i<5;i++) Alpha[i] = 0;
	
	Enabled = true;
	ZRot = 0;
	ShapeIndex = ' ';
	ParameterA=ParameterB=ParameterC=ParameterD=0;
	ParameterE=ParameterF=ParameterG=ParameterH=0;
	ApertureArea = 0;
	Kappa = 0;
	VertexCurvX = 0;
	VertexCurvY = 0;
	AnnularRadius = 0;
	CrossSectionRadius = 0;
	ConeHalfAngle = 0;
	CurvOfRev = 0;
	SurfaceIndex = ' ';
	SurfaceType = 0;
	
	FitOrder = 0;
	
	CubicSplineDYDXbc1 = 0;
	CubicSplineDYDXbcN = 0;
	
	VSHOTRMSSlope = 0;
	VSHOTRMSScale = 0;
	VSHOTRadius = 0;
	VSHOTFocLen = 0;
	VSHOTTarDis = 0;
	
	InteractionType = 0;
	
	ZAperture = 0;

	Optics = NULL;	
    element_number = -1;        //mjw nonsense
}


TSun::TSun()
{
	Reset();
}

void TSun::Reset()
{
	
	int i,j;
	for (i=0;i<3;i++) Origin[i]=Euler[i]=0;
	for (i=0;i<3;i++) for (j=0;j<3;j++) RRefToLoc[i][j]=RLocToRef[i][j]=0;
	
	PointSource = false;	
	ShapeIndex = ' ';
	Sigma = 0;
	
	MaxAngle = 0;
	MaxIntensity = 0;
	MaxRad = 0;
	Xcm = 0;
	Ycm = 0;
	MaxXSun = 0;
	MinXSun = 0;
	MaxYSun = 0;
	MinYSun = 0;
}


/*
 // Small program to test TRayData memory block allocation scheme

int main(int argc, char *argv[])
{
	TRayData r1, r2;
	double pos[3] = { 4, 2, 6 };
	double cos[3] = { 5, 1, 7 };

	r1.Print();
	for (int i=0; i<26; i++)
	{
		printf("appending %d\n", i);
		if (!r1.Append( pos, cos, -i, -i/10, i+1 ))
			break;

		if (i < 13)
		{
			if (!r2.Append( cos, pos, -2*i, -2*i/10, i+1000 ) )
				break;
		}
	}


	r1.Print();
	r2.Print();

	pos[0] = 44; pos[1] = 11, pos[2] = 22;
	cos[0] = 0.1; cos[1] = 0.2; cos[2] = 0.15;
	r1.Overwrite( 25, pos, cos, 2314, 1255, 491057 );
	r1.Overwrite( 13, pos, cos, 214, 55, 49 );

	printf("\nMERGING...\n\n");
	r1.Merge(r2);

	r1.Print();
	r2.Print();

	return 0;
}
*/




TRayData::TRayData()
{
	m_dataCount = 0;
	m_dataCapacity = 0;
}

TRayData::~TRayData()
{
	Clear();
}

TRayData::ray_t *TRayData::Append( double pos[3],
				 double cos[3],
				 int element,
				 int stage,
				 unsigned int raynum )
{
	if (m_dataCount == m_dataCapacity)
	{
		// need to allocate more blocks
		block_t *b = new block_t;
		b->count = 0;
		m_blockList.push_back( b );
		m_dataCapacity += block_size;
	}

	ray_t *r = Index( m_dataCount, true );
	if( r != 0 )
	{
		::memcpy( &r->pos, pos, sizeof(double)*3 );
		::memcpy( &r->cos, cos, sizeof(double)*3 );
		r->element = element;
		r->stage = stage;
		r->raynum = raynum;
		m_dataCount++;
		return r;
	}
	else return 0;
}

bool TRayData::Overwrite( unsigned int idx,
				double pos[3],
				double cos[3],
				int element,
				int stage,
				unsigned int raynum)
{
	ray_t *r = Index( idx, true );
	if ( r != 0 )
	{
		::memcpy( r->pos, pos, sizeof(double)*3 );
		::memcpy( r->cos, cos, sizeof(double)*3 );
		r->element = element;
		r->stage = stage;
		r->raynum = raynum;
		return true;
	}
	else
		return false;
}

bool TRayData::Query( unsigned int idx,
				double pos[3],
				double cos[3],
				int *element,
				int *stage,
				unsigned int *raynum)
{
	ray_t *r = Index( idx, false );
	if ( r != 0 )
	{
		if (pos!=0) ::memcpy( pos, r->pos, sizeof(double)*3 );
		if (cos!=0) ::memcpy( cos, r->cos, sizeof(double)*3 );
		if (element) *element = r->element;
		if (stage) *stage = r->stage;
		if (raynum) *raynum = r->raynum;
		return true;
	}
	else
		return false;

}

void TRayData::Merge( TRayData &src )
{
	std::vector<block_t*> list, partial_blocks;
	size_t i;

	list.reserve( m_blockList.size() + src.m_blockList.size() );

	for (i=0;i<m_blockList.size();i++)
	{
		if (m_blockList[i]->count == block_size)
			list.push_back( m_blockList[i] );
		else
			partial_blocks.push_back( m_blockList[i] );
	}

	for (i=0;i<src.m_blockList.size();i++)
	{
		if (src.m_blockList[i]->count == block_size)
			list.push_back( src.m_blockList[i] );
		else
			partial_blocks.push_back( src.m_blockList[i] );
	}

	src.m_blockList.clear();
	src.m_dataCount = 0;
	src.m_dataCapacity = 0;

	m_blockList = list;
	m_dataCapacity = m_dataCount = m_blockList.size() * block_size;

	// append all the data in the partial blocks

	for (i=0;i<partial_blocks.size();i++)
	{
		block_t *b = partial_blocks[i];
		for (size_t j=0;j<b->count;j++)
		{
			ray_t &r = b->data[j];
			Append( r.pos, r.cos, r.element, r.stage, r.raynum );
		}

		delete b;
	}
	partial_blocks.clear();
}

void TRayData::Clear()
{
	for (size_t i=0;i<m_blockList.size();i++)
		delete m_blockList[i];
	m_blockList.clear();
	m_dataCount = 0;
	m_dataCapacity = 0;
}

st_uint_t TRayData::Count()
{
	return m_dataCount;
}

TRayData::ray_t *TRayData::Index(st_uint_t i, bool write_access)
{
	if (i >= m_dataCapacity)
		return 0;

	size_t block_num = i / block_size;
	size_t block_idx = i % block_size;

	if (block_num >= m_blockList.size()
		 || block_idx >= block_size )
		return 0;

	// update block.count to highest accessed index
	block_t *b = m_blockList[block_num];

	if (write_access && block_idx >= b->count)
		b->count = block_idx+1;

	if (!write_access && block_idx >= b->count)
		return 0;

	return &(b->data[block_idx]);
}

void TRayData::Print()
{
	printf("[ blocks: %d count: %u capacity: %u ]\n",
		m_blockList.size(),
		(unsigned int)m_dataCount,
		(unsigned int)m_dataCapacity );

	size_t n = Count();
	for (size_t i=0;i<n;i++)
	{
		double pos[3],cos[3];
		int elm, stage;
		unsigned int ray;
		if (Query(i, pos, cos, &elm, &stage, &ray))
		{
			printf("   [%u] = { [%lg,%lg,%lg][%lg,%lg,%lg] %d %d %u }\n", i,
				pos[0], pos[1], pos[2],
				cos[0], cos[1], cos[2],
				elm, stage, ray);
		}
	}

	printf("\n");
}

TStage::TStage()
{
	st_uint_t i,j;
	for (i=0;i<3;i++) Origin[i]=AimPoint[i]=Euler[i]=0;
	for (i=0;i<3;i++) for (j=0;j<3;j++) RRefToLoc[i][j]=RLocToRef[i][j]=0;
	
	ZRot = 0;
	MultiHitsPerRay = true;
	Virtual = false;
	TraceThrough = false;
}

TStage::~TStage()
{
	for (st_uint_t i=0;i<ElementList.size();i++)
		delete ElementList[i];
	ElementList.clear();
}

TSystem::TSystem()
{
	SunRayCount = 0;

	sim_raycount=1000;
	sim_raymax=100000;
	sim_errors_sunshape=true;
	sim_errors_optical=true;
}

TSystem::~TSystem()
{
	for (st_uint_t i=0;i<StageList.size();i++)
		delete StageList[i];
	StageList.clear();
}

void TSystem::ClearAll()
{
	for (st_uint_t i=0;i<OpticsList.size();i++)
		delete OpticsList[i];
	OpticsList.clear();

	for (st_uint_t i=0;i<StageList.size();i++)
		delete StageList[i];
	StageList.clear();
}

void TSystem::errlog(const char *fmt, ...)
{
	static char buf[513];
	va_list arglist;
	va_start( arglist, fmt );
#ifdef WIN32
	_vsnprintf(buf,512,fmt,arglist);
#else
	vsnprintf(buf,512,fmt,arglist);
#endif
	va_end( arglist );	
	messages.push_back(buf);
}

