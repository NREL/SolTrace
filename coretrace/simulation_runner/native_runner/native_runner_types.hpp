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

#ifndef SOLTRACE_NATIVE_RUNNER_TYPES_H
#define SOLTRACE_NATIVE_RUNNER_TYPES_H

#include <cstdint>
#include <exception>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "constants.hpp"
#include "element.hpp"
#include "optical_properties.hpp"
#include "ray_source.hpp"
#include "simulation_result.hpp"
#include "surface_intersection_calculator.hpp"

namespace SolTrace::NativeRunner
{

	// #define ACOSM1O180 0.017453292519943295 // acos(-1)/180.0
	// #endif

	// class nanexcept : public std::exception
	// {
	// 	std::string m_text;

	// public:
	// 	nanexcept(const char *text) : m_text(text) {}
	// 	virtual ~nanexcept() throw() {}
	// 	virtual const char *what() const throw() { return m_text.c_str(); }
	// };

	// class FEDataObj : public st_hash_tree
	// {
	// public:
	// 	MatDoub nodes;
	// };

	// class TOpticalProperties
	// {
	// public:
	// 	TOpticalProperties();
	// 	TOpticalProperties &operator=(const TOpticalProperties &rhs);

	// 	char DistributionType;
	// 	// int OpticSurfNumber;
	// 	// int ApertureStopOrGratingType;
	// 	// int DiffractionOrder;
	// 	double Reflectivity;
	// 	double Transmissivity;
	// 	double RMSSlopeError;
	// 	double RMSSpecError;

	// 	double RefractiveIndex[4];
	// 	double AB12[4];

	// 	bool UseReflectivityTable;
	// 	struct refldat
	// 	{
	// 		double angle;
	// 		double refl;
	// 	};
	// 	std::vector<refldat> ReflectivityTable;
	// 	bool UseTransmissivityTable;
	// 	struct transdat
	// 	{
	// 		double angle;
	// 		double trans;
	// 	};
	// 	std::vector<transdat> TransmissivityTable;
	// };

	// Struct for storing runner side only element parameters
	struct ElementParameters
	{
		ElementParameters();
		~ElementParameters();
		// Newton's method controls
		double newton_tolerance;
		uint_fast64_t newton_max_iters;
	};

	class TOpticalPropertySet
	{
	public:
		std::string Name;
		// TOpticalProperties Front;
		// TOpticalProperties Back;
		SolTrace::Data::OpticalProperties Front;
		SolTrace::Data::OpticalProperties Back;
	};

	struct TElement
	{
		TElement();
		~TElement();

		// bool Enabled;

		/////////// ORIENTATION PARAMETERS ///////////////
		double Origin[3];
		double AimPoint[3];
		double ZRot;
		double RRefToLoc[3][3];
		double RLocToRef[3][3];
		double PosSunCoords[3]; // calculated -- position in sun plane coordinates - mw

		/////////// APERTURE PARAMETERS //////////////
		double ZAperture; // calculated
		SolTrace::Data::aperture_ptr aperture;

		/////////// SURFACE PARAMETERS ///////////////
		calculator_ptr icalc;

		// double Kappa;
		// double Alpha[5];
		// double VertexCurvX;
		// double VertexCurvY;
		// double AnnularRadius;
		// double CrossSectionRadius;
		// double ConeHalfAngle;
		// double CurvOfRev;

		/////////// OPTICAL PARAMETERS ///////////////
		TOpticalPropertySet Optics;

		std::string Comment;
		// mjw element number in the stage - unique ID in order
		// of addition to element list
		int_fast64_t element_number;
		SolTrace::Data::element_id sim_data_id;
	};

	using telement_ptr = typename std::shared_ptr<TElement>;
	telement_ptr make_telement(SolTrace::Data::element_ptr el,
							   int_fast64_t el_num,
							   const ElementParameters &eparams);

	struct TSun
	{
		TSun();
		void Reset();
		// void set_values(ray_source_ptr rsrc);

		// char ShapeIndex;
		SolTrace::Data::DistributionType ShapeIndex;
		double Sigma;
		bool PointSource;

		std::vector<double> SunShapeAngle;
		std::vector<double> SunShapeIntensity;
		double MaxAngle;
		double MaxIntensity;

		double Origin[3];

		// calculated
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
			// unsigned int raynum;
			SolTrace::Result::ray_id raynum;
			SolTrace::Result::RayEvent event;
		};

		ray_t *Append(double pos[3],
					  double cos[3],
					  int element,
					  int stage,
					  unsigned int raynum,
					  SolTrace::Result::RayEvent it);

		bool Overwrite(unsigned int idx,
					   double pos[3],
					   double cos[3],
					   int element,
					   int stage,
					   unsigned int raynum,
					   SolTrace::Result::RayEvent it);

		bool Query(unsigned int idx,
				   double pos[3],
				   double cos[3],
				   int *element,
				   int *stage,
				   unsigned int *raynum,
				   SolTrace::Result::RayEvent *it) const;

		void Merge(TRayData &dest);

		void Clear();

		void Print() const;

		uint_fast64_t Count() const;

		ray_t *Index(uint_fast64_t i, bool write_access) const;

	private:
		static const unsigned int block_size = 8192;

		struct block_t
		{
			ray_t data[block_size];
			uint_fast64_t count;
		};

		using block_t_ptr = std::shared_ptr<block_t>;

		std::vector<block_t_ptr> m_blockList;
		uint_fast64_t m_dataCount;
		uint_fast64_t m_dataCapacity;
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

		// std::vector<TElement*> ElementList;
		std::vector<telement_ptr> ElementList;
		// std::map<element_id, telement_ptr> ElementList;

		// calculated
		double Euler[3];
		double RRefToLoc[3][3];
		double RLocToRef[3][3];

		// TRayData RayData;

		int_fast64_t stage_id;
	};

	using tstage_ptr = typename std::shared_ptr<TStage>;
	tstage_ptr make_tstage(const ElementParameters &eparams);
	tstage_ptr make_tstage(SolTrace::Data::element_ptr el, const ElementParameters &eparams);

	struct TSystem
	{
		TSystem();
		~TSystem();

		void ClearAll();
		// void CollectResults();

		TSun Sun;
		std::vector<tstage_ptr> StageList;
		// std::map<element_id, std::pair<int_fast64_t, int_fast64_t> > id_map;

		// system simulation context data
		int sim_raycount;
		int sim_raymax;
		bool sim_dynamic_group; // point-focus heliostat dynamic grouping to reduce stage one computation
		bool sim_errors_sunshape;
		bool sim_errors_optical;

		uint_fast64_t seed;

		// simulation outputs
		// TRayData AllRayData;
		TRayData RayData;
		uint_fast64_t SunRayCount;

		std::vector<std::string> messages;

		void errlog(const char *fmt, ...);
	};

} // namespace SolTrace::NativeRunner

#endif
