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

#include <stdarg.h>

#include <cstring>
#include <memory>
#include <sstream>

#include "calculator_factory.hpp"
#include "matvec.hpp"
#include "native_runner_types.hpp"
#include "simulation_data_export.hpp"
#include "simulation_result_export.hpp"
#include "simulation_runner.hpp"

// TOpticalProperties::TOpticalProperties()
// {
//     for (int i = 0; i < 4; i++)
//         RefractiveIndex[i] = AB12[i] = 0;

//     // OpticSurfNumber = 1;
//     // ApertureStopOrGratingType = 0;
//     // DiffractionOrder = 0;
//     Reflectivity = 0;
//     Transmissivity = 0;
//     RMSSlopeError = 0;
//     RMSSpecError = 0;
//     DistributionType = 'g';
//     UseReflectivityTable = false;
//     UseTransmissivityTable = false;
// }

// TOpticalProperties &TOpticalProperties::operator=(const TOpticalProperties &rhs)
// {
//     DistributionType = rhs.DistributionType;
//     // OpticSurfNumber = rhs.OpticSurfNumber;
//     // ApertureStopOrGratingType = rhs.ApertureStopOrGratingType;
//     // DiffractionOrder = rhs.DiffractionOrder;
//     Reflectivity = rhs.Reflectivity;
//     Transmissivity = rhs.Transmissivity;
//     RMSSlopeError = rhs.RMSSlopeError;
//     RMSSpecError = rhs.RMSSpecError;
//     UseReflectivityTable = rhs.UseReflectivityTable;
//     UseTransmissivityTable = rhs.UseTransmissivityTable;

//     for (int i = 0; i < 4; i++)
//     {
//         RefractiveIndex[i] = rhs.RefractiveIndex[i];
//         AB12[i] = rhs.AB12[i];
//     }

//     ReflectivityTable.resize(rhs.ReflectivityTable.size());
//     for (size_t i = 0; i < rhs.ReflectivityTable.size(); i++)
//     {
//         ReflectivityTable[i].angle = rhs.ReflectivityTable[i].angle;
//         ReflectivityTable[i].refl = rhs.ReflectivityTable[i].refl;
//     }
//     TransmissivityTable.resize(rhs.TransmissivityTable.size());
//     for (size_t i = 0; i < rhs.TransmissivityTable.size(); i++)
//     {
//         TransmissivityTable[i].angle = rhs.TransmissivityTable[i].angle;
//         TransmissivityTable[i].trans = rhs.TransmissivityTable[i].trans;
//     }

//     return *this;
// }

namespace SolTrace::NativeRunner
{

    ElementParameters::ElementParameters()
        : newton_tolerance(1e-6),
          newton_max_iters(20)
    {
    }

    ElementParameters::~ElementParameters()
    {
    }

    TElement::TElement() : aperture(nullptr),
                           icalc(nullptr),
                           Optics()
    {
        int i, j;
        for (i = 0; i < 3; i++)
        {
            // Origin[i] = AimPoint[i] = Euler[i] = PosSunCoords[i] = 0;
            Origin[i] = AimPoint[i] = PosSunCoords[i] = 0;
            // PosSunCoords[i] = 0;
        }
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                RRefToLoc[i][j] = RLocToRef[i][j] = 0;

        // ZRot = 0;
        ZAperture = 0;

        // Optics = nullptr;
        element_number = -1; // mjw nonsense
        sim_data_id = ELEMENT_ID_UNASSIGNED;
        // stage_index = -1;
    }

    TElement::~TElement()
    {
        aperture = nullptr;
        // Optics = nullptr;
    }

    TSun::TSun()
    {
        Reset();
    }

    void TSun::Reset()
    {
        int i, j;
        for (i = 0; i < 3; i++)
            Origin[i] = Euler[i] = 0;
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                RRefToLoc[i][j] = RLocToRef[i][j] = 0;

        PointSource = false;
        // ShapeIndex = ' ';
        ShapeIndex = SunShape::GAUSSIAN;
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

    // void TSun::set_values(ray_source_ptr rsrc)
    // {
    //     CopyVec3(this->Origin, rsrc->get_position().data);

    //     // TODO: Need to get sun parameters here too.
    //     this->PointSource = false;

    //     return;
    // }

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
        : nthreads(1),
          nray_per_thread(0),
          nray_remainder(0)
    {
        this->Clear();
    }

    TRayData::~TRayData()
    {
        this->Clear();
    }

    void TRayData::SetUp(unsigned nthreads, uint_fast64_t nrays)
    {
        this->Clear();
        this->nthreads = nthreads;
        this->nray_per_thread = nrays / nthreads;
        this->nray_remainder = nrays % nthreads;
        for (unsigned k = 0; k < nthreads; ++k)
        {
            this->records[k].clear();
        }
        // std::cout << "nthreads: " << nthreads
        //           << "  nray_per_thread: " << nray_per_thread
        //           << std::endl;
        return;
    }

    TRayData::ray_t_ptr TRayData::Append(unsigned thread_id,
                                         double pos[3],
                                         double cos[3],
                                         int element,
                                         int stage,
                                         uint_fast64_t raynum,
                                         SolTrace::Result::RayEvent rev)
    {
        // std::stringstream ss;
        // ss << "Thread: " << thread_id << " Raynum: " << raynum
        //    << " RayID: " << this->GetRayId(thread_id, raynum)
        //    << std::endl;
        // std::cout << ss.str();

        ray_t_ptr r = this->GetNext(thread_id);

        if (r != nullptr)
        {
            std::memcpy(&r->pos, pos, sizeof(double) * 3);
            std::memcpy(&r->cos, cos, sizeof(double) * 3);
            r->element = element;
            r->stage = stage;
            r->raynum = this->GetRayId(thread_id, raynum);
            r->event = rev;
        }

        return r;
    }

    TRayData::ray_t_ptr TRayData::Index(uint_fast64_t idx) const
    {
        uint_fast64_t n = idx;
        uint_fast64_t nevents = 0;
        unsigned k = 0;
        ray_t_ptr r = nullptr;
        while (k < this->nthreads)
        {
            auto iter = this->records.find(k);
            assert(iter != this->records.end());
            nevents = iter->second.size();
            // std::cout << "k: " << k
            //           << "  n: " << n
            //           << "  nevents: " << nevents
            //           << std::endl;
            if (n >= nevents)
            {
                n -= nevents;
            }
            else
            {
                r = iter->second[n];
                break;
            }
            ++k;
        }
        return r;
    }

    bool TRayData::Query(uint_fast64_t idx,
                         double pos[3],
                         double cos[3],
                         int *element,
                         int *stage,
                         uint_fast64_t *raynum,
                         SolTrace::Result::RayEvent *rev) const
    {

        ray_t_ptr r = this->Index(idx);

        if (r != nullptr)
        {
            if (pos != nullptr)
                std::memcpy(pos, r->pos, sizeof(double) * 3);
            if (cos != nullptr)
                std::memcpy(cos, r->cos, sizeof(double) * 3);
            if (element != nullptr)
                *element = r->element;
            if (stage != nullptr)
                *stage = r->stage;
            if (raynum != nullptr)
                *raynum = r->raynum;
            if (rev != nullptr)
                *rev = r->event;
        }

        return r != nullptr;
    }

    void TRayData::Clear()
    {
        this->records.clear();
        return;
    }

    uint_fast64_t TRayData::Count() const
    {
        uint_fast64_t count = 0;
        for (auto iter = this->records.cbegin();
             iter != this->records.cend();
             ++iter)
        {
            count += iter->second.size();
        }
        return count;
    }

    TRayData::ray_t_ptr TRayData::GetNext(unsigned thread_id)
    {
        ray_t_ptr r = nullptr;
        // auto n = this->records[thread_id].size();
        // this->records[thread_id].push_back(r);
        auto it = this->records.find(thread_id);
        if (it != this->records.end())
        {
            r = std::make_shared<ray_t>();
            it->second.push_back(r);
        }
        return r;
    }

    uint_fast64_t TRayData::GetRayId(unsigned thread_id,
                                     uint_fast64_t raynum)
    {
        uint_fast64_t rayid = thread_id * this->nray_per_thread + raynum;
        rayid += std::min(static_cast<uint_fast64_t>(thread_id),
                          this->nray_remainder);
        return rayid;
    }

    void TRayData::Print() const
    {
        uint_fast64_t n = this->Count();

        for (uint_fast64_t i = 0; i < n; ++i)
        {
            double pos[3], cos[3];
            int elm, stage;
            uint_fast64_t ray;
            SolTrace::Result::RayEvent rev;
            if (Query(i, pos, cos, &elm, &stage, &ray, &rev))
            {
                printf("   [%llu] = { [%lg,%lg,%lg][%lg,%lg,%lg] %d %d %llu %s\(%d) }\n",
                       i,
                       pos[0], pos[1], pos[2],
                       cos[0], cos[1], cos[2],
                       elm,
                       stage,
                       static_cast<long long unsigned>(ray),
                       ray_event_string(rev).c_str(),
                       static_cast<int>(rev));
            }
        }

        return;
    }

    TStage::TStage()
    {
        uint_fast64_t i, j;
        for (i = 0; i < 3; i++)
            Origin[i] = AimPoint[i] = Euler[i] = 0;
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                RRefToLoc[i][j] = RLocToRef[i][j] = 0;

        ZRot = 0;
        MultiHitsPerRay = true;
        Virtual = false;
        TraceThrough = false;

        ElementList.clear();
    }

    TStage::~TStage()
    {
        ElementList.clear();
    }

    TSystem::TSystem()
    {
        SunRayCount = 0;

        sim_raycount = 1000;
        sim_raymax = 100000;
        sim_dynamic_group = true;
        sim_errors_sunshape = true;
        sim_errors_optical = true;

        ClearAll();
    }

    TSystem::~TSystem()
    {
        ClearAll();
    }

    void TSystem::ClearAll()
    {
        StageList.clear();
        Sun.Reset();
        this->RayData.Clear();
    }

    telement_ptr make_telement(element_ptr el,
                               int_fast64_t el_num,
                               const ElementParameters &eparams)
    {
        // std::cout << "Name: " << el->get_name()
        //           << "\nSDID: " << el->get_id()
        //           << "\nNum: " << el_num
        //           << std::endl;

        telement_ptr telem = std::make_shared<TElement>();
        vector_copy(telem->Origin, el->get_origin_stage());
        vector_copy(telem->AimPoint, el->get_aim_vector_stage());
        telem->ZRot = el->get_zrot();

        // vector_copy(telem->Euler, el->get_euler_angles());
        matrix_copy(telem->RRefToLoc, el->get_stage_to_local());
        matrix_copy(telem->RLocToRef, el->get_local_to_stage());

        telem->aperture = el->get_aperture()->make_copy();
        telem->icalc =
            CalculatorFactory::get()->make_calculator(telem->aperture,
                                                      el->get_surface(),
                                                      eparams);

        // How to handle optical properties?
        telem->Optics.Front = *el->get_front_optical_properties();
        telem->Optics.Back = *el->get_back_optical_properties();

        // telem->element_number = el->get_id();
        telem->sim_data_id = el->get_id();
        telem->element_number = el_num;

        return telem;
    }

    tstage_ptr make_tstage(const ElementParameters &eparams)
    {
        tstage_ptr my_stage = std::make_shared<TStage>();

        // Use global coordinates as stage coordinates
        ZeroVec3(my_stage->Origin);
        ZeroVec3(my_stage->AimPoint);
        my_stage->AimPoint[2] = 1.0;
        my_stage->ZRot = 0.0;
        IdentityMat3(my_stage->RRefToLoc);
        IdentityMat3(my_stage->RLocToRef);

        return my_stage;
    }

    tstage_ptr make_tstage(element_ptr el,
                           const ElementParameters &eparams)
    {
        tstage_ptr my_stage = std::make_shared<TStage>();
        auto stage_el = std::dynamic_pointer_cast<StageElement>(el);
        // TODO: Throw an error...
        assert(stage_el != nullptr);

        // TODO: What to do with these fields?
        my_stage->MultiHitsPerRay = true;
        my_stage->Virtual = stage_el->is_virtual();
        my_stage->TraceThrough = true;
        my_stage->stage_id = stage_el->get_stage();

        // Add coordinate stuff
        vector_copy(my_stage->Origin, stage_el->get_origin_global());
        vector_copy(my_stage->AimPoint, stage_el->get_aim_vector_global());
        my_stage->ZRot = stage_el->get_zrot();
        matrix_copy(my_stage->RRefToLoc, stage_el->get_global_to_local());
        matrix_copy(my_stage->RLocToRef, stage_el->get_local_to_global());

        return my_stage;
    }

} // namespace SolTrace::NativeRunner
