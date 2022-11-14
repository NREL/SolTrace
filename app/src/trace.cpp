
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


#include <thread>
#include <mutex>
#include <algorithm>

#include <wx/wx.h>
#include <wx/dirdlg.h>
#include <wx/filename.h>
#include <wx/time.h>
#include <wx/datetime.h>
#include <wx/thread.h>

#include <wex/exttext.h>
#include <wex/numeric.h>
#include <wex/metro.h>
#include <wex/tpdlg.h>
#include <wex/utils.h>

#include "trace.h"
#include "soltrace.h"

static wxThreadProgressDialog *g_currentThreadProgress = NULL;

enum{ ID_NUM_RAYS = wxID_HIGHEST+923 };

BEGIN_EVENT_TABLE( TraceForm, wxPanel )
	EVT_BUTTON( wxID_SETUP, TraceForm::OnCommand )
	EVT_BUTTON( wxID_EXECUTE, TraceForm::OnCommand )
	EVT_NUMERIC( ID_NUM_RAYS, TraceForm::OnCommand )
END_EVENT_TABLE()

TraceForm::TraceForm( wxWindow *parent, Project &prj )
	: wxPanel( parent ), m_prj( prj )
{

	wxStaticBoxSizer *sizer1 = new wxStaticBoxSizer( wxVERTICAL, this, "Parameters" );
	wxFlexGridSizer *flxsizer = new wxFlexGridSizer( 2, wxSize(2,2) );
	flxsizer->Add( new wxStaticText( sizer1->GetStaticBox(), wxID_ANY, "Desired number of ray intersections"), 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 3 );
	flxsizer->Add( m_numRays = new wxNumericCtrl( sizer1->GetStaticBox(), ID_NUM_RAYS, 10000, wxNUMERIC_UNSIGNED ), 0, wxALL, 0 );
	m_numRays->SetFormat( wxNUMERIC_GENERIC, true );
	
	flxsizer->Add( new wxStaticText( sizer1->GetStaticBox(), wxID_ANY, "Maximum number of generated sun rays"), 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 3 );
	flxsizer->Add( m_numMaxSunRays = new wxNumericCtrl( sizer1->GetStaticBox(), wxID_ANY, 1000000, wxNUMERIC_UNSIGNED ), 0, wxALL, 0 );
	m_numMaxSunRays->SetFormat( wxNUMERIC_GENERIC, true );
	
	flxsizer->Add( new wxStaticText( sizer1->GetStaticBox(), wxID_ANY, "Maximum number of CPUs to utilize"), 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 3 );
	flxsizer->Add( m_numCpus = new wxNumericCtrl( sizer1->GetStaticBox(), wxID_ANY, 16, wxNUMERIC_INTEGER ), 0, wxALL, 0 );
	
	flxsizer->Add( new wxStaticText( sizer1->GetStaticBox(), wxID_ANY, "Seed value (-1 for random)"), 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 3 );
	flxsizer->Add( m_seed = new wxNumericCtrl( sizer1->GetStaticBox(), wxID_ANY, 123, wxNUMERIC_INTEGER ), 0, wxALL, 0 );

	flxsizer->Add( m_inclSunShape = new wxCheckBox( sizer1->GetStaticBox(), wxID_ANY, "Include sun shape" ), 0, wxALL|wxALIGN_CENTER_VERTICAL, 3 );
	flxsizer->AddStretchSpacer();
	flxsizer->Add( m_inclOpticalErrors = new wxCheckBox( sizer1->GetStaticBox(), wxID_ANY, "Include optical errors" ), 0, wxALL|wxALIGN_CENTER_VERTICAL, 3 );
	flxsizer->AddStretchSpacer();
	flxsizer->Add( m_asPowerTower      = new wxCheckBox( sizer1->GetStaticBox(), wxID_ANY, "Point-focus system" ), 0, wxALL|wxALIGN_CENTER_VERTICAL, 3 );
	flxsizer->AddStretchSpacer();

	sizer1->Add( flxsizer, 0, wxALL, 5 );

	wxStaticBoxSizer *sizer2 = new wxStaticBoxSizer( wxHORIZONTAL, this, "Working folder" );
	sizer2->Add( m_workDir = new wxExtTextCtrl( sizer2->GetStaticBox(), wxID_ANY ), 1, wxALL|wxALIGN_CENTER_VERTICAL, 4 );
	sizer2->Add( new wxButton( sizer2->GetStaticBox(), wxID_SETUP, "...", wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT ), 0, wxALL, 2 ); 
	
	wxBoxSizer *vsizer1 = new wxBoxSizer( wxVERTICAL );
	vsizer1->Add( sizer1, 0, wxALL, 10 );	
	vsizer1->Add( sizer2, 0, wxEXPAND|wxALL, 10 );


	wxStaticBoxSizer *vsizer2 = new wxStaticBoxSizer( wxVERTICAL, this, "Trace" );
	
	wxMetroButton *tracebtn = new wxMetroButton( vsizer2->GetStaticBox(), wxID_EXECUTE, "Start ray trace", wxNullBitmap, wxDefaultPosition, wxDefaultSize, wxMB_RIGHTARROW );
	wxFont font( tracebtn->GetFont() );
	font.SetPointSize( font.GetPointSize() + 2 );
	tracebtn->SetFont( font );
	vsizer2->Add( tracebtn, 0, wxALL|wxALIGN_CENTER, 10 );

	wxFlexGridSizer *lstsizer = new wxFlexGridSizer( 2, wxSize(2,2) );
	lstsizer->Add( new wxStaticText( vsizer2->GetStaticBox(), wxID_ANY, "Elapsed time for last trace"), 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 3 );
	lstsizer->Add( m_elapsedTime = new wxNumericCtrl( vsizer2->GetStaticBox(), wxID_ANY ), 0, wxALL, 0 );
	m_elapsedTime->SetFormat( 3, false, wxEmptyString, " sec" );
	m_elapsedTime->SetEditable( false );

	wxColour cback( 245, 245, 245 );
	m_elapsedTime->SetBackgroundColour( cback );

	lstsizer->Add( new wxStaticText( vsizer2->GetStaticBox(), wxID_ANY, "Seed value used for last trace"), 0, wxALL|wxALIGN_CENTER_VERTICAL|wxALIGN_RIGHT, 3 );
	lstsizer->Add( m_lastSeed = new wxNumericCtrl( vsizer2->GetStaticBox(), wxID_ANY ), 0, wxALL, 0 );
	m_lastSeed->SetEditable( false );
	m_lastSeed->SetBackgroundColour( cback );

	vsizer2->Add( lstsizer, 0, wxALIGN_RIGHT|wxALL, 10 );


	wxBoxSizer *hsizer = new wxBoxSizer( wxHORIZONTAL );
	hsizer->Add( vsizer1 );
	hsizer->Add( vsizer2, 0, wxALL, 10 );
	SetSizer( hsizer );

}


void TraceForm::SetOptions( size_t nrays, size_t nmaxsunrays, int ncpu, int seed,
	bool sunshape, bool opterr, bool aspowertower )
{
	m_numRays->SetValue( nrays );
	m_numMaxSunRays->SetValue( nmaxsunrays );
	m_numCpus->SetValue( ncpu );
	m_seed->SetValue( seed );
	m_inclSunShape->SetValue( sunshape );
	m_inclOpticalErrors->SetValue( opterr );
    m_asPowerTower->SetValue( aspowertower );
}

void TraceForm::GetOptions( size_t *nrays, size_t *nmaxsunrays, int *ncpu, int *seed,
	bool *sunshape, bool *opterr, bool *aspowertower )
{
	if ( nrays ) *nrays = m_numRays->AsUnsigned();
	if ( nmaxsunrays ) *nmaxsunrays = m_numMaxSunRays->AsUnsigned();
	if ( ncpu ) *ncpu = m_numCpus->AsInteger();
	if ( seed ) *seed = m_seed->AsInteger();
	if ( sunshape ) *sunshape = m_inclSunShape->GetValue();
	if ( opterr ) *opterr = m_inclOpticalErrors->GetValue();
    if ( aspowertower ) *aspowertower = m_asPowerTower->GetValue();
}


void TraceForm::SetWorkDir( const wxString &path )
{
	m_workDir->ChangeValue( path );
}

wxString TraceForm::GetWorkDir()
{
	return m_workDir->GetValue();
}

void TraceForm::OnCommand( wxCommandEvent &evt )
{
	switch( evt.GetId() )
	{
	case ID_NUM_RAYS:
		if ( m_numMaxSunRays->Value() < 100*m_numRays->Value() )
			m_numMaxSunRays->SetValue( 100*m_numRays->Value() );
		break;

	case wxID_EXECUTE:
		StartTrace();
		break;

	case wxID_SETUP:
		{
			wxString dir( wxDirSelector( "Select working folder", m_workDir->GetValue() ) );
			if ( !dir.IsEmpty() )
				m_workDir->ChangeValue( dir );		
		}
		break;
	}
}

int TraceForm::StartTrace( bool , bool quiet, wxArrayString *err )
{
	if ( g_currentThreadProgress != 0 )
	{
		if ( !quiet )
			wxMessageBox("Ray trace in progress.");

		return -999;
	}


	m_lastSeedVal = m_seed->AsInteger();
	
static wxArrayString s_err;
	wxArrayString &ref_errors = err ? *err : s_err;

	ref_errors.clear();

	int msec = RunTraceMultiThreaded( &m_prj, m_numRays->AsInteger(),
			m_numMaxSunRays->AsInteger(),
			m_numCpus->AsInteger(),
			&m_lastSeedVal,
			m_inclSunShape->GetValue(),
			m_inclOpticalErrors->GetValue(),
            m_asPowerTower->GetValue(),
			ref_errors );

	if ( msec < 0 )
		wxShowTextMessageDialog( wxJoin( ref_errors, '\n' ) );

	m_elapsedTime->SetValue( msec*0.001 );
	m_lastSeed->SetValue( m_lastSeedVal );

	MainWindow::Instance().UpdateResults();
	return msec;
}

bool TraceForm::IsRunning()
{
	return g_currentThreadProgress != 0;
}

void TraceForm::CancelTrace()
{
	if ( g_currentThreadProgress )
		g_currentThreadProgress->Cancel();
}

static int LoadSystemIntoContext( Project *System, st_context_t spcxt, wxArrayString &errs )
{
	int errflag = 0;

	/* configure sun */
	char shape = 'i'; // invalid
	double sigma = 0.0;
	if (System->Sun.Shape == SunShape::GAUSSIAN)
	{
		shape = 'g';
		sigma = System->Sun.Sigma;
	}
	else if (System->Sun.Shape == SunShape::PILLBOX)
	{
		shape = 'p';
		sigma = System->Sun.HalfWidth;
	}
	else if (System->Sun.Shape == SunShape::USER_DEFINED)
	{
		shape = 'd';
	}
	else
	{
		errs.Add( wxString("Invalid sun shape: ") + shape);
		return -55;
	}

	//qDebug("Loaded SunShape = %c, shw=%lg", shape, sigma);
	st_sun(spcxt, System->Sun.PointSource?1:0, shape, sigma );

	double x,y,z;
	if (System->Sun.UseLDHSpec)
	{
		double lat, day, hour;
		double Declination, HourAngle, Elevation, Azimuth;

		lat = System->Sun.Latitude;
		day = System->Sun.Day;
		hour = System->Sun.Hour;

		Declination = 180/M_PI*asin(0.39795*cos(0.98563*M_PI/180*(day-173)));
		HourAngle = 15*(hour-12);
		Elevation = 180/M_PI*asin(sin(Declination*M_PI/180)*sin(lat*M_PI/180)+cos(Declination*M_PI/180)*cos(HourAngle*M_PI/180)*cos(lat*M_PI/180));
		Azimuth = 180/M_PI*acos((sin(M_PI/180*Declination)*cos(M_PI/180*lat)-cos(M_PI/180*Declination)*sin(M_PI/180*lat)*cos(M_PI/180*HourAngle))/cos(M_PI/180*Elevation)+0.0000000001);
		if ( sin(HourAngle*M_PI/180) > 0.0 )
			Azimuth = 360 - Azimuth;
		x = -sin(Azimuth*M_PI/180)*cos(Elevation*M_PI/180);
		y = sin(Elevation*M_PI/180);
		z = cos(Azimuth*M_PI/180)*cos(Elevation*M_PI/180);
	}
	else
	{
		x = System->Sun.X;
		y = System->Sun.Y;
		z = System->Sun.Z;
	}

	st_sun_xyz(spcxt, x, y, z );

	int sun_npoints = System->Sun.UserShapeData.size();
	if (sun_npoints > 0)
	{
		double *angle = new double[sun_npoints];
		double *intensity = new double[sun_npoints];

		for (int i=0;i< sun_npoints;i++)
		{
			angle[i] = System->Sun.UserShapeData[i].x;
			intensity[i] = System->Sun.UserShapeData[i].y;
		}

		st_sun_userdata(spcxt, sun_npoints, angle, intensity );

		delete [] angle;
		delete [] intensity;
	}

	st_clear_optics(spcxt);
	for (size_t nopt=0;nopt<System->OpticsList.size();nopt++)
	{
		int idx = st_add_optic(spcxt, (const char*)System->OpticsList[nopt]->Name.c_str() );

		// iterate over the front and back surfaces using a pointer to each
		SurfaceOptic* surfs[2] = { &System->OpticsList[nopt]->Front, &System->OpticsList[nopt]->Back };

		for (size_t si = 0; si < 2; si++)
		{
			//si==0 -> Front, si==1 -> Back
			SurfaceOptic* f = surfs[si];

			double *refl_angles = 0;
			double *refls = 0;
			int refl_npoints = f->ReflectivityTable.size();

			if (f->UseReflectivityTable && refl_npoints > 0)
			{
				refl_angles = new double[refl_npoints];
				refls = new double[refl_npoints];
				for (int k = 0; k < refl_npoints; k++)
				{
					refl_angles[k] = f->ReflectivityTable[k].x;
					refls[k] = f->ReflectivityTable[k].y;
				}
			}

			double* trans_angles = 0;
			double* transs = 0;
			int trans_npoints = f->TransmissivityTable.size();

			if (f->UseTransmissivityTable && trans_npoints > 0)
			{
				trans_angles = new double[trans_npoints];
				transs = new double[trans_npoints];
				for (int k = 0; k < trans_npoints; k++)
				{
					trans_angles[k] = f->TransmissivityTable[k].x;
					transs[k] = f->TransmissivityTable[k].y;
				}

			}

			st_optic(spcxt, idx, si+1, f->ErrorDistribution,
				f->OpticalSurfaceNumber, f->ApertureStopOrGratingType, f->DiffractionOrder,
				f->RefractionIndexReal, f->RefractionIndexImag,
				f->Reflectivity, f->Transmissivity,
				f->GratingCoeffs, f->RMSSlope, f->RMSSpecularity,
				f->UseReflectivityTable ? 1 : 0, refl_npoints,
				refl_angles, refls,
				f->UseTransmissivityTable ? 1 : 0, trans_npoints,
				trans_angles, transs
			);

			if (refl_angles) delete[] refl_angles;
			if (refls) delete[] refls;
			if (trans_angles) delete[] trans_angles;
			if (transs) delete[] transs;

			refl_angles = 0;
			refls = 0;
			trans_angles = 0;
			transs = 0;
		}
	}

	st_clear_stages(spcxt);
	st_add_stages(spcxt, System->StageList.size());

	for (size_t ns=0;ns<System->StageList.size();ns++)
	{
		Stage *stage = System->StageList[ns];

		::st_stage_flags( spcxt, ns, stage->Virtual?1:0, stage->MultiHit?1:0, stage->TraceThrough?1:0);

		st_stage_xyz(spcxt, ns, stage->X, stage->Y, stage->Z );
		st_stage_aim(spcxt, ns, stage->AX, stage->AY, stage->AZ );
		st_stage_zrot(spcxt, ns, stage->ZRot );

		st_clear_elements(spcxt, ns);
		st_add_elements( spcxt, ns, stage->ElementList.size() );

		for (size_t idx=0;idx<stage->ElementList.size();idx++)
		{
			Element *e = stage->ElementList[idx];

			st_element_enabled( spcxt, ns, idx, e->Enabled?1:0 );
			st_element_xyz( spcxt, ns, idx, e->X, e->Y, e->Z );
			st_element_aim( spcxt, ns, idx, e->AX, e->AY, e->AZ );
			st_element_zrot( spcxt, ns, idx, e->ZRot );
			st_element_aperture( spcxt, ns, idx, e->ApertureIndex );
			st_element_aperture_params( spcxt, ns, idx, e->ApertureParams );
			st_element_surface( spcxt, ns, idx, e->SurfaceIndex );
			if (!e->SurfaceFile.IsEmpty())
			{
				wxString sf( e->SurfaceFile );				
				wxFileName fn( sf );
				if (!fn.IsAbsolute())
				{
					sf = wxFileName( MainWindow::Instance().GetFileName() ).GetPath() + "/" + e->SurfaceFile;
					if (!wxFileExists(sf))
					{
						sf = MainWindow::Instance().GetWorkDir() + "/" + e->SurfaceFile;
						if (!wxFileExists(sf))
						{
							errs.Add( "Could not locate surface file: " + e->SurfaceFile );
							errflag = -4;
							break;
						}
					}
				}

				if ( st_element_surface_file( spcxt, ns, idx, (const char*)sf.c_str() ) < 0)
				{
					errflag = -3;
					for (int j=0;j<st_num_messages(spcxt);j++)
							errs.Add( st_message(spcxt, j) );
					break;
				}
			}
			else
				st_element_surface_params( spcxt, ns, idx, e->SurfaceParams );

			int interact = 1; // 1= refract, 2=reflect
			if (e->InteractionType == Element::REFLECTION ) interact = 2;
			st_element_interaction( spcxt, ns, idx, interact );

			st_element_optic( spcxt, ns, idx, (const char*)e->OpticName.c_str() );
		}

	}

	return errflag;
}



static void CheckSeed( int *seed )
{
static bool first_validate = true;
	if( first_validate )
	{
		::srand( ::time( NULL ) );
		first_validate = false;
	}

	if ( *seed < 1 )
		*seed = ( ::rand() % (RAND_MAX/3) + 813 );
}



void CountRayHitsPerElement( Project *System)
{
	for (size_t j=0;j<System->StageList.size();j++)
	{
		System->StageList[j]->RayHits = 0.0;
		for (size_t i=0;i<System->StageList[j]->ElementList.size();i++)
		{
			System->StageList[j]->ElementList[i]->RayHits = 0;
			System->StageList[j]->ElementList[i]->FinalRayHits = 0;
		}
	}

	for (size_t i=0;i<System->Results.Length;i++)
	{
		int stage_num = System->Results.StageMap[i] - 1;
		int element_num = abs(System->Results.ElementMap[i]) - 1;


		if (stage_num >= 0 && stage_num < (int)System->StageList.size()
			&& element_num >= 0 && element_num < (int)System->StageList[stage_num]->ElementList.size() )
		{
			System->StageList[stage_num]->ElementList[element_num]->RayHits++;
			System->StageList[stage_num]->RayHits++;

			// count up final ray hits
			if (System->Results.ElementMap[i] < 0)
				System->StageList[stage_num]->ElementList[element_num]->FinalRayHits++;
		}
	}
}

// forward decl
static int trace_callback_multi_thread( 
	st_uint_t ntracedtotal, st_uint_t ntraced, 
	st_uint_t ntotrace, st_uint_t curstage, 
	st_uint_t nstages, void *data );

class TraceThread : public wxThread
{
private:
	st_context_t m_contextId;
	bool m_cancelFlag;
    bool m_asPowerTower;
	bool m_sunshape;
	bool m_opterrs;
	int m_nrays;
	int m_nmaxrays;

	size_t m_nTraceTotal;
	size_t m_nTraced;
	size_t m_nToTrace;
	size_t m_curStage;
	size_t m_nStages;
	int m_iThread;
	int m_seedVal;
	int m_resultCode;

	Project* m_system;
	wxArrayString* m_errmsg;

	wxMutex m_statusLock;
public:
	TraceThread( Project* system, st_context_t spcxt, wxArrayString* errmsg, int ithread, int seed, bool aspowertower, int nrays, int nmaxrays, bool sunshape, bool opterrs )
		: wxThread( wxTHREAD_JOINABLE ), m_cancelFlag( false )
	{
		m_iThread = ithread;
		m_contextId = spcxt;
		m_seedVal = seed;
		m_resultCode = -1;
        m_asPowerTower = aspowertower;
		m_nrays = nrays;
		m_nmaxrays = nmaxrays;
		m_sunshape = sunshape;
		m_opterrs = opterrs;
		m_errmsg = errmsg;

		m_system = system;

		m_nTraceTotal = m_nTraced = m_nToTrace = m_curStage = m_nStages = 0;
	}

	virtual ~TraceThread()
	{
		::st_free_context( m_contextId );
	}

	void cancelTrace()
	{
		m_cancelFlag = 1;
	}

	bool isTraceCanceled()
	{
		return (m_cancelFlag != 0);
	}

	st_context_t contextId()
	{
		return m_contextId;
	}

	int resultCode()
	{
		return m_resultCode;
	}

	void updateStatus( size_t ntracedtotal, size_t ntraced, size_t ntotrace, size_t curstage, size_t nstages )
	{
		wxMutexLocker lock( m_statusLock );

		m_nTraceTotal = ntracedtotal;
		m_nTraced = ntraced;
		m_nToTrace = ntotrace;
		m_curStage = curstage;
		m_nStages = nstages;
	}

	void status(size_t *total, size_t *traced, size_t *totrace, size_t *stage, size_t *nstages )
	{
		wxMutexLocker lock( m_statusLock );

		*total = m_nTraceTotal;
		*traced = m_nTraced;
		*totrace = m_nToTrace;
		*stage = m_curStage;
		*nstages = m_nStages;
	}
	
	virtual ExitCode Entry()
	{
		::st_sim_errors(m_contextId, m_sunshape ? 1 : 0, m_opterrs ? 1 : 0);
		::st_sim_params(m_contextId, m_nrays, m_nmaxrays);

		m_resultCode = ::st_sim_run( m_contextId, 
			(unsigned int) m_seedVal,
            m_asPowerTower,
			trace_callback_multi_thread, this );

		return 0;
	}
};


int trace_callback_multi_thread(st_uint_t ntracedtotal, st_uint_t ntraced, st_uint_t ntotrace, st_uint_t curstage, st_uint_t nstages, void *data)
{
	TraceThread *t = (TraceThread*)data;
	t->updateStatus( ntracedtotal, ntraced, ntotrace, curstage, nstages );
	return t->isTraceCanceled() ? 0 : 1; // if it was canceled, stop processing by return 0
}

int RunTraceMultiThreaded( Project *System, int nrays, int nmaxrays,
						  int nmaxthreads, int *seed, bool sunshape, bool opterrs, bool aspowertower,
						  wxArrayString &errors, bool is_cmd )
{
	if (nmaxthreads < 1)
	{
		errors.Add( "invalid number of cpu threads" );
		return -888;
	}

	CheckSeed( seed );

	size_t ncpus = wxThread::GetCPUCount();
	if (nmaxthreads >= 1 && ncpus > (size_t)nmaxthreads) ncpus = (size_t)nmaxthreads;

    wxThreadProgressDialog *tpd = 0;

    if( is_cmd )
    {
        wxPrintf("\nRunning simulation with %d threads...", (int)ncpus);
    }
    else
    {
        tpd = new wxThreadProgressDialog( &MainWindow::Instance(), ncpus, true );
	    tpd->CenterOnParent();
	    tpd->Show();
    }
		
	std::vector<TraceThread*> ThreadList;
	int SeedVal = *seed;
	wxStopWatch sw;

	bool ok = true;

	for (size_t i = 0; i < ncpus && ok == true; i++)
	{
		st_context_t spcxt = ::st_create_context();

		int result = LoadSystemIntoContext(System, spcxt, errors);
		if (result < 0)
		{
			errors.Add("error loading system into simulation context");
			ok = false;
			continue;
		}

		int rays_this_thread = nrays / ncpus;
		if (i == 0) rays_this_thread += (nrays % ncpus);

		::st_sim_errors(spcxt, sunshape ? 1 : 0, opterrs ? 1 : 0);
		::st_sim_params(spcxt, rays_this_thread, nmaxrays);
		SeedVal += i * 123;

		ThreadList.push_back(new TraceThread(System, spcxt, &errors, i, SeedVal, aspowertower, rays_this_thread, nmaxrays, sunshape, opterrs));
	}

	g_currentThreadProgress = tpd;
	for (size_t i=0;i<ThreadList.size();i++)
	{
		ThreadList[i]->Create();
		ThreadList[i]->Run();
	}

	size_t ntotal, ntraced, ntotrace, stagenum, nstages;

	while (1)
	{
		size_t num_finished = 0;
		for (size_t i=0;i<ThreadList.size();i++)
			if ( !ThreadList[i]->IsRunning() )
				num_finished++;

		if (num_finished == ThreadList.size())
			break;

		int ntotaltraces = 0;
        wxString cmd_threadstate; cmd_threadstate.clear();
		
        for (size_t i=0;i<ThreadList.size();i++)
		{
			ThreadList[i]->status(&ntotal, &ntraced, &ntotrace, &stagenum, &nstages);
			if( is_cmd )
            {
                cmd_threadstate += wxString::Format("%5.1f", 100.0f*((float)ntraced)/((float)std::max(1,(int)ntotrace)) );

                if( i==ThreadList.size()-1 )
                    wxPrintf( "\nStage %d of %d: %s", (int)stagenum, (int)nstages, cmd_threadstate.c_str() );
            }
            else
            {
                tpd->Update( i, 100.0f*((float)ntraced)/((float)std::max(1,(int)ntotrace)), wxString::Format("Stage %d of %d", (int)stagenum, (int)nstages) );
            }
			ntotaltraces += ntotal;
		}
		
		// need to process events
		wxYield();

        if( !is_cmd )
        {
		    if (tpd->IsCanceled())
		    {
			    for (size_t i=0;i<ThreadList.size();i++)
				    ThreadList[i]->cancelTrace();
		    }
        }

		// sleep a little
		wxMilliSleep( is_cmd ? 200 : 50 );
	}
	
	// wait on the joinable threads
	// make sure all have finished before continuing
	for (size_t i=0;i<ThreadList.size();i++)
		ThreadList[i]->Wait();

	bool errors_found = false;
	for (size_t i=0;i<ThreadList.size();i++)
	{
		st_context_t cxt = ThreadList[i]->contextId();
		int code = ThreadList[i]->resultCode();

		if (code < 0)
		{
			errors_found = true;
			for ( int j=0;j<st_num_messages(cxt);j++ )
				errors.Add( st_message(cxt, j ) );

			errors.Add( wxString::Format("error in trace thread %d, code %d", (int)(i+1), code ));
		}
	}

	int millisec = (int)sw.Time();

	// recompute transforms for results viewing
	System->RecomputeTransforms();

	// aggregate all the output data from the threads
	if (!errors_found)
	{
		std::vector<st_context_t> ContextList;
		for (size_t i=0;i<ThreadList.size();i++)
			ContextList.push_back( ThreadList[i]->contextId() );

		if (!System->Results.ReadResultsFromContextList( ContextList ))
		{
			errors.Add( "Allocation error reading results from trace context - (mt)" );
			errors_found = true;
		}

		CountRayHitsPerElement( System );
	}
	else
		System->Results.FreeMemory();

	for (size_t i=0;i<ThreadList.size();i++)
		delete ThreadList[i];

	ThreadList.clear();
	
	g_currentThreadProgress = NULL;
	
	if(! is_cmd )
    {
	    if ( tpd->IsCanceled() )
	    {
            delete tpd;
		    errors.Add("ray trace canceled by user");
		    return -5;
	    }

        delete tpd;
    }
		
	return errors_found ? -2 : millisec;
}


/*
int trace_callback_single_thread(st_uint_t ntracedtotal, st_uint_t ntraced,
											st_uint_t ntotrace, st_uint_t curstage, st_uint_t nstages,
											void *data)
{
	QProgressDialog *dlg = (QProgressDialog*)data;
	dlg->setRange( 0, ntotrace );
	dlg->setValue( ntraced );
	dlg->setLabelText( QString("Tracing Stage %1 of %2").arg(curstage).arg(nstages) );

	if (dlg->wasCanceled()) return 0;
	else return 1;
}



int RunTraceSingleThread( StSystem *System, int nrays, int nmaxrays,
						  int *seed, bool sunshape, bool opterrs,
						  QWidget *dialog_parent, QStringList &errors )
{
	ValidateSeed( seed );
	dialog_parent->setEnabled( false );

	QProgressDialog progress("Tracing...", "Abort", 0, nrays, dialog_parent);
	progress.setWindowModality(Qt::WindowModal);

	st_context_t spcxt = ::st_create_context();

	int result = LoadSystemIntoContext( System, spcxt, errors );
	if (result < 0)
	{
		::st_free_context(spcxt);
		errors.append( QObject::tr("Error loading system into simulation context") );
		dialog_parent->setEnabled(true);
		return -1;
	}

	::st_sim_errors( spcxt, sunshape?1:0, opterrs?1:0 );
	::st_sim_params( spcxt, nrays, nmaxrays );

	QTime sw;
	sw.start();

	int code = ::st_sim_run( spcxt, (unsigned int) *seed, trace_callback_single_thread, &progress );
	int millisec = (int)sw.elapsed();

	if (code < 0)
	{
		for ( int j=0;j<st_num_messages(spcxt);j++ )
			errors.append( st_message(spcxt, j ) );

		errors.append( QObject::tr("Error in trace, code %2").arg(code));

		System->Results.freeMemory();
	}
	else
	{
		if (!System->Results.readResultsFromContext(spcxt))
			errors.append(QObject::tr("Allocation error reading results from trace context - (mt)"));

		CountRayHitsPerElement( System );
	}


	// recompute transforms for results viewing
	System->recomputeTransforms();

	dialog_parent->setEnabled( true );

	return (code<0) ? -2 : millisec;
}
*/
