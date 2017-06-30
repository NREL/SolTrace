#ifndef __trace_h
#define __trace_h

#include <wx/wx.h>
#include <wx/panel.h>

#include "project.h"

class wxExtTextCtrl;
class wxNumericCtrl;
class wxCheckBox;

int RunTraceMultiThreaded( Project *System, int nrays, int nmaxrays,
						int nmaxthreads, int *seed, bool sunshape, bool opterrs, bool aspowertower,
						wxArrayString &errors );

class TraceForm : public wxPanel
{
	Project &m_prj;

public:
	TraceForm( wxWindow *parent, Project &m_prj );

	void SetOptions( size_t nrays, size_t nmaxsunrays, int ncpu, int seed,
		bool sunshape, bool opterr, bool aspowertower );
	void GetOptions( size_t *nrays, size_t *nmaxsunrays, int *ncpu, int *seed,
		bool *sunshape, bool *opterr, bool *aspowertower );

	void SetWorkDir( const wxString &path );
	wxString GetWorkDir();
	
	bool IsRunning();
	// returns milliseconds elapsed, or negative number indicating error
	int StartTrace( bool mt = true, bool quiet = false, wxArrayString *err = 0 );
	void CancelTrace();
	int GetLastSeedVal() { return m_lastSeedVal; }

private:
	
	void OnCommand( wxCommandEvent &evt );

	wxNumericCtrl *m_numRays, *m_numMaxSunRays, *m_numCpus, *m_seed;
	wxCheckBox *m_inclSunShape, *m_inclOpticalErrors, *m_asPowerTower;
	wxExtTextCtrl *m_workDir;

	int m_lastSeedVal;
	wxNumericCtrl *m_elapsedTime, *m_lastSeed;

	DECLARE_EVENT_TABLE();
};


#endif
