#ifndef __script_h
#define __script_h

#include <vector>
#include <wx/frame.h>
#include <wx/stc/stc.h>

#include <wex/lkscript.h>

class SolTraceScriptWindowFactory : public wxLKScriptWindowFactory
{
public:
	SolTraceScriptWindowFactory();
	virtual ~SolTraceScriptWindowFactory();
	virtual wxLKScriptWindow *Create();
};


class SolTraceScriptWindow  : public wxLKScriptWindow
{
public:
	SolTraceScriptWindow( wxWindow *parent, int id = wxID_ANY );
	
protected:
	virtual void OnScriptStarted();
	virtual void OnScriptStopped();
	virtual void OnHelp();

	DECLARE_EVENT_TABLE();
};


#endif

