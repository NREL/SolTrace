#ifndef __geometry_h
#define __geometry_h

#include <wx/wx.h>
#include <wx/panel.h>
#include <wx/grid.h>

#include "project.h"

class wxNotebook;

class wxExtGridCtrl;
class wxExtTextCtrl;
class wxNumericCtrl;

class ElementTable;
class GeometryForm;

class StageForm : public wxPanel
{
public:
	StageForm( wxWindow *parent, GeometryForm *geo, Project &prj, Stage *s );

	void UpdateFromData();

	void Import(const wxString &file);  //TJW 6-29-17 adding capability to import file w/o the dialog, but by directly passing the file as argumen
	//void Import();                    //TJW 6-29-17
	void Export();
	void Clear();
	void Delete(int idx=-1);
	void Append(int n=-1);
	void Insert();

	void EditStageZRot();
	void EditZRot(int idx=-1);
	void EditAperture(int idx=-1);
	void EditSurface(int idx=-1);
	void EditOptics(int idx=-1);

	void Modified();
	void NameModified();

private:

	void UpdateProperties();
	Project &m_prj;
	Stage *m_stage;

	GeometryForm *m_geoForm;

	wxExtTextCtrl *m_stageName;
	wxCheckBox *m_virtualStage, *m_multipleHits, *m_traceThrough;
	wxNumericCtrl *m_x, *m_y, *m_z, *m_ax, *m_ay, *m_az, *m_zrot;

	ElementTable *m_gridTable;
	wxExtGridCtrl *m_grid;


	void OnGridCellChange(wxGridEvent &evt);
	void OnGridCellSelect(wxGridEvent &evt);
	void OnGridRangeSelect(wxGridRangeSelectEvent &evt);
	void OnGridCellDClick(wxGridEvent &evt);
	void OnGridCellRightClick(wxGridEvent &evt);
	
	void OnCommand( wxCommandEvent &evt );

	DECLARE_EVENT_TABLE()
};

class GeometryForm : public wxPanel
{
public:
	GeometryForm( wxWindow *parent, Project &prj );

	void NewStage(const wxString &name = "New stage", int pos = -1 );
	void DeleteStage(int sel, bool quiet=false);
	void ClearStages(bool quiet=false);
	StageForm *GetStageForm(int idx);
	StageForm *GetStageForm( Stage *stg );

	void UpdateForm();
	void UpdateStageNames();
	void UpdateStageFlagsAndGeom(int stage_num);
	
	void Modified();
private:
	Project &m_prj;
	wxNotebook *m_stageTabs;

	void OnButton(wxCommandEvent &evt);

	DECLARE_EVENT_TABLE()
};

#endif
