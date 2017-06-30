#ifndef __elementlist_h
#define __elementlist_h

#include <unordered_map>
#include <wx/scrolwin.h>

#include "project.h"


struct ElementItem
{
	unsigned int stage, element;

	ElementItem() { stage=element=0; }
	ElementItem( unsigned int s, unsigned int e ) { stage=s; element=e; }

	wxString Label( Project &prj );
};

class ElementListBox : public wxScrolledWindow
{
	Project &m_prj;
	int m_scrollRate;
	int m_space;
	int m_itemHeight;
	wxSize m_chkBoxSize, m_bestSize;
	std::vector<ElementItem> m_items;
	std::unordered_map<unsigned int,bool> m_sel;
	
public:
	ElementListBox( Project &prj, wxWindow *parent, int id, 
		const wxPoint &pos=wxDefaultPosition, 
		const wxSize &size=wxDefaultSize);
	bool IsSelected( size_t idx );
	bool IsSelected( unsigned int s, unsigned int e );
	bool GetStageElementIndices( size_t i, unsigned int *stage, unsigned int *elem );
	unsigned int Code( unsigned int s, unsigned int e );
	void Select( unsigned int s, unsigned int e, bool tf = true );
	void SelectAll();
	void UnselectAllStage( unsigned int stage );
	void ClearSelections();
	void Clear();
	void Reserve( size_t n );
	void Add( unsigned int stage, unsigned int element );
	size_t Count();

	void RecalculateBestSize();
	void OnPaint( wxPaintEvent &evt );
	void OnLeftDown(wxMouseEvent &evt);
	void OnErase(wxEraseEvent &evt);
	void OnResize(wxSizeEvent &evt);		
	void Invalidate();
	void ResetScrollbars();
	virtual wxSize DoGetBestSize() const;

	DECLARE_EVENT_TABLE();
};
#endif
