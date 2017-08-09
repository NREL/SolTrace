#include "soltrace.h"

#ifdef ST_CONSOLE_APP

#define _CONSOLE

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
    #pragma hdrstop
#endif

// for all others, include the necessary headers (this file is usually all you
// need because it includes almost all "standard" wxWidgets headers)
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include <wx/app.h>
#include <wx/cmdline.h>

#include "project.h"
#include "trace.h"

// ============================================================================
// implementation
// ============================================================================

static const wxCmdLineEntryDesc cmdLineDesc[] =
{
    { wxCMD_LINE_SWITCH, "h", "help", "show this help message", wxCMD_LINE_VAL_NONE, wxCMD_LINE_OPTION_HELP },
    { wxCMD_LINE_OPTION, "f", "file", "Run an stinput file", wxCMD_LINE_VAL_STRING, wxCMD_LINE_OPTION_MANDATORY },
    { wxCMD_LINE_OPTION, "r", "rays", "Required ray intersections (=1e4)", wxCMD_LINE_VAL_NUMBER },
    { wxCMD_LINE_OPTION, "m", "maxrays", "Max. number of sun rays (=100*rays)", wxCMD_LINE_VAL_NUMBER },
    { wxCMD_LINE_OPTION, "c", "cpu", "Number of CPU's (=16)", wxCMD_LINE_VAL_NUMBER},
    { wxCMD_LINE_OPTION, "d", "seed", "Random seed (=-1)", wxCMD_LINE_VAL_NUMBER},
    { wxCMD_LINE_OPTION, "p", "sunshape", "Enable sunshape (=1)", wxCMD_LINE_VAL_NUMBER},
    { wxCMD_LINE_OPTION, "e", "error", "Enable optical error (=1)", wxCMD_LINE_VAL_NUMBER},
    { wxCMD_LINE_OPTION, "t", "tower", "Run as power tower (=1)", wxCMD_LINE_VAL_NUMBER},
    { wxCMD_LINE_OPTION, "o", "out", "File to write output (=trace.out)", wxCMD_LINE_VAL_STRING},
    //{ wxCMD_LINE_OPTION, "s", "script", "Run an lk script file", wxCMD_LINE_VAL_STRING},
    //{ wxCMD_LINE_SWITCH, "i", "interactive", "Interactive mode"},
    { wxCMD_LINE_NONE }
};

int main(int argc, char **argv)
{
    wxApp::CheckBuildOptions(WX_BUILD_OPTIONS_SIGNATURE, "program");

    wxInitializer initializer;
    if ( !initializer )
    {
        fprintf(stderr, "Failed to initialize the wxWidgets library, aborting.");
        return -1;
    }
    
    //default values
    wxString fname = "C:/Users/mwagner/Documents/NREL/projects/SolTrace-git/app/deploy/x64/A1.stinput";
    wxString fnout = "trace.out";
    long rays = (int)1e4;
    long maxrays = 100*rays;
    long threads = 16;
    long l_seed = -1;
    long l_sunshape = 1;
    long l_error = 1;
    long l_tower = 1;

    wxCmdLineParser parser(cmdLineDesc, argc, argv);

    parser.Parse();
    
    if( parser.Found("f", &fname) )
    {
        //wxPrintf( fname.mb_str() );

        //check that file exists
        if(! ::wxFileExists( fname ) )
        {
            wxPrintf( "\nInput file not found! Invalid path." );
            return 0;
        }
    }

    parser.Found("r", &rays);
    
    parser.Found("m", &maxrays);
    
    parser.Found("c", &threads);

    parser.Found("d", &l_seed);

    parser.Found("p", &l_sunshape);

    parser.Found("e", &l_error);

    parser.Found("t", &l_tower);
            
    if( parser.Found("o", &fnout) )
    {
        //check that file exists
        if(! ::wxFileExists( fnout ) )
        {
            wxPrintf( "\nOutput file not found! Invalid path." );
            return 0;
        }
    }


    if ( argc == 1 )
    {
        // If there were no command-line options supplied, emit a message
        // otherwise it's not obvious that the sample ran successfully
        wxPrintf("Welcome to the SolTrace 'console' interface!\n");
        wxPrintf("For more information, run again with the --help option\n");
        return 0;
    }

    // create and execute according to the commands
    FILE *fp_in = fopen( fname.c_str(), "r" );
	
	Project project;
    project.Read( fp_in );
    
    fclose(fp_in); 
    
    //wxPrintf( "\nAperture 0: %f", project.StageList.at(0)->ElementList.at(0)->ApertureParams[0] );

    //type conversion
    int seed = (int)l_seed;
    bool sunshape = l_sunshape == 1L;
    bool error = l_error == 1L;
    bool tower = l_tower == 1L;

    wxArrayString ref_errors;

    int msec = RunTraceMultiThreaded( &project, 
            (int)rays,
			(int)maxrays,
			threads,
			&seed,
			sunshape,
			error,
            tower,
			ref_errors,
            true);


    //------------------------------------------------------------------------------
    FILE *fp_out = fopen( fnout.c_str(), "w" );
	if ( !fp_out )
	{
		wxPrintf("Could not open file for writing:\n\n" + fnout);
		return 0;
	}

    //create output
    if ( msec < 0 )
		wxPrintf( wxJoin( ref_errors, '\n' ) );

    // by default, all ray data is in stage coordinates
	RayData &rd = project.Results;
	
	int nwrite = rd.PrepareExport( RayData::COORD_GLOBAL, 0 );
	
	// progress update every 0.5 %
	int nupdate = nwrite / 200;
	size_t nwr = 0;
	size_t bytes = 0;
	
	double Pos[3], Cos[3];
	int Elm, Stg, Ray;

	fputs( "Pos X,Pos Y,Pos Z,Cos X,Cos Y,Cos Z,Element,Stage,Ray Number\n", fp_out );
	while( rd.GetNextExport( project, Pos, Cos, Elm, Stg, Ray ) )
	{
		/*if ( nwr % nupdate == 0 )
		{
			bool proceed = pd.Update( (int)( 100*((double)nwr)/((double)nwrite) ), 
					"Writing data to " + dlg.GetPath() + wxString::Format(" (%.2lf MB)", bytes*0.000001 ) );
			if( !proceed )
				break;
		}*/

		int nb = fprintf( fp_out, "%lg,%lg,%lg,%lg,%lg,%lg,%d,%d,%d\n",
			Pos[0], Pos[1], Pos[2],
			Cos[0], Cos[1], Cos[2],
			Elm, Stg, Ray );

		if ( nb < 0 ) break; // file error, disk space issue?

		bytes += nb;
		nwr++;
	}

	if ( ferror( fp_out ) )
		wxPrintf(wxString("An error occurred exporting the CSV data file: ") + strerror(errno));
    
    fclose(fp_out);

    return 0;
}



#endif