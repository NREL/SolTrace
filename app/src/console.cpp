#include "soltrace.h"

#ifdef ST_CONSOLE_APP

#define _CRT_SECURE_NO_DEPRECATE 1
#define _CRT_NON_CONFORMING_SWPRINTFS 1
#define _SCL_SECURE_NO_WARNINGS 1
//#define NDEBUG
#define _CONSOLE
#define wxUSE_GUI 0

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

// ============================================================================
// implementation
// ============================================================================

static const wxCmdLineEntryDesc cmdLineDesc[] =
{
    { wxCMD_LINE_SWITCH, "h", "help", "show this help message",
        wxCMD_LINE_VAL_NONE, wxCMD_LINE_OPTION_HELP },
    { wxCMD_LINE_SWITCH, "d", "dummy", "a dummy switch" },
    // ... your other command line options here...

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

    wxCmdLineParser parser(cmdLineDesc, argc, argv);
    switch ( parser.Parse() )
    {
        case -1:
            // help was given, terminating
            break;

        case 0:
            // everything is ok; proceed
            if (parser.Found("d"))
            {
                wxPrintf("Dummy switch was given...\n");

                while (1)
                {
                    wxChar input[128];
                    wxPrintf("Try to guess the magic number (type 'quit' to escape): ");
                    if ( !wxFgets(input, WXSIZEOF(input), stdin) )
                        break;

                    // kill the last '\n'
                    input[wxStrlen(input) - 1] = 0;

                    if (wxStrcmp(input, "quit") == 0)
                        break;

                    long val;
                    if (!wxString(input).ToLong(&val))
                    {
                        wxPrintf("Invalid number...\n");
                        continue;
                    }

                    if (val == 42)
                        wxPrintf("You guessed!\n");
                    else
                        wxPrintf("Bad luck!\n");
                }
            }
            break;

        default:
            break;
    }

    if ( argc == 1 )
    {
        // If there were no command-line options supplied, emit a message
        // otherwise it's not obvious that the sample ran successfully
        wxPrintf("Welcome to the wxWidgets 'console' sample!\n");
        wxPrintf("For more information, run it again with the --help option\n");
    }

    // do something useful here

    return 0;
}



#endif