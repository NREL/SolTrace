#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <string>
#include <cstring>
#include <math.h>
#include <vector>
#include <exception>


#ifndef M_PI
	#define M_PI 3.141592653589793238462643
#endif

std::vector< std::string > split( const std::string &str, const std::string &delim, bool ret_empty, bool ret_delim )
{
	std::vector< std::string > list;

	char cur_delim[2] = {0,0};
	std::string::size_type m_pos = 0;
	std::string token;
	
	while (m_pos < str.length())
	{
		std::string::size_type pos = str.find_first_of(delim, m_pos);
		if (pos == std::string::npos)
		{
			cur_delim[0] = 0;
			token.assign(str, m_pos, std::string::npos);
			m_pos = str.length();
		}
		else
		{
			cur_delim[0] = str[pos];
			std::string::size_type len = pos - m_pos;			
			token.assign(str, m_pos, len);
			m_pos = pos + 1;
		}
		
		if (token.empty() && !ret_empty)
			continue;

		list.push_back( token );
		
		if ( ret_delim && cur_delim[0] != 0 && m_pos < str.length() )
			list.push_back( std::string( cur_delim ) );
	}
	
	return list;
}


#include "../stapi.h"

static void read_line(char *buf, int len, FILE *fp )
{
	fgets(buf, len, fp);
	int nch = strlen(buf);
	if (nch > 0 && buf[nch-1] == '\n')
		buf[nch-1] = 0;
	if (nch-1 > 0 && buf[nch-2] == '\r')
		buf[nch-2] = 0;
}

static void ZeroTransform(double ref[3][3], double loc[3][3], double eu[3])
{
	for (int i=0;i<3;i++)
	{
		eu[i]=0;
		for (int j=0;j<3;j++)
			ref[i][j]=loc[i][j]=0;
	}
}

static bool read_sun( FILE *fp, st_context_t cxt )
{
	if (!fp) return false;

	char buf[1024];
	int bi = 0, count = 0;
	char cshape = 'g';
	double Sigma, HalfWidth;
	bool PointSource;

	read_line( buf, 1023, fp );

	sscanf(buf, "SUN\tPTSRC\t%d\tSHAPE\t%c\tSIGMA\t%lg\tHALFWIDTH\t%lg",
		&bi, &cshape, &Sigma, &HalfWidth);
	PointSource = (bi!=0);
	cshape = tolower(cshape);
	
	st_sun( cxt, PointSource?1:0, cshape, cshape=='g' ? Sigma : HalfWidth );

	read_line( buf, 1023, fp );
	double X, Y, Z, Latitude, Day, Hour;
	bool UseLDHSpec;
	sscanf(buf, "XYZ\t%lg\t%lg\t%lg\tUSELDH\t%d\tLDH\t%lg\t%lg\t%lg",
		&X, &Y, &Z, &bi, &Latitude, &Day, &Hour);
	UseLDHSpec = (bi!=0);
		
	if ( UseLDHSpec )
	{
		double Declination, HourAngle, Elevation, Azimuth;

		Declination = 180/M_PI*asin(0.39795*cos(0.98563*M_PI/180*(Day-173)));
		HourAngle = 15*(Hour-12);
		Elevation = 180/M_PI*asin(sin(Declination*M_PI/180)*sin(Latitude*M_PI/180)+cos(Declination*M_PI/180)*cos(HourAngle*M_PI/180)*cos(Latitude*M_PI/180));
		Azimuth = 180/M_PI*acos((sin(M_PI/180*Declination)*cos(M_PI/180*Latitude)-cos(M_PI/180*Declination)*sin(M_PI/180*Latitude)*cos(M_PI/180*HourAngle))/cos(M_PI/180*Elevation)+0.0000000001);
		if ( sin(HourAngle*M_PI/180) > 0.0 )
			Azimuth = 360 - Azimuth;
		X = -sin(Azimuth*M_PI/180)*cos(Elevation*M_PI/180);
		Y = sin(Elevation*M_PI/180);
		Z = cos(Azimuth*M_PI/180)*cos(Elevation*M_PI/180);
	}

	st_sun_xyz( cxt, X, Y, Z );

	printf("sun ps? %d cs: %c  %lg %lg %lg\n", PointSource?1:0, cshape, X, Y, Z);

	read_line( buf, 1023, fp );
	sscanf(buf, "USER SHAPE DATA\t%d", &count);
	if (count > 0)
	{
		double *angle = new double[count];
		double *intensity = new double[count];

		for (int i=0;i<count;i++)
		{
			double x, y;
			read_line( buf, 1023, fp );
			sscanf(buf, "%lg\t%lg", &x, &y);
			angle[i] = x;
			intensity[i] = y;
		}

		st_sun_userdata(cxt, count, angle, intensity );

		delete [] angle;
		delete [] intensity;	
	}

	return true;
}

bool read_optic_surface(FILE *fp, st_context_t cxt, int iopt, int fb)
{
	if (!fp) return false;
	char buf[1024];
	read_line(buf, 1023, fp);
	std::vector<std::string> parts  = split( std::string(buf), "\t", true, false );
	if (parts.size() < 15)
	{
		printf("too few tokens for optical surface: %d\n", parts.size());
		printf("\t>> %s\n", buf);
		return false;
	}

	char ErrorDistribution = 'g';
	if (parts[1].length() > 0)
		ErrorDistribution = parts[1][0];

	int ApertureStopOrGratingType = atoi( parts[2].c_str() );
	int OpticalSurfaceNumber = atoi( parts[3].c_str() );
	int DiffractionOrder = atoi( parts[4].c_str() );
	double Reflectivity = atof( parts[5].c_str() );
	double Transmissivity = atof( parts[6].c_str() );
	double RMSSlope = atof( parts[7].c_str() );
	double RMSSpecularity = atof( parts[8].c_str() );
	double RefractionIndexReal = atof( parts[9].c_str() );
	double RefractionIndexImag = atof( parts[10].c_str() );
	double GratingCoeffs[4];
	GratingCoeffs[0] = atof( parts[11].c_str() );
	GratingCoeffs[1] = atof( parts[12].c_str() );
	GratingCoeffs[2] = atof( parts[13].c_str() );
	GratingCoeffs[3] = atof( parts[14].c_str() );

	bool UseReflectivityTable = false;
	int refl_npoints = 0;
	double *refl_angles = 0;
	double *refls = 0;

	bool UseTransmissivityTable = false;
	int trans_npoints = 0;
	double* trans_angles = 0;
	double* transs = 0;

	if (parts.size() >= 17)
	{
		UseReflectivityTable = (atoi( parts[15].c_str() ) > 0);
		refl_npoints = atoi( parts[16].c_str() );
		if (parts.size() >= 19)
		{
			UseTransmissivityTable = (atoi(parts[17].c_str()) > 0);
			trans_npoints = atoi(parts[18].c_str());
		}
	}

	if (UseReflectivityTable)
	{
		refl_angles = new double[refl_npoints];
		refls = new double[refl_npoints];

		for (int i=0;i<refl_npoints;i++)
		{
			read_line(buf,1023,fp);
			sscanf(buf, "%lg %lg", &refl_angles[i], &refls[i]);
		}
	}
	if (UseTransmissivityTable)
	{
		trans_angles = new double[trans_npoints];
		transs = new double[trans_npoints];

		for (int i = 0; i < trans_npoints; i++)
		{
			read_line(buf, 1023, fp);
			sscanf(buf, "%lg %lg", &trans_angles[i], &transs[i]);
		}
	}

	
	st_optic( cxt, iopt, fb, ErrorDistribution,
		OpticalSurfaceNumber, ApertureStopOrGratingType, DiffractionOrder,
		RefractionIndexReal, RefractionIndexImag,
		Reflectivity, Transmissivity,
		GratingCoeffs, RMSSlope, RMSSpecularity,
		UseReflectivityTable ? 1 : 0, refl_npoints,
		refl_angles, refls,
		UseTransmissivityTable? 1 : 0, trans_npoints,
		trans_angles, transs
		);

	if (refl_angles != 0) delete [] refl_angles;
	if (refls != 0) delete [] refls;
	if (trans_angles != 0) delete[] trans_angles;
	if (transs != 0) delete[] transs;
	return true;
}


bool read_optic(FILE *fp, st_context_t cxt)
{
	if (!fp) return false;
	char buf[1024];
	read_line( buf, 1023, fp );

	if (strncmp( buf, "OPTICAL PAIR", 12) == 0)
	{
		int iopt = st_add_optic( cxt, (const char*)(buf+13) );
		read_optic_surface( fp, cxt, iopt, 1 );
		read_optic_surface( fp, cxt, iopt, 2 );
		return true;
	}
	else return false;
}


bool read_element(FILE *fp, st_context_t cxt, int istage)
{
	int ielm = ::st_add_element( cxt, istage );

	char buf[1024];
	read_line(buf, 1023, fp);

	std::vector<std::string> tok = split( buf, "\t", true, false );
	if (tok.size() < 29)
	{
		printf("too few tokens for element: %d\n", tok.size());
		printf("\t>> %s\n", buf);
		return false;
	}

	
	st_element_enabled( cxt, istage, ielm,  atoi( tok[0].c_str() ) ? 1 : 0 );
	st_element_xyz( cxt, istage, ielm,  
		atof( tok[1].c_str() ), 
		atof( tok[2].c_str() ), 
		atof( tok[3].c_str() ) );
	st_element_aim( cxt, istage, ielm,  
		atof( tok[4].c_str() ), 
		atof( tok[5].c_str() ), 
		atof( tok[6].c_str() ) );
	st_element_zrot( cxt, istage, ielm,  atof( tok[7].c_str() ) );
	if (tok[8].length() > 0) st_element_aperture( cxt, istage, ielm, tok[8][0] );
	
	double Params[8];

	for (int i=0;i<8;i++)
		Params[i] = atof( tok[i+9].c_str() );

	st_element_aperture_params( cxt, istage, ielm,  Params );



	if (tok[17].length() > 0) st_element_surface( cxt, istage, ielm, tok[17][0] );
	
	for (int i=0;i<8;i++)
		Params[i] = atof( tok[i+18].c_str() );

	std::string SurfaceFile = tok[26];
	if (!SurfaceFile.empty())
	{
		if ( st_element_surface_file( cxt, istage, ielm,  SurfaceFile.c_str() ) < 0)
		{
			for (int j=0;j<st_num_messages(cxt);j++)
				printf( "error: %s\n", st_message(cxt, j) );
			return false;
		}
	}
	else
		st_element_surface_params( cxt, istage, ielm, Params );

	
	st_element_optic( cxt, istage, ielm,  tok[27].c_str() );
	st_element_interaction( cxt, istage, ielm,  atoi( tok[28].c_str()) );

	return true;
}


bool read_stage(FILE *fp, st_context_t cxt)
{
	if (!fp) return false;

	char buf[1024];
	read_line( buf, 1023, fp );

	int virt=0,multi=1,count=0,tr=0;
	double X, Y, Z, AX, AY, AZ, ZRot;


	sscanf(buf, "STAGE\tXYZ\t%lg\t%lg\t%lg\tAIM\t%lg\t%lg\t%lg\tZROT\t%lg\tVIRTUAL\t%d\tMULTIHIT\t%d\tELEMENTS\t%d\tTRACETHROUGH\t%d",
		&X, &Y, &Z,
		&AX, &AY, &AZ,
		&ZRot,
		&virt,
		&multi,
		&count,
		&tr );

	read_line( buf, 1023, fp ); // read name

	int istage = st_add_stage( cxt );
	
	::st_stage_flags( cxt, istage, virt, multi, tr );

	st_stage_xyz(cxt, istage, X, Y, Z );
	st_stage_aim(cxt, istage, AX, AY, AZ );
	st_stage_zrot(cxt, istage, ZRot );

	printf("stage '%s': [%d] %lg %lg %lg   %lg %lg %lg   %lg   %d %d %d\n",
		buf, count, X, Y, Z, AX, AY, AZ, ZRot, virt, multi, tr );

	st_clear_elements(cxt, istage);
	
	for (int i=0;i<count;i++)
		if (!read_element( fp, cxt, istage )) 
		{ printf("error in element %d of %d in stage %d\n", i, count, istage ); return false; }

	return true;
}

bool read_system(FILE *fp, st_context_t cxt)
{
	if (!fp) return false;

	char buf[1024];

	char c = fgetc(fp);
	if ( c == '#' )
	{
		int vmaj = 0, vmin = 0, vmic = 0;
		read_line( buf, 1023, fp ); sscanf( buf, " SOLTRACE VERSION %d.%d.%d INPUT FILE", &vmaj, &vmin, &vmic);

		unsigned int file_version = vmaj*10000 + vmin*100 + vmic;
		
		printf( "loading input file version %d.%d.%d\n", vmaj, vmin, vmic );
	}
	else
	{
		ungetc( c, fp );
		printf("input file must start with '#'\n");
		return false;
	}

	if ( !read_sun( fp, cxt ) ) return false;
	
	int count = 0;

	count = 0;
	read_line( buf, 1023, fp ); sscanf(buf, "OPTICS LIST COUNT\t%d", &count);
	
	::st_clear_optics( cxt );
	for (int i=0;i<count;i++)
		if (!read_optic( fp, cxt )) return false;

	count = 0;
	read_line( buf, 1023, fp ); sscanf(buf, "STAGE LIST COUNT\t%d", &count);
	::st_clear_stages( cxt );
	for (int i=0;i<count;i++)
		if (!read_stage( fp, cxt )) return false;

	return true;
}

bool write_data_file(const char *file, st_context_t cxt)
{
	std::ofstream fout;
	try
	{
		fout.open(file);
	}
	catch (...)
	{
		return false;
	}

	double buf[9];
	double SunXMin, SunXMax, SunYMin, SunYMax;
	int SunRayCount;
	::st_sun_stats( cxt, &SunXMin, &SunXMax, &SunYMin, &SunYMax, &SunRayCount );

	int Length = ::st_num_intersections( cxt );
	
	double *Xi = new double[Length];
	double *Yi = new double[Length];
	double *Zi = new double[Length];
	double *Xc = new double[Length];
	double *Yc = new double[Length];
	double *Zc = new double[Length];
	int *Em = new int[Length];
	int *Sm = new int[Length];
	int *Rn = new int[Length];

	::st_locations( cxt, Xi, Yi, Zi );
	::st_cosines( cxt, Xc, Yc, Zc );
	::st_elementmap( cxt, Em );
	::st_stagemap( cxt, Sm );
	::st_raynumbers( cxt, Rn );
	
	buf[0] = SunXMin;
	buf[1] = SunXMax;
	buf[2] = SunYMin;
	buf[3] = SunYMax;
	buf[4] = SunRayCount;
	buf[5] = Length;

	fout << "#       SunXMin        SunXMax        SunYMin        SunYMax    SunRayCount         Length\n";
	for (int i = 0; i < 6; i++)
		fout << std::scientific << std::setw(15) << std::setprecision(6) << buf[i];
	fout << "\n";
	
	fout << "#             X              Y              Z           Xcos           Ycos           Zcos       Stagemap     Elementmap      Raynumber\n";
	for (size_t i=0;i<Length;i++)
	{
		buf[0] = Xi[i];
		buf[1] = Yi[i];
		buf[2] = Zi[i];
		buf[3] = Xc[i];
		buf[4] = Yc[i];
		buf[5] = Zc[i];
		buf[6] = Sm[i];
		buf[7] = Em[i];
		buf[8] = Rn[i];

		//fwrite( (void*)buf, sizeof(double), sizeof(buf), fp );
		for (int j = 0; j < 6; j++)
			fout << std::scientific << std::setw(15) << std::setprecision(6) << buf[j];
		for (int j = 6; j < 9; j++)
			fout << std::fixed << std::setw(15) << std::setprecision(0) << (int)buf[j];
		fout << "\n";
	}
	
	delete [] Xi;
	delete [] Yi;
	delete [] Zi;
	
	delete [] Xc;
	delete [] Yc;
	delete [] Zc;

	delete [] Em;
	delete [] Sm;
	delete [] Rn;

	//fclose(fp);
	fout.close();

	return true;
}

int trace_progress(st_uint_t ntracedtotal, st_uint_t ntraced,
											st_uint_t ntotrace, st_uint_t curstage, st_uint_t nstages,
											void *data)
{
	static int last_stage = -1;
	static double last_percent = -1;

	if (curstage != last_stage)
	{
		printf("Stage %2d:\n", (int)curstage);
		last_stage = curstage;
	}
	
	return 1;
}

int main(int argc, char *argv[])
{
	// strace.exe <system.stinput> <number of rays> <max rays>  <seed> <sunshape? 0/1> <errors? 0/1>
	// example strace.exe supernova.stinput 1000000 100000000 8454 1 1
	if (argc < 7)
	{
		printf("strace: too few arguments.  usage:\n\t"
			"strace.exe <system.stinput> <number of rays> <max rays> <seed> <sunshape? 0/1> <errors? 0/1> <pointfocus? 0/1>\n");
		return -1;
	}

	//printf("Hit enter to continue...\n");
	//getchar();
	
	const char *file = argv[1];
	int nrays = atoi( argv[2] );
	int maxrays = atoi( argv[3] );
	int seed = atoi( argv[4] );
	int sunshape = atoi( argv[5] );
	int errors = atoi( argv[6] );
    int aspointfocus = atoi( argv[7] );

	st_context_t cxt = ::st_create_context();

	FILE *fp = fopen(file, "r");
	if (!fp)
	{
		printf("failed to open system input file\n");
		return -1;
	}

	printf("input file: %s\n", file);
	if ( !read_system( fp, cxt ) )
	{
		printf("error in input file.\n");
		fclose(fp);
		return -1;
	}
	
	fclose(fp);

	::st_sim_params( cxt, nrays, maxrays );
	::st_sim_errors( cxt, sunshape, errors );
	int code = ::st_sim_run( cxt, (unsigned int)seed, aspointfocus==aspointfocus, trace_progress, 0 );
	
	if ( code >= 0 )
	{
		printf("trace finished with code %d\n", code );
		std::string output(file);
		output += ".rays";
		write_data_file( output.c_str(), cxt );
	}
	else
	{
		for (int j=0;j<st_num_messages(cxt);j++)
			printf( "error: %s\n", st_message(cxt, j) );
	}
	
	::st_free_context( cxt );
	return 0;
}