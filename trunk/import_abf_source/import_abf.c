#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mex.h"
#include "abf_interface.h"

/* import data from an Axon Instruments binary (.abf) file */
/* see import_abf.m for more information */
/* JAB 6/15/07 jbender@caltech.edu */

#define FUNC_NAME		"import_abf"	/* Function name string */
#define ERR_MSG_LEN	512			/* char array size for an output error message */
#define ARG_BUF_SIZE	512			/* buffer size for input argument (filename) */

/* Get input argument */
void get_input_arg_string( const mxArray *arg_0, char *string );
int get_input_arg_int( const mxArray *arg_0 );

/* give informational error strings when crashing */
void parse_error_code( void *data, const int ec );

/* Mex function ************************************************************/
void mexFunction(
	 int nlhs,		/* # of left-hand side arguments */ 
	 mxArray *plhs[],	/* Left-hand side arguments */ 
	 int nrhs,		/* # of right-hand side arguments */
	 const mxArray *prhs[]	/* Right-hand side arguments */
)
{
	char filename[ARG_BUF_SIZE];
	char err_msg[ERR_MSG_LEN];
        void *data = mxMalloc( 64 ); /* malloc so it can be freed upon exit */
        //float *samplePeriod = (float*)mxMalloc( sizeof( float ) );
	float *samplePeriod;
    int *numEpisodes;
	int rows, cols;
	int episode=-1;
    long offset=-1, samples=-1;
	int error_code = 0;
    long n_chan, n_epi, n_samp;
    int i;
#ifdef _abf_DEBUG_
    FILE *logfile = fopen( "c:\logfile.txt", "w" );

    if( logfile == 0 ) mexErrMsgTxt( "error opening logfile\n" );
#endif

	/* check number of input args */
	if( nrhs < 1 || nrhs > 3 ) {
		sprintf( err_msg, "%s: wrong number of input arguments ", FUNC_NAME );
		mexErrMsgTxt( err_msg );
	}
	/* check number of output arguments */
	if ( nlhs > 3 ) {
		sprintf( err_msg, "%s: too many output arguments.", FUNC_NAME );
		mexErrMsgTxt( err_msg );
	}

	/* extract input arguments */
	get_input_arg_string( prhs[0], filename );
  if( nrhs == 2 ) /* episode number */
    episode = get_input_arg_int( prhs[1] );
  if( nrhs == 3 ) { /* samples and offset */
    offset = get_input_arg_int( prhs[1] );
    samples = get_input_arg_int( prhs[2] );
  }

    /* compare arguments with file data */
    error_code = get_filesize( filename, &n_chan, &n_epi, &n_samp );
    parse_error_code( data, error_code );

    if( offset == -1 ) offset = 0;

 //   if( offset == -1 ) {
 //     if( episode != -1 ) /* episodic file */
 //       parse_error_code( data, _ABF_OFFSET_EPISODIC_ERROR );
 //     else offset = 0;
 //   } 
    if( episode > n_epi ) {
        sprintf( err_msg, "%s: file only has %d episodes", FUNC_NAME, n_epi );
        mxFree( data );
        mexErrMsgTxt( err_msg );
    }
    if( samples == -1 ) samples = n_samp; /* try to read all the samples */
    if( offset > n_samp ) {
        sprintf( err_msg, "%s: offset is greater than number of samples in file (%d)", FUNC_NAME, n_samp );
        mxFree( data );
        mexErrMsgTxt( err_msg );
    }
    if( (offset == -1 && samples > n_samp) || (samples + offset > n_samp) ) {
        sprintf( err_msg, "%s: file only has %d samples (per episode, if episodic)", FUNC_NAME, n_samp );
        mexWarnMsgTxt( err_msg );
        if( offset == -1 ) samples = n_samp;
        else samples = n_samp - offset;
    }
#ifdef _abf_DEBUG_
    fprintf( logfile, "n_chan %d n_epi %d n_samp %d\n", n_chan, n_epi, n_samp );
    fprintf( logfile, "reading epi %d, offset %d, samples %d\n", episode, offset, samples );
#endif

    /* allocate memory -- could fail/crash/etc. if user is asking for too many samples */
    rows = n_chan;
    cols = samples;
    mxFree( data ); /* allocated at init; free before reallocating */
    data = mxMalloc( sizeof( float ) * rows * cols );
    plhs[1] = mxCreateNumericMatrix( 1, 1, mxSINGLE_CLASS, 0 );
    plhs[2] = mxCreateNumericMatrix( 1, 1, mxINT16_CLASS, 0);
    samplePeriod = (float*) mxGetData( plhs[1] );
    numEpisodes = (int*) mxGetData( plhs[2] );
    *numEpisodes = n_epi;
#ifdef _abf_DEBUG_
    fprintf( logfile, "data allocated size %d at %d\n", sizeof( float )*rows*cols, data );
#endif

    /* read data from file, channel by channel */
    for( i = 0; i < n_chan; i++ ) {
#ifdef _abf_DEBUG_
      fprintf( logfile, "reading channel %d\n", i );
      error_code = read_channel( logfile, filename, i, offset, samples, episode, data, samplePeriod );
#else
      error_code = read_channel( filename, i, offset, samples, episode, data, samplePeriod );
#endif
      parse_error_code( data, error_code );
    }
#ifdef _abf_DEBUG_
    fprintf( logfile, "success -- returning\ndata: " );
    for( i = 0; i < 20; i++ )
      fprintf( logfile, "%d:%.2f ", i, *((float*)data + i) );
    fprintf( logfile, "\n" );
#endif

  //mexPrintf("sample rate: %f\n", *samplePeriod );


    /* format answer for return */
    if( nlhs > 0 ) {
      if( sizeof( float ) == 4 )
      {
        plhs[0] = mxCreateNumericMatrix( cols, rows, mxSINGLE_CLASS, 0 );
      } else if( sizeof( float ) == 8 )
      {
        plhs[0] = mxCreateNumericMatrix( cols, rows, mxDOUBLE_CLASS, 0 );
      } else {
        sprintf( err_msg, "%s: C type 'float' on this machine is neither 32- nor 64-bit", FUNC_NAME );
        mexErrMsgTxt( err_msg );
      }
      mxFree( mxGetPr( plhs[0] ) ); /* free memory allocated by creation */
      mxSetData( plhs[0], data ); /* replace with data */
      mxSetData( plhs[1],  samplePeriod );
      mxSetData( plhs[2],  numEpisodes );
    }
#ifdef _abf_DEBUG_
    fclose( logfile );
#endif
}

/*******************************************************************************
 *
 * Function get_input_arg_string
 *
 *******************************************************************************/
void get_input_arg_string( const mxArray *arg_0, char *string )
{
	int name_len;
	char err_msg[ERR_MSG_LEN];  /* Error message string*/
	
	/* Check argument type */
	if ( !mxIsChar( arg_0 ) ) {
		sprintf( err_msg, "%s: input argument of incorrect type.", FUNC_NAME );
		mexErrMsgTxt( err_msg  );
	}

	/* Check that argument has only one row */
	if ( mxGetM( arg_0 ) != 1 ) {
		sprintf( err_msg, "%s: input argument has too many rows.", FUNC_NAME ); 
		mexErrMsgTxt( err_msg );
	}

	/* Get length of filename string */
	name_len = mxGetN( arg_0 );

	/* Make sure command is less than max possible */
	if ( name_len > ARG_BUF_SIZE ) {
		sprintf( err_msg, "%s: filename longer than maximum allowed size %d chars.", FUNC_NAME, ARG_BUF_SIZE );
		mexErrMsgTxt( err_msg );
	}

	/* Turn argument into C style string */
	if( mxGetString( arg_0, string, name_len + 1 ) ) {
		sprintf( err_msg, "%s: unable to convert Matlab char array to C style string. ", FUNC_NAME );
		mexErrMsgTxt( err_msg );
	}
}

/*******************************************************************************
 *
 * Function get_input_arg_int
 *
 *******************************************************************************/
int get_input_arg_int( const mxArray *arg_0 )
{
    double argval;
	char err_msg[ERR_MSG_LEN];  /* Error message string*/
	
	/* Check argument type */
	if ( !mxIsNumeric( arg_0 ) ) {
		sprintf( err_msg, "%s: input argument of incorrect type.", FUNC_NAME );
		mexErrMsgTxt( err_msg  );
	}

	/* Check that argument has only one row */
	if ( mxGetM( arg_0 ) != 1 ) {
		sprintf( err_msg, "%s: input argument has too many rows.", FUNC_NAME ); 
		mexErrMsgTxt( err_msg );
	}

	/* Check that argument has only one column */
	if ( mxGetN( arg_0 ) != 1 ) {
		sprintf( err_msg, "%s: input argument has too many columns.", FUNC_NAME ); 
		mexErrMsgTxt( err_msg );
	}

	/* Turn argument into C style int */
	argval = mxGetScalar( arg_0 );
	return (int)argval;
}

/*******************************************************************************
 *
 * Function parse_error_code
 *
 *******************************************************************************/
void parse_error_code( void *data, const int ec )
{
	char err_msg[ERR_MSG_LEN];
	
    if( ec != _ABF_READ_SUCCESS ) {
        if( ec == _ABF_OPEN_ERROR )
          sprintf( err_msg, "%s: error opening file (file not found?)", FUNC_NAME );
        else if( ec == _ABF_CHANNEL_ERROR )
          sprintf( err_msg, "%s: bad channel number", FUNC_NAME );
        else if( ec == _ABF_EPISODE_ERROR )
          sprintf( err_msg, "%s: bad episode number", FUNC_NAME );
        else if( ec == _ABF_N_SAMPLES_ERROR )
          sprintf( err_msg, "%s: wrong number of samples read", FUNC_NAME );
        else if( ec == _ABF_READ_ERROR )
          sprintf( err_msg, "%s: read error", FUNC_NAME );
        else if( ec == _ABF_UNSUPPORTED_ERROR )
          sprintf( err_msg, "%s: file mode unsupported (only gap-free or episodic files are supported)", FUNC_NAME );
        else if( ec == _ABF_OFFSET_EPISODIC_ERROR )
          sprintf( err_msg, "%s: cannot read episodic data with an offset", FUNC_NAME );
        else sprintf( err_msg, "%s: unknown error %d", FUNC_NAME, ec );
        mxFree( data );
        mexErrMsgTxt( err_msg );
    }
}
