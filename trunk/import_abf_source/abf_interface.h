#ifndef _abf_INTERFACE_h_
#define _abf_INTERFACE_h_

/* import data from an Axon Instruments binary (.abf) file */
/* header for static C code to wrap Axon's C++ DLL */
/* JAB 6/15/07 jbender@caltech.edu */

/* use Matlab's memory functions if running in a MEX environment */
#ifdef _abf_USING_MEX_
#define malloc( x ) mxMalloc( x )
#define free( x ) mxFree( x )
#endif

#define _ABF_READ_SUCCESS 0
#define _ABF_OPEN_ERROR 1
#define _ABF_CHANNEL_ERROR 2
#define _ABF_EPISODE_ERROR 3
#define _ABF_N_SAMPLES_ERROR 4
#define _ABF_READ_ERROR 5
#define _ABF_UNSUPPORTED_ERROR 6
#define _ABF_OFFSET_EPISODIC_ERROR 7

/* returns the number of channels, episodes, and samples in this file
   for gap-free data, send N_EPISODES = -1 */
int get_filesize( char *filename, long *n_channels, long *n_episodes, long *n_samples );

/* reads SAMPLES beginning at OFFSET from one CHANNEL in FILENAME, stores in DATA
   CHANNEL is 0-based and logical, not physical
   use EPISODE = -1 for gap-free recordings, otherwise reads the specified episode
   this function does no error checking, so use get_filesize() and check your
     indices before calling */
int read_channel(
#ifdef _abf_DEBUG_
       FILE *logfile,
#endif
       char *filename, int channel, long offset, long samples,
       int episode, void *data, float *sampleInterval );

/* internal function for clean exiting */
int cleanup( void *memory, int fp, int return_code );

#endif
