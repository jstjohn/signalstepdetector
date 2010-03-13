#include <stdio.h>
#include <stdlib.h>
#include "abf_interface.h"
#include "abffiles.h"
#ifdef _abf_USING_MEX_
#include "mex.h"
#endif

/* import data from an Axon Instruments binary (.abf) file */
/* source for static C code to wrap Axon's C++ DLL */
/* JAB 6/15/07 jbender@caltech.edu */

#define MIN( x, y ) (x > y? y : x)

/* returns the number of channels, episodes, and samples in this file
   for gap-free data, send N_EPISODES = -1 */
int get_filesize( char *filename, long *n_channels, long *n_episodes, long *n_samples )
{
  int fp; /* file pointer */
  int error_code = 0;
  ABFFileHeader header;
  DWORD max_episodes = 0; /* returns number of sweeps */
  UINT max_samples = 0; /* returns number of samples */

  /* try to open file */  
  if( !ABF_ReadOpen( filename, &fp, ABF_DATAFILE,
      &header, &max_samples, &max_episodes, &error_code ) )
    return _ABF_OPEN_ERROR;

  if( header.nOperationMode != 3 && /* gap-free */
      header.nOperationMode != 5 && /* episodic */
      header.nOperationMode != 2) { /* Fixed-Length */     
    return _ABF_UNSUPPORTED_ERROR;
  }

  *n_channels = header.nADCNumChannels;
  if( header.nOperationMode == 5 || header.nOperationMode == 2 ) { /* episodic */
    *n_samples = max_samples;
    *n_episodes = max_episodes;
  }
  else if( header.nOperationMode == 3 ) { /* gap-free */
    /* "episodes" are artificial; just return total samples */
    *n_episodes = -1;
    *n_samples = max_samples * max_episodes;
  }

  ABF_Close( fp, &error_code );
  return _ABF_READ_SUCCESS;
}

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
       int episode, void *data, float *sampleInterval )
{
  int fp; /* file pointer */
  int error_code = 0;
  ABFFileHeader header;
  DWORD max_episodes = 0; /* returns number of sweeps */
  UINT max_samples = (UINT)samples; /* requests/returns number of samples */
  UINT samples_read = 0, n_samp_to_copy;
#ifdef _abf_DEBUG_
  int i;
#endif
  int phys_channel, cur_episode;
  long residual = 0;
  float *rawdata;
  void *data_ptr;

#ifdef _abf_DEBUG_
#ifdef _abf_USING_MEX_
  fprintf( logfile, "  using mex\n" );
#endif
  fprintf( logfile, "  opening file, chan %d samp %d, epi %d\n", channel, samples, episode );
#endif

  /* try to open file */
  if( !ABF_ReadOpen( filename, &fp, ABF_DATAFILE,
      &header, &max_samples, &max_episodes, &error_code ) )
    return _ABF_OPEN_ERROR;
    
  /* translate logical to physical channel number */
  phys_channel = header.nADCSamplingSeq[channel];
#ifdef _abf_DEBUG_
  fprintf( logfile, "  phys chan %d max_samp %d max_epi %d\n", phys_channel, max_samples, max_episodes );
#endif

  //mexPrintf("sample rate: %f\n", header.fADCSequenceInterval );
  *sampleInterval = header.fADCSequenceInterval;

  /* select starting episode number */
  if( header.nOperationMode == 5 || header.nOperationMode == 2) { /* episodic */
    if( episode == -1 ) /* user did not specify an episode, just choose first */
      cur_episode = 1;
    else cur_episode = episode;
  } else if( header.nOperationMode == 3 ) { /* gap-free */
    cur_episode = (int)(offset / max_samples + 1);
    residual = offset - (cur_episode-1) * max_samples;
  }
#ifdef _abf_DEBUG_
  fprintf( logfile, "  chose episode %d, residual %d\n", cur_episode, residual );
#endif

  /* allocate data for entire channel */
  rawdata = (float*)malloc( sizeof( float ) * max_samples );
#ifdef _abf_DEBUG_
  fprintf( logfile, "  rawdata size %d at %d\n", sizeof( float ) * max_samples, rawdata );
#endif
  /* may not be able to read all samples simultaneously, so loop until done */
  samples_read = -residual; /* correct for offset */
  while( (long)samples_read < samples && cur_episode <= max_episodes ) {
    /* read MAX_SAMPLES on PHYS_CHANNEL in EPISODE, store in RAWDATA */
      //mexPrintf("phys_channel: %d, cur_episode: %d, max_sample: %d\n", phys_channel, cur_episode, max_samples);
    if( !ABF_ReadChannel( fp, &header, phys_channel, cur_episode, rawdata, &max_samples, &error_code ) ) {
      if( error_code == ABF_EINVALIDCHANNEL ) return cleanup( rawdata, fp, _ABF_CHANNEL_ERROR );
      if( error_code == ABF_EEPISODERANGE ) return cleanup( rawdata, fp, _ABF_EPISODE_ERROR );
      return cleanup( rawdata, fp, _ABF_READ_ERROR );
    }
    
    /* take SAMPLES and copy them to DATA */
#ifdef _abf_DEBUG_
    fprintf( logfile, "  %d samples read to rawdata from episode %d:\n   ", max_samples, cur_episode );
    for( i = 0; i < 20; i++ )
      fprintf( logfile, "%d:%.2f ", i, *(rawdata+i) );
    fprintf( logfile, "\n  total samples read %d, trying episode %d\n", max_samples + samples_read, cur_episode+1 );
#endif
    /* find position in data array */
    data_ptr = (float*)data + channel*samples + samples_read + residual;
    /* copy rawdata into data */
    n_samp_to_copy = MIN( samples - samples_read,
                          MIN( max_samples, max_samples + samples_read )
                         );
#ifdef _abf_DEBUG_
    fprintf( logfile, "  data ptr at %d\n", data_ptr );
    fprintf( logfile, "  copying %d samples from %d\n", n_samp_to_copy, rawdata+residual );
#endif
    memcpy( data_ptr, rawdata + residual, sizeof( float ) * n_samp_to_copy );

    /* increment loop counters */
    samples_read += max_samples;
    cur_episode++;
    residual = 0; /* don't offset rawdata after first iteration */
  }

  return cleanup( rawdata, fp, _ABF_READ_SUCCESS );
}

int cleanup( void *memory, int fp, int return_code )
{
  int error_code;
  free( memory );
  ABF_Close( fp, &error_code );
  return return_code;
}
