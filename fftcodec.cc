#include <math.h>
//#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include "fftcodec.h"
#include "font.h"

FreeTypeFont g_font;

fftcodec_t* fftcodec_create( int p_nfftsize, int p_fsamplerate )
{
	fftcodec_t* l_pfftcodec;
	l_pfftcodec = (fftcodec_t*)malloc(sizeof(fftcodec_t));
	if(!l_pfftcodec){
		return NULL;
	}
	l_pfftcodec->m_nfftsize = p_nfftsize;
	l_pfftcodec->m_fsamplerate = p_fsamplerate;
	//
	// codec initialization
	//
	l_pfftcodec->m_codec_mode = CODEC_MODE_IDLE;
	l_pfftcodec->m_codec_flags = 0;
	pthread_mutex_init( &l_pfftcodec->m_codec_mutex, NULL );
	//
	// receiver initialization
	//
	l_pfftcodec->m_rx_mode = RX_MODE_IDLE;
	l_pfftcodec->m_pvsample = vsample_create( NWINDOW, NSUPERSAMPLE, p_fsamplerate );
	l_pfftcodec->m_irx = 0;
	l_pfftcodec->m_irxframe = 0;
	l_pfftcodec->m_rx_raw = (double*)malloc( sizeof(double) * p_nfftsize );
	l_pfftcodec->m_rx_x = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * p_nfftsize );
	l_pfftcodec->m_rx_X = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * p_nfftsize );
	l_pfftcodec->m_rx_Xframes = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * p_nfftsize * NFRAMES );
    l_pfftcodec->m_rx_impulse = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * p_nfftsize );
    l_pfftcodec->m_rx_Impulse = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * p_nfftsize );
    l_pfftcodec->m_rx_Sinv = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * (p_nfftsize/2+1) );
    l_pfftcodec->m_rx_Sinvd = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * (p_nfftsize + 1) );
    l_pfftcodec->m_rx_sinv = (double*)malloc( sizeof(double) * p_nfftsize * 2 );
	l_pfftcodec->m_rx_plan = fftw_plan_dft_1d( p_nfftsize,
							reinterpret_cast<fftw_complex*>(l_pfftcodec->m_rx_x),
							reinterpret_cast<fftw_complex*>(l_pfftcodec->m_rx_X),
                            FFTW_FORWARD, FFTW_MEASURE );
    l_pfftcodec->m_rx_impulse_plan = fftw_plan_dft_1d( p_nfftsize,
                            reinterpret_cast<fftw_complex*>(l_pfftcodec->m_rx_Impulse),
                            reinterpret_cast<fftw_complex*>(l_pfftcodec->m_rx_impulse),
                            FFTW_BACKWARD, FFTW_MEASURE );
    l_pfftcodec->m_rx_sinv_plan = fftw_plan_dft_c2r_1d(
                p_nfftsize,
                reinterpret_cast<fftw_complex*>(l_pfftcodec->m_rx_Sinv),
                l_pfftcodec->m_rx_sinv,
                FFTW_MEASURE);
    l_pfftcodec->m_rx_Sinvd_plan = fftw_plan_dft_r2c_1d(
                p_nfftsize*2,
                l_pfftcodec->m_rx_sinv,
                reinterpret_cast<fftw_complex*>(l_pfftcodec->m_rx_Sinvd),
                FFTW_MEASURE);
    spectrum_set( l_pfftcodec->m_rx_impulse, p_nfftsize, 0.0);
    bzero( l_pfftcodec->m_rx_sinv, sizeof(double) * p_nfftsize * 2 );
    l_pfftcodec->m_rx_nf_mean = (double*)malloc( sizeof(double) * p_nfftsize );
	l_pfftcodec->m_rx_nf_rms = (double*)malloc( sizeof(double) * p_nfftsize );
	l_pfftcodec->m_rx_ir_mean = (double*)malloc( sizeof(double) * p_nfftsize );
	l_pfftcodec->m_rx_ir_rms = (double*)malloc( sizeof(double) * p_nfftsize );
	l_pfftcodec->m_rx_Xarg = (double*)malloc( sizeof(double) * p_nfftsize );
	//
	// transmitter initialization
	//
	l_pfftcodec->m_tx_mode = TX_MODE_SILENCE;
	l_pfftcodec->m_itx = 0;
	l_pfftcodec->m_itxframe = 0;
	l_pfftcodec->m_itxsweep = 0;
	l_pfftcodec->m_tx_x = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * p_nfftsize );
	l_pfftcodec->m_tx_X = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * p_nfftsize );
	l_pfftcodec->m_tx_plan = fftw_plan_dft_1d( p_nfftsize,
							reinterpret_cast<fftw_complex*>(l_pfftcodec->m_tx_X),
							reinterpret_cast<fftw_complex*>(l_pfftcodec->m_tx_x),
                            FFTW_BACKWARD, FFTW_MEASURE );
	//
	// renderer initialization
	//
	l_pfftcodec->m_draw_mode = DRAW_MODE_TIME;
	l_pfftcodec->m_draw_src = DRAW_SRC_RX;
	l_pfftcodec->m_draw_flags = DRAW_FLAG_NF|DRAW_FLAG_IR;
	l_pfftcodec->m_iconstellation0 = 0;
	l_pfftcodec->m_iconstellation1 = p_nfftsize-1;
	l_pfftcodec->m_draw_constellation_max = 1.0;
	l_pfftcodec->m_draw_db_min = -180.0;
	l_pfftcodec->m_draw_db_max =  5.0;
	l_pfftcodec->m_draw_time_max = 1.0;
	l_pfftcodec->m_draw_data_x = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * p_nfftsize );
    l_pfftcodec->m_draw_data_X = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * p_nfftsize );
    l_pfftcodec->m_draw_data_impulse = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * p_nfftsize );

	return l_pfftcodec;
}

void fftcodec_store_sample( fftcodec_t* p_pfftcodec, double p_val )
{
	pthread_mutex_lock( &p_pfftcodec->m_codec_mutex );
	vsample_sample_in( p_pfftcodec->m_pvsample, p_val );
	while( p_pfftcodec->m_pvsample->m_ndata_out ){
		double l_vval = vsample_read_virtual( p_pfftcodec->m_pvsample );
		p_pfftcodec->m_rx_raw[p_pfftcodec->m_irx] = l_vval;
		if( ++p_pfftcodec->m_irx == p_pfftcodec->m_nfftsize ){
			// received a complete frame
			double* l_psrc = p_pfftcodec->m_rx_raw;
			std::complex<double>* l_pdst = p_pfftcodec->m_rx_x;
			long l_i;
			for( l_i=p_pfftcodec->m_nfftsize; l_i; l_i-- ){
				*(l_pdst++) = *(l_psrc++);
			}
			fftw_execute( p_pfftcodec->m_rx_plan );
			fftcodec_decode_spectrum( p_pfftcodec );
			p_pfftcodec->m_irx = 0;
		}
	}
	pthread_mutex_unlock( &p_pfftcodec->m_codec_mutex );
}

double fftcodec_read_sample( fftcodec_t* p_pfftcodec )
{
	pthread_mutex_lock( &p_pfftcodec->m_codec_mutex );
	if( p_pfftcodec->m_itx == 0 ){
		fftcodec_encode_spectrum( p_pfftcodec );
		fftw_execute( p_pfftcodec->m_tx_plan );
	}
	double l_val = std::real(p_pfftcodec->m_tx_x[p_pfftcodec->m_itx]);
	if(++p_pfftcodec->m_itx == p_pfftcodec->m_nfftsize){
		p_pfftcodec->m_itx = 0;
	}
	pthread_mutex_unlock( &p_pfftcodec->m_codec_mutex );
	return l_val;
}

void fftcodec_encode_spectrum( fftcodec_t* p_pfftcodec )
{
	switch( p_pfftcodec->m_tx_mode ){
		case TX_MODE_SILENCE:
			spectrum_set( p_pfftcodec->m_tx_X, p_pfftcodec->m_nfftsize, 0.0 );
			break;

		case TX_MODE_IMPULSE:
			spectrum_set( p_pfftcodec->m_tx_X, p_pfftcodec->m_nfftsize, 1.0/p_pfftcodec->m_nfftsize );
			break;
			
		case TX_MODE_COSINE:
        printf("encoding cosine spectrum at %ld\n",p_pfftcodec->m_itxsweep);
			spectrum_set( p_pfftcodec->m_tx_X, p_pfftcodec->m_nfftsize, 0.0 );
			if( p_pfftcodec->m_itxsweep == 0 ){
				p_pfftcodec->m_tx_X[0] = 1.0;
			}else if( p_pfftcodec->m_itxsweep == p_pfftcodec->m_nfftsize/2 ){
				p_pfftcodec->m_tx_X[ p_pfftcodec->m_nfftsize/2 ] = 1.0;
			}else{
				p_pfftcodec->m_tx_X[p_pfftcodec->m_itxsweep] = 0.5;
				p_pfftcodec->m_tx_X[p_pfftcodec->m_nfftsize-p_pfftcodec->m_itxsweep]=0.5;
			}
			break;
			
		case TX_MODE_SYNC_ONLY:
			spectrum_set( p_pfftcodec->m_tx_X, p_pfftcodec->m_nfftsize, 1.0/p_pfftcodec->m_nfftsize );
// 			spectrum_set( p_pfftcodec->m_tx_X, p_pfftcodec->m_nfftsize, 0.0 );
// 			p_pfftcodec->m_tx_X[1] = 1.0/p_pfftcodec->m_nfftsize;
// 			p_pfftcodec->m_tx_X[p_pfftcodec->m_nfftsize-1] = 1.0/p_pfftcodec->m_nfftsize;
			break;

		
	}
	p_pfftcodec->m_itxframe++;
	fftcodec_update_codec_mode( p_pfftcodec );
}

void fftcodec_decode_spectrum( fftcodec_t* p_pfftcodec )
{
	long n;
	std::complex<double>* l_pdata=p_pfftcodec->m_rx_X;
	double* l_parg0 = p_pfftcodec->m_rx_Xarg;
    long i=0;
	for(n=p_pfftcodec->m_nfftsize;n;--n){
        *l_pdata /= p_pfftcodec->m_nfftsize;
        if(i<4){
            //printf("<%lf,%lf> ", std::real(*l_pdata), std::imag(*l_pdata));
        }
		*l_parg0 = std::arg( *l_pdata );
		l_pdata++;
		l_parg0++;
        i++;
	}
    //printf("\n");
	l_parg0 = p_pfftcodec->m_rx_Xarg;
	double* l_parg1 = l_parg0+1;
	double  l_argmean = 0.0;
	for(n=p_pfftcodec->m_nfftsize-1;n;--n){
		double l_darg = *l_parg1 - *l_parg0;
		if( l_darg > M_PI ) {
			l_darg = l_darg-2*M_PI;
		}
		if( l_darg <-M_PI ) {
			l_darg = l_darg+2*M_PI;
		}
		l_darg = fabs(l_darg);
		l_argmean += l_darg;
		if( p_pfftcodec->m_codec_flags&CODEC_FLAG_PRINT ){
			printf("arg0:%f arg1:%f darg:%f\n", *l_parg0, *l_parg1, l_darg );
		}
		l_parg0++;
		l_parg1++;
	}
	if( p_pfftcodec->m_codec_flags&CODEC_FLAG_PRINT ){
		p_pfftcodec->m_codec_flags ^= CODEC_FLAG_PRINT;
	}
	
	l_argmean /= p_pfftcodec->m_nfftsize;
	p_pfftcodec->m_rx_Xargmean = l_argmean;
	switch( p_pfftcodec->m_rx_mode ){
    case RX_MODE_IDLE:
        break;
    case RX_MODE_SWEEP_TEST:
        if(p_pfftcodec->m_irxframe == NRXSWEEPFRAMES){
            int l_itxsweep = p_pfftcodec->m_itxsweep;
            std::complex<double> Ii = p_pfftcodec->m_rx_X[l_itxsweep];
            printf("sampling %i :<%lg,%lg> abs:%lg\n", l_itxsweep, std::real(Ii), std::imag(Ii), std::abs(Ii) );
            p_pfftcodec->m_rx_Impulse[l_itxsweep] = Ii;
            p_pfftcodec->m_rx_Impulse[p_pfftcodec->m_nfftsize-l_itxsweep] = std::conj(Ii);
        }
        if(++p_pfftcodec->m_irxframe == NSWEEPFRAMES){
            p_pfftcodec->m_irxframe=0;
        }
        break;

    case RX_MODE_NOISE_FLOOR:
        spectrum_copy(
            p_pfftcodec->m_rx_X,
            &p_pfftcodec->m_rx_Xframes[p_pfftcodec->m_irxframe*p_pfftcodec->m_nfftsize],
            p_pfftcodec->m_nfftsize );
        if(++p_pfftcodec->m_irxframe == NFRAMES){
            spectrum_stats(
                p_pfftcodec->m_rx_Xframes,
                p_pfftcodec->m_rx_nf_mean,
                p_pfftcodec->m_rx_nf_rms,
                p_pfftcodec->m_nfftsize,
                NFRAMES );
            p_pfftcodec->m_rx_mode = RX_MODE_IDLE;
        }
        break;

    case RX_MODE_IMPULSE_RESPONSE:
        spectrum_copy(
            p_pfftcodec->m_rx_X,
            &p_pfftcodec->m_rx_Xframes[p_pfftcodec->m_irxframe*p_pfftcodec->m_nfftsize],
            p_pfftcodec->m_nfftsize );
        if(++p_pfftcodec->m_irxframe == NFRAMES){
            spectrum_stats(
                p_pfftcodec->m_rx_Xframes,
                p_pfftcodec->m_rx_ir_mean,
                p_pfftcodec->m_rx_ir_rms,
                p_pfftcodec->m_nfftsize,
                NFRAMES );
            p_pfftcodec->m_rx_mode = RX_MODE_IDLE;
        }
        break;

    case RX_MODE_SYNC_ONLY:
        fftcodec_sync( p_pfftcodec );
        break;
	}
	fftcodec_update_codec_mode( p_pfftcodec );
}

void fftcodec_sync( fftcodec_t* p_pfftcodec )
{
	// estimate the sample error from the pilot tone phase
	double l_error = p_pfftcodec->m_rx_Xarg[1]/2/M_PI*p_pfftcodec->m_nfftsize;
	double l_delta_fv,l_delta_error;
	
	switch( p_pfftcodec->m_sync_mode ){
		case SYNC_MODE_OFFSET_ANALYZE:
			// use the error to adjust the virtual sampler frequency
			l_delta_fv = l_error*p_pfftcodec->m_fsamplerate/p_pfftcodec->m_nfftsize;
			p_pfftcodec->m_pvsample->m_fv += l_delta_fv;
			p_pfftcodec->m_delta_fv = l_delta_fv;
			p_pfftcodec->m_sync_mode = SYNC_MODE_OFFSET_RESET;
			break;
			
		case SYNC_MODE_OFFSET_RESET:
			// revert the virtual sampler frequency back
			p_pfftcodec->m_pvsample->m_fv -= p_pfftcodec->m_delta_fv;
			p_pfftcodec->m_sync_mode = SYNC_MODE_FREQUENCY_ANALYZE;
			break;
			
		case SYNC_MODE_FREQUENCY_ANALYZE:
			// store the error
			p_pfftcodec->m_error = l_error;
			p_pfftcodec->m_sync_mode = SYNC_MODE_FREQUENCY_APPLY;
			break;
			
		case SYNC_MODE_FREQUENCY_APPLY:
			// calculate the error differential
			l_delta_error = l_error - p_pfftcodec->m_error;
			p_pfftcodec->m_pvsample->m_fv += l_delta_error*p_pfftcodec->m_fsamplerate/p_pfftcodec->m_nfftsize;
			p_pfftcodec->m_sync_mode = SYNC_MODE_OFFSET_ANALYZE;
			break;
	}
}

void fftcodec_codec_mode( fftcodec_t* p_pfftcodec, codec_mode_t p_mode )
{
	pthread_mutex_lock( &p_pfftcodec->m_codec_mutex );
	if( p_pfftcodec->m_codec_flags&CODEC_FLAG_MODE_QUEUED ){
		pthread_mutex_unlock( &p_pfftcodec->m_codec_mutex );
		return;
	}
	p_pfftcodec->m_codec_flags |= CODEC_FLAG_MODE_QUEUED;
	p_pfftcodec->m_codec_mode_next = p_mode;
	pthread_mutex_unlock( &p_pfftcodec->m_codec_mutex );
}


void fftcodec_init_codec_mode( fftcodec_t* p_pfftcodec )
{
	if( (p_pfftcodec->m_codec_flags&CODEC_FLAG_MODE_QUEUED) == 0 )
		return;
	p_pfftcodec->m_codec_flags ^= CODEC_FLAG_MODE_QUEUED;
	switch( p_pfftcodec->m_codec_mode_next ){
		case CODEC_MODE_IDLE:
			fftcodec_init_tx_mode( p_pfftcodec, TX_MODE_SILENCE );
			fftcodec_init_rx_mode( p_pfftcodec, RX_MODE_IDLE );
			p_pfftcodec->m_codec_mode = CODEC_MODE_IDLE;
			break;
		case CODEC_MODE_IMPULSE_TEST:
			fftcodec_init_tx_mode( p_pfftcodec, TX_MODE_IMPULSE );
			fftcodec_init_rx_mode( p_pfftcodec, RX_MODE_IDLE );
			p_pfftcodec->m_codec_mode = CODEC_MODE_IMPULSE_TEST;
			break;
		case CODEC_MODE_SWEEP_TEST:
			fftcodec_init_tx_mode( p_pfftcodec, TX_MODE_COSINE );
            fftcodec_init_rx_mode( p_pfftcodec, RX_MODE_SWEEP_TEST );
			p_pfftcodec->m_codec_mode = CODEC_MODE_SWEEP_TEST;
			p_pfftcodec->m_itxsweep = 0;
			break;
		case CODEC_MODE_NOISE_FLOOR:
		case CODEC_MODE_NOISE_FLOOR_INIT:
			fftcodec_init_tx_mode( p_pfftcodec, TX_MODE_SILENCE );
			fftcodec_init_rx_mode( p_pfftcodec, RX_MODE_IDLE );
			p_pfftcodec->m_codec_mode = CODEC_MODE_NOISE_FLOOR_INIT;
			p_pfftcodec->m_codec_mode_next = CODEC_MODE_IDLE;
			p_pfftcodec->m_codec_flags |= CODEC_FLAG_MODE_QUEUED;
			break;
		case CODEC_MODE_IMPULSE_RESPONSE:
		case CODEC_MODE_IMPULSE_RESPONSE_INIT:
			fftcodec_init_tx_mode( p_pfftcodec, TX_MODE_IMPULSE );
			fftcodec_init_rx_mode( p_pfftcodec, RX_MODE_IDLE );
			p_pfftcodec->m_codec_mode = CODEC_MODE_IMPULSE_RESPONSE_INIT;
			p_pfftcodec->m_codec_mode_next = CODEC_MODE_IDLE;
			p_pfftcodec->m_codec_flags |= CODEC_FLAG_MODE_QUEUED;
			break;
		case CODEC_MODE_SYNC_ONLY:
		case CODEC_MODE_SYNC_ONLY_INIT:
			fftcodec_init_tx_mode( p_pfftcodec, TX_MODE_SYNC_ONLY );
			fftcodec_init_rx_mode( p_pfftcodec, RX_MODE_IDLE );
			p_pfftcodec->m_codec_mode = CODEC_MODE_SYNC_ONLY_INIT;
			p_pfftcodec->m_pvsample->m_fv = p_pfftcodec->m_pvsample->m_fs;
			break;
	}
}

void fftcodec_init_rx_mode( fftcodec_t* p_pfftcodec, rx_mode_t p_rxmode )
{
	p_pfftcodec->m_rx_mode = p_rxmode;
	p_pfftcodec->m_irxframe = 0;
	if( p_rxmode == RX_MODE_SYNC_ONLY )
		p_pfftcodec->m_sync_mode = SYNC_MODE_OFFSET_ANALYZE;
    if( p_rxmode == RX_MODE_SWEEP_TEST ){
        spectrum_set( p_pfftcodec->m_rx_Impulse, p_pfftcodec->m_nfftsize, 0.0 );
    }
}

void fftcodec_init_tx_mode( fftcodec_t* p_pfftcodec, tx_mode_t p_txmode )
{
	p_pfftcodec->m_tx_mode = p_txmode;
	p_pfftcodec->m_itxframe = 0;
}

void fftcodec_update_codec_mode( fftcodec_t* p_pfftcodec )
{
	switch( p_pfftcodec->m_codec_mode ){
		case CODEC_MODE_IDLE:
			fftcodec_init_codec_mode( p_pfftcodec );
			break;
		case CODEC_MODE_IMPULSE_TEST:
			fftcodec_init_codec_mode( p_pfftcodec );
			break;
		case CODEC_MODE_SWEEP_TEST:
			if( p_pfftcodec->m_itxframe == NSWEEPFRAMES ){
				if( p_pfftcodec->m_itxsweep == p_pfftcodec->m_nfftsize/2 ){
                    printf("Sweep test complete. Saving results.\n");
                    fftw_execute( p_pfftcodec->m_rx_impulse_plan );
                    fftcodec_init_tx_mode( p_pfftcodec, TX_MODE_SILENCE );
                    fftcodec_init_rx_mode( p_pfftcodec, RX_MODE_IDLE );
                    p_pfftcodec->m_codec_mode = CODEC_MODE_IDLE;
					p_pfftcodec->m_itxsweep = 0;
                    FILE* rfile = fopen("response.txt","w");
                    if(rfile){
                        fprintf(rfile, "FFT size:%d Sample rate:%lf\n",
                                p_pfftcodec->m_nfftsize,
                                p_pfftcodec->m_fsamplerate);
                        std::complex<double> *l_pI = p_pfftcodec->m_rx_Impulse;
                        for(int i=0;i<=p_pfftcodec->m_nfftsize/2;i++){
                            fprintf(rfile, "%20.19lg, %20.19lg\n",
                                    std::real(*l_pI),std::imag(*l_pI));
                            l_pI++;
                        }
                        fclose(rfile);
                    }
                    // compute the geometric mean of the impulse spectrum
                    double gm=0;
                    for(int i=1;i<=p_pfftcodec->m_nfftsize/2;i++){
                        gm += log(abs(p_pfftcodec->m_rx_Impulse[i]));
                    }
                    gm /= p_pfftcodec->m_nfftsize/2;
                    gm = exp(gm);
                    // scale the impulse spectrum by the geometric mean
                    for(int i=0;i<=p_pfftcodec->m_nfftsize/2;i++){
                        p_pfftcodec->m_rx_Sinv[i] = std::complex<double>( abs(p_pfftcodec->m_rx_Impulse[i])/gm, 0.0 );
                    }
                    // compute the inverse and cull certain frequencies;
                    for(int i=0;i<=p_pfftcodec->m_nfftsize/2;i++){
                        if(i==0 || i==p_pfftcodec->m_nfftsize/2){
                            p_pfftcodec->m_rx_Sinv[i] = std::complex<double>(1.0, 0.0);
                        }else{
                            double a = abs(p_pfftcodec->m_rx_Sinv[i]);
                            if(a<0.1){
                                p_pfftcodec->m_rx_Sinv[i] = std::complex<double>(1.0, 0.0);
                            }else{
                                p_pfftcodec->m_rx_Sinv[i] = std::complex<double>(1.0/a, 0.0);
                            }
                        }
                    }
                    // write out the single inverse spectrum
                    FILE* ifile = fopen("Sinv.txt","w");
                    if(ifile){
                        fprintf(ifile, "FFTsize:%d\n",
                                p_pfftcodec->m_nfftsize);
                        std::complex<double> *l_pI = p_pfftcodec->m_rx_Sinv;
                        for(int i=0;i<=p_pfftcodec->m_nfftsize/2;i++){
                            fprintf(ifile, "%20.19lg, %20.19lg\n",
                                    std::real(*l_pI),std::imag(*l_pI));
                            l_pI++;
                        }
                        fclose(ifile);
                    }

                    // shift in the time domain half way into the buffer
                    for(int i=0;i<=p_pfftcodec->m_nfftsize/2;i++){
                        std::complex<double> k = exp( std::complex<double>(0.0, (double)i*2.0*M_PI*0.5) );
                        p_pfftcodec->m_rx_Sinv[i] *= k;
                    }
                    // convert to the time domain
                    fftw_execute( p_pfftcodec->m_rx_sinv_plan );
                    // convert back to the double size frequency domain
                    fftw_execute( p_pfftcodec->m_rx_Sinvd_plan );
                    // write out the double inverse spectrum
                    ifile = fopen("Sinvd.txt","w");
                    if(ifile){
                        fprintf(ifile, "FFTsize:%d\n",
                                p_pfftcodec->m_nfftsize);
                        std::complex<double> *l_pI = p_pfftcodec->m_rx_Sinvd;
                        for(int i=0;i<=p_pfftcodec->m_nfftsize;i++){
                            *l_pI /= p_pfftcodec->m_nfftsize * 2;
                            fprintf(ifile, "%20.19lg, %20.19lg\n",
                                    std::real(*l_pI),std::imag(*l_pI));
                            l_pI++;
                        }
                        fclose(ifile);
                    }
				}else{
                    fftcodec_init_tx_mode( p_pfftcodec, TX_MODE_COSINE );
                    p_pfftcodec->m_itxsweep++;
				}
			}
			fftcodec_init_codec_mode( p_pfftcodec );
			break;
		case CODEC_MODE_NOISE_FLOOR_INIT:
			if( p_pfftcodec->m_itxframe == NINITFRAMES ){
				fftcodec_init_rx_mode( p_pfftcodec, RX_MODE_NOISE_FLOOR );
				p_pfftcodec->m_codec_mode = CODEC_MODE_NOISE_FLOOR;
			}
			break;
		case CODEC_MODE_NOISE_FLOOR:
			if( p_pfftcodec->m_rx_mode == RX_MODE_IDLE ){
				p_pfftcodec->m_codec_flags |= CODEC_FLAG_NF_COMPLETE;
				fftcodec_init_codec_mode( p_pfftcodec );
			}
			break;
		case CODEC_MODE_IMPULSE_RESPONSE_INIT:
			if( p_pfftcodec->m_itxframe == NINITFRAMES ){
				fftcodec_init_rx_mode( p_pfftcodec, RX_MODE_IMPULSE_RESPONSE );
				p_pfftcodec->m_codec_mode = CODEC_MODE_IMPULSE_RESPONSE;
			}
			break;
		case CODEC_MODE_IMPULSE_RESPONSE:
			if( p_pfftcodec->m_rx_mode == RX_MODE_IDLE ){
				p_pfftcodec->m_codec_flags |= CODEC_FLAG_IR_COMPLETE;
				fftcodec_init_codec_mode( p_pfftcodec );
			}
			break;
		case CODEC_MODE_SYNC_ONLY_INIT:
			if( p_pfftcodec->m_itxframe == NINITFRAMES ){
				fftcodec_init_rx_mode( p_pfftcodec, RX_MODE_SYNC_ONLY );
				p_pfftcodec->m_codec_mode = CODEC_MODE_SYNC_ONLY;
			}
			break;
		case CODEC_MODE_SYNC_ONLY:
			fftcodec_init_codec_mode( p_pfftcodec );
			break;
	}
}

void fftcodec_draw_init( fftcodec_t* p_pfftcodec, Display* p_pDisplay )
{
// 	p_pfftcodec->m_pfont = XLoadQueryFont( p_pDisplay, "-misc-fixed-*-*-*-*-*-140-*-*-*-*-*-*" );
// 	if( !p_pfftcodec->m_pfont ){
// 		printf("font not found\n");
// 		return;
// 	}
// 	Font l_fid;
// 	unsigned int l_first, l_last;
// 	l_fid = p_pfftcodec->m_pfont->fid;
// 	l_first = p_pfftcodec->m_pfont->min_char_or_byte2;
// 	l_last  = p_pfftcodec->m_pfont->max_char_or_byte2;
// 	g_font_base = glGenLists( (GLuint)l_last + 1 );
// 	if( !g_font_base ){
// 		printf("unable to allocate display lists for font\n");
// 		return;
// 	}
// 	GLenum l_error;
// 	do{
// 	}while( (l_error=glGetError())!=GL_NO_ERROR );
// 	glXUseXFont( l_fid, l_first, l_last - l_first + 1, g_font_base + l_first );
// 	l_error = glGetError();
// 	if( l_error!=GL_NO_ERROR ){
// 		printf( "glXUseXFont error: %d\n", l_error );
// 	}
// 	p_pfftcodec->m_text_height = p_pfftcodec->m_pfont->ascent + p_pfftcodec->m_pfont->descent;

// 	g_font.LoadOutline(  "/usr/share/fonts/truetype/ubuntu-font-family/UbuntuMono-R.ttf",
// 										12,
// 									  Pixel32( 0,0,0 ),
// 										Pixel32( 255, 255, 255 ),
// 										1.5 );
	g_font.LoadOutline("UbuntuMono-R.ttf", 12, Pixel32(0,0,0), Pixel32(128,128,128), 1.5);
	p_pfftcodec->m_text_height = 12;
	
	
	
	GLenum l_error;
	do{
	}while( (l_error=glGetError())!=GL_NO_ERROR );
	double l_range[2];
	glGetDoublev( GL_SMOOTH_POINT_SIZE_RANGE, l_range );
	l_error = glGetError();
	if( l_error == GL_NO_ERROR ){
		printf("Smooth Point size range:%f  %f\n",l_range[0],l_range[1]);
	}else{
		printf("glGetDoublev error:%d\n",l_error);
	}

}

const char* fftcodec_draw_mode_string( fftcodec_t* p_pfftcodec )
{
	switch( p_pfftcodec->m_draw_mode ){
    case DRAW_MODE_TIME:
        return "time";
    case DRAW_MODE_IMPULSE:
        return "impulse";
    case DRAW_MODE_FREQUENCY_ABS_LINEAR:
        return "frequency(linear)";
    case DRAW_MODE_FREQUENCY_ABS_LOG:
        return "frequency(log)";
    case DRAW_MODE_FREQUENCY_CONSTELLATION:
        return "constellation";
	}
}

const char* fftcodec_codec_mode_string( fftcodec_t* p_pfftcodec )
{
	switch( p_pfftcodec->m_codec_mode ){
		case CODEC_MODE_IDLE:
			return "Idle";
		case CODEC_MODE_IMPULSE_TEST:
			return "Impulse Test";
		case CODEC_MODE_SWEEP_TEST:
			return "Sweep Test";
		case CODEC_MODE_NOISE_FLOOR_INIT:
			return "Noise Floor Start";
		case CODEC_MODE_NOISE_FLOOR:
			return "Noise Floor";
		case CODEC_MODE_IMPULSE_RESPONSE_INIT:
			return "Impulse Response Start";
		case CODEC_MODE_IMPULSE_RESPONSE:
			return "Impulse Response";
		case CODEC_MODE_SYNC_ONLY_INIT:
			return "Sync Only Start";
		case CODEC_MODE_SYNC_ONLY:
			return "Sync Only";
	}
}

const char* fftcodec_rx_mode_string( fftcodec_t* p_pfftcodec )
{
	switch( p_pfftcodec->m_rx_mode ){
    case RX_MODE_IDLE:
        return "Idle";
    case RX_MODE_SWEEP_TEST:
        return "Sweep Test";
    case RX_MODE_NOISE_FLOOR:
        return "Noise Floor";
    case RX_MODE_IMPULSE_RESPONSE:
        return "Impulse Response";
    case RX_MODE_SYNC_ONLY:
        return "Sync Only";
	}
}

const char* fftcodec_tx_mode_string( fftcodec_t* p_pfftcodec )
{
	switch( p_pfftcodec->m_tx_mode ){
		case TX_MODE_SILENCE:
			return "Silence";
		case TX_MODE_IMPULSE:
			return "Impulse";
		case TX_MODE_COSINE:
			return "Cosine";
		case TX_MODE_SYNC_ONLY:
			return "Sync Only";
	}
}

void fftcodec_draw( fftcodec_t* p_pfftcodec, int p_width, int p_height )
{
	char l_str[1024];
	std::complex<double>* l_psrcx;
	std::complex<double>* l_psrcX;
	switch( p_pfftcodec->m_draw_src ){
		case DRAW_SRC_RX:
			l_psrcx = p_pfftcodec->m_rx_x;
			l_psrcX = p_pfftcodec->m_rx_X;
			break;
		case DRAW_SRC_TX:
			l_psrcx = p_pfftcodec->m_tx_x;
			l_psrcX = p_pfftcodec->m_tx_X;
			break;
	}
	pthread_mutex_lock( &p_pfftcodec->m_codec_mutex );
	spectrum_copy( l_psrcx, p_pfftcodec->m_draw_data_x, p_pfftcodec->m_nfftsize );
	spectrum_copy( l_psrcX, p_pfftcodec->m_draw_data_X, p_pfftcodec->m_nfftsize );
    spectrum_copy( p_pfftcodec->m_rx_impulse, p_pfftcodec->m_draw_data_impulse, p_pfftcodec->m_nfftsize );
	pthread_mutex_unlock( &p_pfftcodec->m_codec_mutex );
	l_psrcx = p_pfftcodec->m_draw_data_x;
	l_psrcX = p_pfftcodec->m_draw_data_X;
	std::complex<double>* l_pval;
	long n,i,N;
	GLdouble l_color[3];
	long l_npoints = p_pfftcodec->m_nfftsize/2 + 1;
	glDrawBuffer(GL_BACK);
	glClearColor( 0.0, 0.0, 0.0, 0.0 );
	glClear( GL_COLOR_BUFFER_BIT );
	glEnable( GL_LINE_SMOOTH );
	glEnable( GL_BLEND );
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	glViewport( 0, 0, (GLsizei)p_width, (GLsizei)p_height );
	glColor3d( 1.0, 1.0, 1.0 );
	glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
	glHint( GL_POINT_SMOOTH_HINT, GL_NICEST );
	switch( p_pfftcodec->m_draw_mode ){
    case DRAW_MODE_TIME:
        glMatrixMode( GL_PROJECTION );
        glLoadIdentity();
        gluOrtho2D( 0.0, 1.0, -p_pfftcodec->m_draw_time_max, p_pfftcodec->m_draw_time_max );
        glMatrixMode( GL_MODELVIEW );
        glLoadIdentity();
        l_pval = l_psrcx;
        glColor4d( 0.5, 1.0, 0.5, 1.0 );
        glBegin( GL_LINE_STRIP );
        for(n=0;n<p_pfftcodec->m_nfftsize;n++){
            double l_x = (double)n/(p_pfftcodec->m_nfftsize-1);
            double l_y = std::real(*l_pval);
            glVertex2d( l_x, l_y );
            l_pval++;
        }
        glEnd();
        glColor3d(1.0,1.0,1.0);
        g_font.Printf( 0, p_height - p_pfftcodec->m_text_height, "Scale:%7.5f", p_pfftcodec->m_draw_time_max );
        break;

    case DRAW_MODE_IMPULSE:
        glMatrixMode( GL_PROJECTION );
        glLoadIdentity();
        gluOrtho2D( 0.0, 1.0, -p_pfftcodec->m_draw_time_max, p_pfftcodec->m_draw_time_max );
        glMatrixMode( GL_MODELVIEW );
        glLoadIdentity();
        l_pval = p_pfftcodec->m_draw_data_impulse;
        glColor4d( 0.5, 1.0, 0.5, 0.5 );
        glBegin( GL_LINE_STRIP );
        for(n=0;n<p_pfftcodec->m_nfftsize;n++){
            double l_x = (double)n/(p_pfftcodec->m_nfftsize-1);
            double l_y = std::real(*l_pval);
            glVertex2d( l_x, l_y );
            l_pval++;
        }
        glEnd();
        glColor3d(1.0,1.0,1.0);
// 			sprintf( l_str, "Scale:%7.5f", p_pfftcodec->m_draw_time_max );
// 			glWindowPos2i( (GLint)0, (GLint)0 );
// 			glListBase( g_font_base );
// 			glCallLists( strlen(l_str), GL_UNSIGNED_BYTE, (GLubyte*) l_str );
        g_font.Printf( 0, p_height - p_pfftcodec->m_text_height, "Scale:%7.5f", p_pfftcodec->m_draw_time_max );
        break;

    case DRAW_MODE_FREQUENCY_ABS_LINEAR:
        glMatrixMode( GL_PROJECTION );
        glLoadIdentity();
        gluOrtho2D( 0.0, 1.0, 0.0, 1.0 );
        glMatrixMode( GL_MODELVIEW );
        glLoadIdentity();
        // draw the spectrum
        spectrum_draw_abs( l_psrcX, l_npoints );
        // draw the noise floor
        if( p_pfftcodec->m_codec_flags&CODEC_FLAG_NF_COMPLETE ){
            if( p_pfftcodec->m_draw_flags&DRAW_FLAG_NF ){
                // draw the noise floor
                glColor3d( 1.0, 0.5, 1.0 );
                dspectrum_draw_abs( p_pfftcodec->m_rx_nf_mean, l_npoints );
                dspectrum_draw_abs( p_pfftcodec->m_rx_nf_rms, l_npoints );
            }
        }
        // draw the impulse response
        if( p_pfftcodec->m_codec_flags&CODEC_FLAG_IR_COMPLETE ){
            if( p_pfftcodec->m_draw_flags&DRAW_FLAG_IR ){
                // draw the impulse response
                glColor3d( 0.5, 1.0, 1.0 );
                dspectrum_draw_abs( p_pfftcodec->m_rx_ir_mean, l_npoints );
                dspectrum_draw_abs( p_pfftcodec->m_rx_ir_rms, l_npoints );
            }
        }
        break;

    case DRAW_MODE_FREQUENCY_ABS_LOG:
        glMatrixMode( GL_PROJECTION );
        glLoadIdentity();
        gluOrtho2D( 0.0, 1.0, p_pfftcodec->m_draw_db_min, p_pfftcodec->m_draw_db_max );
        glMatrixMode( GL_MODELVIEW );
        glLoadIdentity();
        spectrum_draw_log( l_psrcX, l_npoints );
        // draw the noise floor
        if( p_pfftcodec->m_codec_flags&CODEC_FLAG_NF_COMPLETE ){
            if( p_pfftcodec->m_draw_flags&DRAW_FLAG_NF ){
                // draw the noise floor
                glColor3d( 1.0, 0.5, 1.0 );
                dspectrum_draw_log_bounds( p_pfftcodec->m_rx_nf_mean, p_pfftcodec->m_rx_nf_rms, l_npoints );
            }
        }
        // draw the impulse response
        if( p_pfftcodec->m_codec_flags&CODEC_FLAG_IR_COMPLETE ){
            if( p_pfftcodec->m_draw_flags&DRAW_FLAG_IR ){
                // draw the impulse response
                glColor3d( 0.5, 1.0, 1.0 );
                dspectrum_draw_log_bounds( p_pfftcodec->m_rx_ir_mean, p_pfftcodec->m_rx_ir_rms, l_npoints );
            }
        }
        glColor3d(1.0,1.0,1.0);
// 			sprintf( l_str, "dB max:%5.1f dB min:%5.1f", p_pfftcodec->m_draw_db_max, p_pfftcodec->m_draw_db_min );
// 			glWindowPos2i( (GLint)0, (GLint)0 );
// 			glListBase( g_font_base );
// 			glCallLists( strlen(l_str), GL_UNSIGNED_BYTE, (GLubyte*) l_str );
        g_font.Printf( 0, p_height - p_pfftcodec->m_text_height,
                                    "dB max:%5.1f dB min:%5.1f", p_pfftcodec->m_draw_db_max,
                                    p_pfftcodec->m_draw_db_min );
        break;

    case DRAW_MODE_FREQUENCY_CONSTELLATION:
        glEnable( GL_POINT_SMOOTH );
        glPointSize( 3.0 );
        glMatrixMode( GL_PROJECTION );
        glLoadIdentity();
        gluOrtho2D(
            -p_pfftcodec->m_draw_constellation_max*p_width/p_height,
            p_pfftcodec->m_draw_constellation_max*p_width/p_height,
            -p_pfftcodec->m_draw_constellation_max,
            p_pfftcodec->m_draw_constellation_max );
        glMatrixMode( GL_MODELVIEW );
        glLoadIdentity();
        l_pval = &l_psrcX[p_pfftcodec->m_iconstellation0];
        N = p_pfftcodec->m_iconstellation1 - p_pfftcodec->m_iconstellation0 + 1;
        glBegin( GL_POINTS );
        for(n=p_pfftcodec->m_iconstellation0,i=0;n<=p_pfftcodec->m_iconstellation1;n++,i++){
            double l_hue = (double)i/N;
            color_set_hsv( l_color, l_hue, 0.75, 1.0 );
            glColor3dv( l_color );
            glVertex2d( std::real(*l_pval), std::imag(*l_pval) );
            l_pval++;
        }
        glEnd();
        glColor3d(1.0,1.0,1.0);
// 			sprintf( l_str, "Scale:%8.6f MRA:%+8.6f", p_pfftcodec->m_draw_constellation_max, p_pfftcodec->m_rx_Xargmean );
// 			glWindowPos2i( (GLint)0, (GLint)0 );
// 			glListBase( g_font_base );
// 			glCallLists( strlen(l_str), GL_UNSIGNED_BYTE, (GLubyte*) l_str );
        g_font.Printf( 0, p_height - p_pfftcodec->m_text_height, "Scale:%8.6f MRA:%+8.6f",
                                    p_pfftcodec->m_draw_constellation_max,
                           p_pfftcodec->m_rx_Xargmean );
        break;
	}
	// draw the title information
	glColor3d(1.0,1.0,1.0);
// 	sprintf( l_str, "Modes Draw:%s CODEC:%s RX:%s TX:%s",
// 					 fftcodec_draw_mode_string( p_pfftcodec ),
// 					 fftcodec_codec_mode_string( p_pfftcodec ),
// 					 fftcodec_rx_mode_string( p_pfftcodec ),
// 					 fftcodec_tx_mode_string( p_pfftcodec ) );
// 	glWindowPos2i( (GLint)0, (GLint)(p_height - p_pfftcodec->m_text_height) );
// 	glListBase( g_font_base );
// 	glCallLists( strlen(l_str), GL_UNSIGNED_BYTE, (GLubyte*) l_str );
	g_font.Printf( 0, p_pfftcodec->m_text_height,
								"Modes Draw:%s CODEC:%s RX:%s TX:%s",
							 fftcodec_draw_mode_string( p_pfftcodec ),
								fftcodec_codec_mode_string( p_pfftcodec ),
								fftcodec_rx_mode_string( p_pfftcodec ),
								fftcodec_tx_mode_string( p_pfftcodec ) );
}

double color_func( double p_x )
{
	while( p_x <-0.5 ) p_x += 1.0;
	while( p_x > 0.5 ) p_x -= 1.0;
	if( p_x < 0.0 ) p_x = -p_x;
	p_x*=6;
	if( p_x <= 1.0 )
		return 1.0;
	else if( p_x <= 2.0 )
		return 2.0 - p_x;
	else
		return 0.0;
}

void color_set_hsv( GLdouble* p_c, double p_h, double p_s, double p_v )
{
	double l_r = color_func( p_h );
	double l_g = color_func( p_h - 1.0/3.0 );
	double l_b = color_func( p_h - 2.0/3.0 );
	p_s = 1.0 - p_s;
	l_r += (1.0 - l_r)*p_s;
	l_g += (1.0 - l_g)*p_s;
	l_b += (1.0 - l_b)*p_s;
	l_r*=p_v;
	l_g*=p_v;
	l_b*=p_v;
	p_c[0] = l_r;
	p_c[1] = l_g;
	p_c[2] = l_b;
}

void color_set( GLdouble* p_v, double p_r, double p_g, double p_b )
{
	p_v[0] = p_r;
	p_v[1] = p_g;
	p_v[2] = p_b;
}

void spectrum_draw_abs( std::complex<double>* p_psrc, long p_npoints )
{
	std::complex<double>* l_pdata = p_psrc;
	long n;
	glBegin( GL_LINE_STRIP );
	for(n=0;n<p_npoints;n++){
		double l_x = (double)n/(p_npoints-1);
		double l_y = std::abs(*l_pdata);
		glVertex2d( l_x, l_y );
		l_pdata++;
	}
	glEnd();
}

void spectrum_draw_log( std::complex<double>* p_psrc, long p_npoints )
{
	std::complex<double>* l_pdata = p_psrc;
	long n;
	glBegin( GL_LINE_STRIP );
	for(n=0;n<p_npoints;n++){
		double l_x = (double)n/(p_npoints-1);
		double l_y = std::abs(*l_pdata);
        
		if( l_y < 1e-9 ) l_y=1e-9;
		l_y = 20.0*log10(l_y);
		glVertex2d( l_x, l_y );
		l_pdata++;
	}
	glEnd();
}

void dspectrum_draw_abs( double* p_psrc, long p_npoints )
{
	double* l_pdval = p_psrc;
	long n;
	glBegin( GL_LINE_STRIP );
	for(n=0;n<p_npoints;n++){
		double l_x = (double)n/(p_npoints-1);
		glVertex2d( l_x, *l_pdval);
		l_pdval++;
	}
	glEnd();
}

void dspectrum_draw_log( double* p_psrc, long p_npoints )
{
	double* l_pdval = p_psrc;
	long n;
	glBegin( GL_LINE_STRIP );
	for(n=0;n<p_npoints;n++){
		double l_x = (double)n/(p_npoints-1);
		double l_y = *l_pdval;
		if( l_y < 1e-9 ) l_y=1e-9;
		l_y = 20.0*log10(l_y);
		glVertex2d( l_x, l_y );
		l_pdval++;
	}
	glEnd();
}

void dspectrum_draw_log_bounds( double* p_psrc, double* p_prms, long p_npoints )
{
	double* l_pdval = p_psrc;
	double* l_prms = p_prms;
	long n;
	glBegin( GL_LINE_STRIP );
	for(n=0;n<p_npoints;n++){
		double l_x = (double)n/(p_npoints-1);
		double l_y = *l_pdval;
		if( l_y < 1e-9 ) l_y=1e-9;
		l_y = 20.0*log10(l_y);
		glVertex2d( l_x, l_y );
		l_pdval++;
	}
	glEnd();
	l_pdval = p_psrc;
	glBegin( GL_LINE_STRIP );
	for(n=0;n<p_npoints;n++){
		double l_x = (double)n/(p_npoints-1);
		double l_y = *l_pdval + *l_prms;
		if( l_y < 1e-9 ) l_y=1e-9;
		l_y = 20.0*log10(l_y);
		glVertex2d( l_x, l_y );
		l_pdval++;
		l_prms++;
	}
	glEnd();
	l_pdval = p_psrc;
	l_prms = p_prms;
	glBegin( GL_LINE_STRIP );
	for(n=0;n<p_npoints;n++){
		double l_x = (double)n/(p_npoints-1);
		double l_y = *l_pdval - *l_prms;
		if( l_y < 1e-9 ) l_y=1e-9;
		l_y = 20.0*log10(l_y);
		glVertex2d( l_x, l_y );
		l_pdval++;
		l_prms++;
	}
	glEnd();
}



void spectrum_stats( std::complex<double>* p_psrc, double* p_pmean, double* p_prms, long p_nfftsize, long p_nframes )
{
	long n;
	// allocate memory for the absolute values
	double* l_pcabs0 = (double*)malloc(sizeof(double)*p_nfftsize*p_nframes);
	if(!l_pcabs0)return;
	// compute the absolute values
	double* l_pcabs = l_pcabs0;
	std::complex<double>* l_pdata = p_psrc;
	for(n=p_nfftsize*p_nframes;n;--n){
		*l_pcabs = std::abs(*l_pdata);
		l_pcabs++;
		l_pdata++;
	}
	// calculate the mean
	double* l_pmean = p_pmean;
	double* l_pcabs_start = l_pcabs0;
	for(n=p_nfftsize;n;--n){
		*l_pmean = 0.0;
		long f;
		l_pcabs = l_pcabs_start;
		for(f=p_nframes;f;--f){
			*l_pmean += *l_pcabs;
			l_pcabs += p_nfftsize;
		}
		*l_pmean /= p_nframes;
		l_pmean++;
		l_pcabs_start++;
	}
	// calculate the rms
	l_pmean = p_pmean;
	double* l_prms = p_prms;
	l_pcabs_start = l_pcabs0;
	for(n=p_nfftsize;n;--n){
		*l_prms = 0.0;
		long f;
		l_pcabs = l_pcabs_start;
		for(f=p_nframes;f;--f){
			double diff = *l_pcabs - *l_pmean;
			*l_prms += diff*diff;
			l_pcabs += p_nfftsize;
		}
		*l_prms /= p_nframes;
		*l_prms = sqrt(*l_prms);
		l_prms++;
		l_pmean++;
		l_pcabs_start++;
	}
	
	free(l_pcabs0);
}

void spectrum_copy( std::complex<double>* p_psrc, std::complex<double>* p_pdst, long p_npoints )
{
	std::complex<double>* l_psrc = p_psrc;
	std::complex<double>* l_pdst = p_pdst;
	long n;
	for(n=p_npoints;n;--n){
		*p_pdst = *p_psrc;
		p_pdst++;
		p_psrc++;
	}
}

void spectrum_set( std::complex<double>* p_pdst, long p_npoints, std::complex<double> p_val )
{
	std::complex<double>* l_pdst = p_pdst;
	long n;
	for(n=p_npoints;n;--n){
		*p_pdst = p_val;
		p_pdst++;
	}
}
