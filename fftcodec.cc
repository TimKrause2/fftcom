#include <math.h>
//#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include "fftcodec.h"
#include "font.h"

FreeTypeFont g_font;

fftcodec::fftcodec( int p_nfftsize, int p_fsamplerate )
{
    m_nfftsize = p_nfftsize;
    m_fsamplerate = p_fsamplerate;
	//
	// codec initialization
	//
    m_codec_mode = CODEC_MODE_IDLE;
    m_codec_flags = 0;
    pthread_mutex_init( &m_codec_mutex, NULL );
	//
	// receiver initialization
	//
    m_rx_mode = RX_MODE_IDLE;
    m_pvsample = new vsample( NWINDOW, NSUPERSAMPLE, p_fsamplerate );
    m_irx = 0;
    m_irxframe = 0;
    m_rx_raw = (double*)malloc( sizeof(double) * p_nfftsize );
    m_rx_x = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * p_nfftsize );
    m_rx_X = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * p_nfftsize );
    m_rx_Xframes = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * p_nfftsize * NFRAMES );
    m_rx_impulse = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * p_nfftsize );
    m_rx_Impulse = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * p_nfftsize );
    m_rx_Sinv = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * (p_nfftsize/2+1) );
    m_rx_Sinvd = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * (p_nfftsize + 1) );
    m_rx_sinv = (double*)malloc( sizeof(double) * p_nfftsize * 2 );
    m_rx_plan = fftw_plan_dft_1d( p_nfftsize,
                            reinterpret_cast<fftw_complex*>(m_rx_x),
                            reinterpret_cast<fftw_complex*>(m_rx_X),
                            FFTW_FORWARD, FFTW_MEASURE );
    m_rx_impulse_plan = fftw_plan_dft_1d( p_nfftsize,
                            reinterpret_cast<fftw_complex*>(m_rx_Impulse),
                            reinterpret_cast<fftw_complex*>(m_rx_impulse),
                            FFTW_BACKWARD, FFTW_MEASURE );
    m_rx_sinv_plan = fftw_plan_dft_c2r_1d(
                p_nfftsize,
                reinterpret_cast<fftw_complex*>(m_rx_Sinv),
                m_rx_sinv,
                FFTW_MEASURE);
    m_rx_Sinvd_plan = fftw_plan_dft_r2c_1d(
                p_nfftsize*2,
                m_rx_sinv,
                reinterpret_cast<fftw_complex*>(m_rx_Sinvd),
                FFTW_MEASURE);
    spectrum_set( m_rx_impulse, p_nfftsize, 0.0);
    bzero( m_rx_sinv, sizeof(double) * p_nfftsize * 2 );
    m_rx_nf_mean = new double[ p_nfftsize ];
    m_rx_nf_rms = new double[ p_nfftsize ];
    m_rx_ir_mean = new double[ p_nfftsize ];
    m_rx_ir_rms = new double[ p_nfftsize ];
    m_rx_Xarg = new double[ p_nfftsize ];
	//
	// transmitter initialization
	//
    m_tx_mode = TX_MODE_SILENCE;
    m_itx = 0;
    m_itxframe = 0;
    m_itxsweep = 0;
    m_tx_x = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * p_nfftsize );
    m_tx_X = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * p_nfftsize );
    m_tx_plan = fftw_plan_dft_1d( p_nfftsize,
                            reinterpret_cast<fftw_complex*>(m_tx_X),
                            reinterpret_cast<fftw_complex*>(m_tx_x),
                            FFTW_BACKWARD, FFTW_MEASURE );
	//
	// renderer initialization
	//
    m_draw_mode = DRAW_MODE_TIME;
    m_draw_src = DRAW_SRC_RX;
    m_draw_flags = DRAW_FLAG_NF|DRAW_FLAG_IR;
    m_iconstellation0 = 0;
    m_iconstellation1 = p_nfftsize-1;
    m_draw_constellation_max = 1.0;
    m_draw_db_min = -180.0;
    m_draw_db_max =  5.0;
    m_draw_time_max = 1.0;
    m_draw_data_x = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * p_nfftsize );
    m_draw_data_X = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * p_nfftsize );
    m_draw_data_impulse = (std::complex<double>*)fftw_malloc( sizeof(std::complex<double>) * p_nfftsize );
}

void fftcodec::store_sample( double p_val )
{
    pthread_mutex_lock( &m_codec_mutex );
    m_pvsample->sample_in( p_val );
    while( m_pvsample->m_ndata_out ){
        double l_vval = m_pvsample->read_virtual();
        m_rx_raw[m_irx] = l_vval;
        if( ++m_irx == m_nfftsize ){
			// received a complete frame
            double* l_psrc = m_rx_raw;
            std::complex<double>* l_pdst = m_rx_x;
			long l_i;
            for( l_i=m_nfftsize; l_i; l_i-- ){
				*(l_pdst++) = *(l_psrc++);
			}
            fftw_execute( m_rx_plan );
            decode_spectrum( );
            m_irx = 0;
		}
	}
    pthread_mutex_unlock( &m_codec_mutex );
}

double fftcodec::read_sample(  )
{
    pthread_mutex_lock( &m_codec_mutex );
    if( m_itx == 0 ){
        encode_spectrum( );
        fftw_execute( m_tx_plan );
	}
    double l_val = std::real(m_tx_x[m_itx]);
    if(++m_itx == m_nfftsize){
        m_itx = 0;
	}
    pthread_mutex_unlock( &m_codec_mutex );
	return l_val;
}

void fftcodec::encode_spectrum(  )
{
    switch( m_tx_mode ){
		case TX_MODE_SILENCE:
            spectrum_set( m_tx_X, m_nfftsize, 0.0 );
			break;

		case TX_MODE_IMPULSE:
            spectrum_set( m_tx_X, m_nfftsize, 1.0/m_nfftsize );
			break;
			
		case TX_MODE_COSINE:
        printf("encoding cosine spectrum at %ld\n",m_itxsweep);
            spectrum_set( m_tx_X, m_nfftsize, 0.0 );
            if( m_itxsweep == 0 ){
                m_tx_X[0] = 1.0;
            }else if( m_itxsweep == m_nfftsize/2 ){
                m_tx_X[ m_nfftsize/2 ] = 1.0;
			}else{
                m_tx_X[m_itxsweep] = 0.5;
                m_tx_X[m_nfftsize-m_itxsweep]=0.5;
			}
			break;
			
		case TX_MODE_SYNC_ONLY:
            spectrum_set( m_tx_X, m_nfftsize, 1.0/m_nfftsize );
			break;
	}
    m_itxframe++;
    update_codec_mode( );
}

void fftcodec::decode_spectrum(  )
{
	long n;
    std::complex<double>* l_pdata=m_rx_X;
    double* l_parg0 = m_rx_Xarg;
    long i=0;
    for(n=m_nfftsize;n;--n){
        *l_pdata /= m_nfftsize;
		*l_parg0 = std::arg( *l_pdata );
		l_pdata++;
		l_parg0++;
        i++;
	}
    l_parg0 = m_rx_Xarg;
	double* l_parg1 = l_parg0+1;
	double  l_argmean = 0.0;
    for(n=m_nfftsize-1;n;--n){
		double l_darg = *l_parg1 - *l_parg0;
		if( l_darg > M_PI ) {
			l_darg = l_darg-2*M_PI;
		}
		if( l_darg <-M_PI ) {
			l_darg = l_darg+2*M_PI;
		}
		l_darg = fabs(l_darg);
		l_argmean += l_darg;
        if( m_codec_flags&CODEC_FLAG_PRINT ){
			printf("arg0:%f arg1:%f darg:%f\n", *l_parg0, *l_parg1, l_darg );
		}
		l_parg0++;
		l_parg1++;
	}
    if( m_codec_flags&CODEC_FLAG_PRINT ){
        m_codec_flags ^= CODEC_FLAG_PRINT;
	}
	
    l_argmean /= m_nfftsize;
    m_rx_Xargmean = l_argmean;
    switch( m_rx_mode ){
    case RX_MODE_IDLE:
        break;
    case RX_MODE_SWEEP_TEST:
        if(m_irxframe == NRXSWEEPFRAMES){
            int l_itxsweep = m_itxsweep;
            std::complex<double> Ii = m_rx_X[l_itxsweep];
            printf("sampling %i :<%lg,%lg> abs:%lg\n", l_itxsweep, std::real(Ii), std::imag(Ii), std::abs(Ii) );
            m_rx_Impulse[l_itxsweep] = Ii;
            m_rx_Impulse[m_nfftsize-l_itxsweep] = std::conj(Ii);
        }
        if(++m_irxframe == NSWEEPFRAMES){
            m_irxframe=0;
        }
        break;

    case RX_MODE_NOISE_FLOOR:
        spectrum_copy(
            m_rx_X,
            &m_rx_Xframes[m_irxframe*m_nfftsize],
            m_nfftsize );
        if(++m_irxframe == NFRAMES){
            spectrum_stats(
                m_rx_Xframes,
                m_rx_nf_mean,
                m_rx_nf_rms,
                m_nfftsize,
                NFRAMES );
            m_rx_mode = RX_MODE_IDLE;
        }
        break;

    case RX_MODE_IMPULSE_RESPONSE:
        spectrum_copy(
            m_rx_X,
            &m_rx_Xframes[m_irxframe*m_nfftsize],
            m_nfftsize );
        if(++m_irxframe == NFRAMES){
            spectrum_stats(
                m_rx_Xframes,
                m_rx_ir_mean,
                m_rx_ir_rms,
                m_nfftsize,
                NFRAMES );
            m_rx_mode = RX_MODE_IDLE;
        }
        break;

    case RX_MODE_SYNC_ONLY:
        sync( );
        break;
	}
    update_codec_mode( );
}

void fftcodec::sync(  )
{
	// estimate the sample error from the pilot tone phase
    double l_error = m_rx_Xarg[1]/2/M_PI*m_nfftsize;
	double l_delta_fv,l_delta_error;
	
    switch( m_sync_mode ){
		case SYNC_MODE_OFFSET_ANALYZE:
			// use the error to adjust the virtual sampler frequency
            l_delta_fv = l_error*m_fsamplerate/m_nfftsize;
            m_pvsample->m_fv += l_delta_fv;
            m_delta_fv = l_delta_fv;
            m_sync_mode = SYNC_MODE_OFFSET_RESET;
			break;
			
		case SYNC_MODE_OFFSET_RESET:
			// revert the virtual sampler frequency back
            m_pvsample->m_fv -= m_delta_fv;
            m_sync_mode = SYNC_MODE_FREQUENCY_ANALYZE;
			break;
			
		case SYNC_MODE_FREQUENCY_ANALYZE:
			// store the error
            m_error = l_error;
            m_sync_mode = SYNC_MODE_FREQUENCY_APPLY;
			break;
			
		case SYNC_MODE_FREQUENCY_APPLY:
			// calculate the error differential
            l_delta_error = l_error - m_error;
            m_pvsample->m_fv += l_delta_error*m_fsamplerate/m_nfftsize;
            m_sync_mode = SYNC_MODE_OFFSET_ANALYZE;
			break;
	}
}

void fftcodec::codec_mode( codec_mode_t p_mode )
{
    pthread_mutex_lock( &m_codec_mutex );
    if( m_codec_flags&CODEC_FLAG_MODE_QUEUED ){
        pthread_mutex_unlock( &m_codec_mutex );
		return;
	}
    m_codec_flags |= CODEC_FLAG_MODE_QUEUED;
    m_codec_mode_next = p_mode;
    pthread_mutex_unlock( &m_codec_mutex );
}


void fftcodec::init_codec_mode(  )
{
    if( (m_codec_flags&CODEC_FLAG_MODE_QUEUED) == 0 )
		return;
    m_codec_flags ^= CODEC_FLAG_MODE_QUEUED;
    switch( m_codec_mode_next ){
		case CODEC_MODE_IDLE:
            init_tx_mode( TX_MODE_SILENCE );
            init_rx_mode( RX_MODE_IDLE );
            m_codec_mode = CODEC_MODE_IDLE;
			break;
		case CODEC_MODE_IMPULSE_TEST:
            init_tx_mode( TX_MODE_IMPULSE );
            init_rx_mode( RX_MODE_IDLE );
            m_codec_mode = CODEC_MODE_IMPULSE_TEST;
			break;
		case CODEC_MODE_SWEEP_TEST:
            init_tx_mode( TX_MODE_COSINE );
            init_rx_mode( RX_MODE_SWEEP_TEST );
            m_codec_mode = CODEC_MODE_SWEEP_TEST;
            m_itxsweep = 0;
			break;
		case CODEC_MODE_NOISE_FLOOR:
		case CODEC_MODE_NOISE_FLOOR_INIT:
            init_tx_mode( TX_MODE_SILENCE );
            init_rx_mode( RX_MODE_IDLE );
            m_codec_mode = CODEC_MODE_NOISE_FLOOR_INIT;
            m_codec_mode_next = CODEC_MODE_IDLE;
            m_codec_flags |= CODEC_FLAG_MODE_QUEUED;
			break;
		case CODEC_MODE_IMPULSE_RESPONSE:
		case CODEC_MODE_IMPULSE_RESPONSE_INIT:
            init_tx_mode( TX_MODE_IMPULSE );
            init_rx_mode( RX_MODE_IDLE );
            m_codec_mode = CODEC_MODE_IMPULSE_RESPONSE_INIT;
            m_codec_mode_next = CODEC_MODE_IDLE;
            m_codec_flags |= CODEC_FLAG_MODE_QUEUED;
			break;
		case CODEC_MODE_SYNC_ONLY:
		case CODEC_MODE_SYNC_ONLY_INIT:
            init_tx_mode( TX_MODE_SYNC_ONLY );
            init_rx_mode( RX_MODE_IDLE );
            m_codec_mode = CODEC_MODE_SYNC_ONLY_INIT;
            m_pvsample->m_fv = m_pvsample->m_fs;
			break;
	}
}

void fftcodec::init_rx_mode( rx_mode_t p_rxmode )
{
    m_rx_mode = p_rxmode;
    m_irxframe = 0;
	if( p_rxmode == RX_MODE_SYNC_ONLY )
        m_sync_mode = SYNC_MODE_OFFSET_ANALYZE;
    if( p_rxmode == RX_MODE_SWEEP_TEST ){
        spectrum_set( m_rx_Impulse, m_nfftsize, 0.0 );
    }
}

void fftcodec::init_tx_mode( tx_mode_t p_txmode )
{
    m_tx_mode = p_txmode;
    m_itxframe = 0;
}

void fftcodec::update_codec_mode(  )
{
    switch( m_codec_mode ){
		case CODEC_MODE_IDLE:
            init_codec_mode();
			break;
		case CODEC_MODE_IMPULSE_TEST:
            init_codec_mode( );
			break;
		case CODEC_MODE_SWEEP_TEST:
            if( m_itxframe == NSWEEPFRAMES ){
                if( m_itxsweep == m_nfftsize/2 ){
                    printf("Sweep test complete. Saving results.\n");
                    fftw_execute( m_rx_impulse_plan );
                    init_tx_mode( TX_MODE_SILENCE );
                    init_rx_mode( RX_MODE_IDLE );
                    m_codec_mode = CODEC_MODE_IDLE;
                    m_itxsweep = 0;
                    FILE* rfile = fopen("response.txt","w");
                    if(rfile){
                        fprintf(rfile, "FFT size:%d Sample rate:%lf\n",
                                m_nfftsize,
                                m_fsamplerate);
                        std::complex<double> *l_pI = m_rx_Impulse;
                        for(int i=0;i<=m_nfftsize/2;i++){
                            fprintf(rfile, "%20.19lg, %20.19lg\n",
                                    std::real(*l_pI),std::imag(*l_pI));
                            l_pI++;
                        }
                        fclose(rfile);
                    }
                    // compute the geometric mean of the impulse spectrum
                    double gm=0;
                    for(int i=1;i<=m_nfftsize/2;i++){
                        gm += log(abs(m_rx_Impulse[i]));
                    }
                    gm /= m_nfftsize/2;
                    gm = exp(gm);
                    // scale the impulse spectrum by the geometric mean
                    for(int i=0;i<=m_nfftsize/2;i++){
                        m_rx_Sinv[i] = std::complex<double>( abs(m_rx_Impulse[i])/gm, 0.0 );
                    }
                    // compute the inverse and cull certain frequencies;
                    for(int i=0;i<=m_nfftsize/2;i++){
                        if(i==0 || i==m_nfftsize/2){
                            m_rx_Sinv[i] = std::complex<double>(1.0, 0.0);
                        }else{
                            double a = abs(m_rx_Sinv[i]);
                            if(a<0.1){
                                m_rx_Sinv[i] = std::complex<double>(1.0, 0.0);
                            }else{
                                m_rx_Sinv[i] = std::complex<double>(1.0/a, 0.0);
                            }
                        }
                    }
                    // write out the single inverse spectrum
                    FILE* ifile = fopen("Sinv.txt","w");
                    if(ifile){
                        fprintf(ifile, "FFTsize:%d\n",
                                m_nfftsize);
                        std::complex<double> *l_pI = m_rx_Sinv;
                        for(int i=0;i<=m_nfftsize/2;i++){
                            fprintf(ifile, "%20.19lg, %20.19lg\n",
                                    std::real(*l_pI),std::imag(*l_pI));
                            l_pI++;
                        }
                        fclose(ifile);
                    }

                    // shift in the time domain half way into the buffer
                    for(int i=0;i<=m_nfftsize/2;i++){
                        std::complex<double> k = exp( std::complex<double>(0.0, (double)i*2.0*M_PI*0.5) );
                        m_rx_Sinv[i] *= k;
                    }
                    // convert to the time domain
                    fftw_execute( m_rx_sinv_plan );
                    // convert back to the double size frequency domain
                    fftw_execute( m_rx_Sinvd_plan );
                    // write out the double inverse spectrum
                    ifile = fopen("Sinvd.txt","w");
                    if(ifile){
                        fprintf(ifile, "FFTsize:%d\n",
                                m_nfftsize);
                        std::complex<double> *l_pI = m_rx_Sinvd;
                        for(int i=0;i<=m_nfftsize;i++){
                            *l_pI /= m_nfftsize * 2;
                            fprintf(ifile, "%20.19lg, %20.19lg\n",
                                    std::real(*l_pI),std::imag(*l_pI));
                            l_pI++;
                        }
                        fclose(ifile);
                    }
				}else{
                    init_tx_mode( TX_MODE_COSINE );
                    m_itxsweep++;
				}
			}
            init_codec_mode();
			break;
		case CODEC_MODE_NOISE_FLOOR_INIT:
            if( m_itxframe == NINITFRAMES ){
                init_rx_mode( RX_MODE_NOISE_FLOOR );
                m_codec_mode = CODEC_MODE_NOISE_FLOOR;
			}
			break;
		case CODEC_MODE_NOISE_FLOOR:
            if( m_rx_mode == RX_MODE_IDLE ){
                m_codec_flags |= CODEC_FLAG_NF_COMPLETE;
                init_codec_mode( );
			}
			break;
		case CODEC_MODE_IMPULSE_RESPONSE_INIT:
            if( m_itxframe == NINITFRAMES ){
                init_rx_mode( RX_MODE_IMPULSE_RESPONSE );
                m_codec_mode = CODEC_MODE_IMPULSE_RESPONSE;
			}
			break;
		case CODEC_MODE_IMPULSE_RESPONSE:
            if( m_rx_mode == RX_MODE_IDLE ){
                m_codec_flags |= CODEC_FLAG_IR_COMPLETE;
                init_codec_mode( );
			}
			break;
		case CODEC_MODE_SYNC_ONLY_INIT:
            if( m_itxframe == NINITFRAMES ){
                init_rx_mode( RX_MODE_SYNC_ONLY );
                m_codec_mode = CODEC_MODE_SYNC_ONLY;
			}
			break;
		case CODEC_MODE_SYNC_ONLY:
            init_codec_mode( );
			break;
	}
}

void fftcodec::draw_init( Display* p_pDisplay )
{
	g_font.LoadOutline("UbuntuMono-R.ttf", 12, Pixel32(0,0,0), Pixel32(128,128,128), 1.5);
    m_text_height = 12;
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

const char* fftcodec::draw_mode_string(  )
{
    switch( m_draw_mode ){
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
    return "";
}

const char* fftcodec::codec_mode_string(  )
{
    switch( m_codec_mode ){
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
    return "";
}

const char* fftcodec::rx_mode_string(  )
{
    switch( m_rx_mode ){
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
    return "";
}

const char* fftcodec::tx_mode_string(  )
{
    switch( m_tx_mode ){
		case TX_MODE_SILENCE:
			return "Silence";
		case TX_MODE_IMPULSE:
			return "Impulse";
		case TX_MODE_COSINE:
			return "Cosine";
		case TX_MODE_SYNC_ONLY:
			return "Sync Only";
	}
    return "";
}

void fftcodec::draw( int p_width, int p_height )
{
	char l_str[1024];
	std::complex<double>* l_psrcx;
	std::complex<double>* l_psrcX;
    switch( m_draw_src ){
		case DRAW_SRC_RX:
            l_psrcx = m_rx_x;
            l_psrcX = m_rx_X;
			break;
		case DRAW_SRC_TX:
            l_psrcx = m_tx_x;
            l_psrcX = m_tx_X;
			break;
	}
    pthread_mutex_lock( &m_codec_mutex );
    spectrum_copy( l_psrcx, m_draw_data_x, m_nfftsize );
    spectrum_copy( l_psrcX, m_draw_data_X, m_nfftsize );
    spectrum_copy( m_rx_impulse, m_draw_data_impulse, m_nfftsize );
    pthread_mutex_unlock( &m_codec_mutex );
    l_psrcx = m_draw_data_x;
    l_psrcX = m_draw_data_X;
	std::complex<double>* l_pval;
	long n,i,N;
	GLdouble l_color[3];
    long l_npoints = m_nfftsize/2 + 1;
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
    switch( m_draw_mode ){
    case DRAW_MODE_TIME:
        glMatrixMode( GL_PROJECTION );
        glLoadIdentity();
        gluOrtho2D( 0.0, 1.0, -m_draw_time_max, m_draw_time_max );
        glMatrixMode( GL_MODELVIEW );
        glLoadIdentity();
        l_pval = l_psrcx;
        glColor4d( 0.5, 1.0, 0.5, 1.0 );
        glBegin( GL_LINE_STRIP );
        for(n=0;n<m_nfftsize;n++){
            double l_x = (double)n/(m_nfftsize-1);
            double l_y = std::real(*l_pval);
            glVertex2d( l_x, l_y );
            l_pval++;
        }
        glEnd();
        glColor3d(1.0,1.0,1.0);
        g_font.Printf( 0, p_height - m_text_height, "Scale:%7.5f", m_draw_time_max );
        break;

    case DRAW_MODE_IMPULSE:
        glMatrixMode( GL_PROJECTION );
        glLoadIdentity();
        gluOrtho2D( 0.0, 1.0, -m_draw_time_max, m_draw_time_max );
        glMatrixMode( GL_MODELVIEW );
        glLoadIdentity();
        l_pval = m_draw_data_impulse;
        glColor4d( 0.5, 1.0, 0.5, 0.5 );
        glBegin( GL_LINE_STRIP );
        for(n=0;n<m_nfftsize;n++){
            double l_x = (double)n/(m_nfftsize-1);
            double l_y = std::real(*l_pval);
            glVertex2d( l_x, l_y );
            l_pval++;
        }
        glEnd();
        glColor3d(1.0,1.0,1.0);
// 			sprintf( l_str, "Scale:%7.5f", m_draw_time_max );
// 			glWindowPos2i( (GLint)0, (GLint)0 );
// 			glListBase( g_font_base );
// 			glCallLists( strlen(l_str), GL_UNSIGNED_BYTE, (GLubyte*) l_str );
        g_font.Printf( 0, p_height - m_text_height, "Scale:%7.5f", m_draw_time_max );
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
        if( m_codec_flags&CODEC_FLAG_NF_COMPLETE ){
            if( m_draw_flags&DRAW_FLAG_NF ){
                // draw the noise floor
                glColor3d( 1.0, 0.5, 1.0 );
                dspectrum_draw_abs( m_rx_nf_mean, l_npoints );
                dspectrum_draw_abs( m_rx_nf_rms, l_npoints );
            }
        }
        // draw the impulse response
        if( m_codec_flags&CODEC_FLAG_IR_COMPLETE ){
            if( m_draw_flags&DRAW_FLAG_IR ){
                // draw the impulse response
                glColor3d( 0.5, 1.0, 1.0 );
                dspectrum_draw_abs( m_rx_ir_mean, l_npoints );
                dspectrum_draw_abs( m_rx_ir_rms, l_npoints );
            }
        }
        break;

    case DRAW_MODE_FREQUENCY_ABS_LOG:
        glMatrixMode( GL_PROJECTION );
        glLoadIdentity();
        gluOrtho2D( 0.0, 1.0, m_draw_db_min, m_draw_db_max );
        glMatrixMode( GL_MODELVIEW );
        glLoadIdentity();
        spectrum_draw_log( l_psrcX, l_npoints );
        // draw the noise floor
        if( m_codec_flags&CODEC_FLAG_NF_COMPLETE ){
            if( m_draw_flags&DRAW_FLAG_NF ){
                // draw the noise floor
                glColor3d( 1.0, 0.5, 1.0 );
                dspectrum_draw_log_bounds( m_rx_nf_mean, m_rx_nf_rms, l_npoints );
            }
        }
        // draw the impulse response
        if( m_codec_flags&CODEC_FLAG_IR_COMPLETE ){
            if( m_draw_flags&DRAW_FLAG_IR ){
                // draw the impulse response
                glColor3d( 0.5, 1.0, 1.0 );
                dspectrum_draw_log_bounds( m_rx_ir_mean, m_rx_ir_rms, l_npoints );
            }
        }
        glColor3d(1.0,1.0,1.0);
// 			sprintf( l_str, "dB max:%5.1f dB min:%5.1f", m_draw_db_max, m_draw_db_min );
// 			glWindowPos2i( (GLint)0, (GLint)0 );
// 			glListBase( g_font_base );
// 			glCallLists( strlen(l_str), GL_UNSIGNED_BYTE, (GLubyte*) l_str );
        g_font.Printf( 0, p_height - m_text_height,
                                    "dB max:%5.1f dB min:%5.1f", m_draw_db_max,
                                    m_draw_db_min );
        break;

    case DRAW_MODE_FREQUENCY_CONSTELLATION:
        glEnable( GL_POINT_SMOOTH );
        glPointSize( 3.0 );
        glMatrixMode( GL_PROJECTION );
        glLoadIdentity();
        gluOrtho2D(
            -m_draw_constellation_max*p_width/p_height,
            m_draw_constellation_max*p_width/p_height,
            -m_draw_constellation_max,
            m_draw_constellation_max );
        glMatrixMode( GL_MODELVIEW );
        glLoadIdentity();
        l_pval = &l_psrcX[m_iconstellation0];
        N = m_iconstellation1 - m_iconstellation0 + 1;
        glBegin( GL_POINTS );
        for(n=m_iconstellation0,i=0;n<=m_iconstellation1;n++,i++){
            double l_hue = (double)i/N;
            color_set_hsv( l_color, l_hue, 0.75, 1.0 );
            glColor3dv( l_color );
            glVertex2d( std::real(*l_pval), std::imag(*l_pval) );
            l_pval++;
        }
        glEnd();
        glColor3d(1.0,1.0,1.0);
// 			sprintf( l_str, "Scale:%8.6f MRA:%+8.6f", m_draw_constellation_max, m_rx_Xargmean );
// 			glWindowPos2i( (GLint)0, (GLint)0 );
// 			glListBase( g_font_base );
// 			glCallLists( strlen(l_str), GL_UNSIGNED_BYTE, (GLubyte*) l_str );
        g_font.Printf( 0, p_height - m_text_height, "Scale:%8.6f MRA:%+8.6f",
                                    m_draw_constellation_max,
                           m_rx_Xargmean );
        break;
	}
	// draw the title information
	glColor3d(1.0,1.0,1.0);
// 	sprintf( l_str, "Modes Draw:%s CODEC:%s RX:%s TX:%s",
// 					 fftcodec_draw_mode_string( p_pfftcodec ),
// 					 fftcodec_codec_mode_string( p_pfftcodec ),
// 					 fftcodec_rx_mode_string( p_pfftcodec ),
// 					 fftcodec_tx_mode_string( p_pfftcodec ) );
// 	glWindowPos2i( (GLint)0, (GLint)(p_height - m_text_height) );
// 	glListBase( g_font_base );
// 	glCallLists( strlen(l_str), GL_UNSIGNED_BYTE, (GLubyte*) l_str );
    g_font.Printf( 0, m_text_height,
        "Modes Draw:%s CODEC:%s RX:%s TX:%s",
        draw_mode_string( ),
        codec_mode_string( ),
        rx_mode_string( ),
        tx_mode_string( ) );
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
