#include <complex>
#include <fftw3.h>
#include <pthread.h>
#include <GL/glew.h>
#include <X11/Xlib.h>
#include <glm/glm.hpp>

#include "vsample.h"
#include <lgraph.h>
#include <esfont.h>
#include <pgraph.h>

#define NWINDOW 256
#define NSUPERSAMPLE 4096
#define NFRAMES 1024
#define NSWEEPFRAMES 5
#define NRXSWEEPFRAMES 4
#define NINITFRAMES 8

#define HUE_TIME1 155.5
#define SAT_TIME1 1.0
#define VAL_TIME1 1.0
#define HUE_TIME0 155.5
#define SAT_TIME0 0.580
#define VAL_TIME0 0.602

#define HUE_FREQ1 0.0
#define SAT_FREQ1 0.0
#define VAL_FREQ1 1.0
#define HUE_FREQ0 0.0
#define SAT_FREQ0 0.0
#define VAL_FREQ0 0.455

#define HUE_NF1 45.7
#define SAT_NF1 0.534
#define VAL_NF1 1.0
#define HUE_NF0 45.7
#define SAT_NF0 0.534
#define VAL_NF0 0.450

#define HUE_IR1 189.4
#define SAT_IR1 0.545
#define VAL_IR1 1.0
#define HUE_IR0 189.4
#define SAT_IR0 0.545
#define VAL_IR0 0.450

typedef enum {
	CODEC_MODE_IDLE,
	CODEC_MODE_IMPULSE_TEST,
	CODEC_MODE_SWEEP_TEST,
	CODEC_MODE_NOISE_FLOOR_INIT,
	CODEC_MODE_NOISE_FLOOR,
	CODEC_MODE_IMPULSE_RESPONSE_INIT,
	CODEC_MODE_IMPULSE_RESPONSE,
	CODEC_MODE_SYNC_ONLY_INIT,
	CODEC_MODE_SYNC_ONLY,
	CODEC_NMODES
} codec_mode_t;

typedef enum {
	RX_MODE_IDLE,
    RX_MODE_SWEEP_TEST,
	RX_MODE_NOISE_FLOOR,
	RX_MODE_IMPULSE_RESPONSE,
	RX_MODE_SYNC_ONLY,
	RX_NMODES
} rx_mode_t;

typedef enum {
	SYNC_MODE_OFFSET_ANALYZE,
	SYNC_MODE_OFFSET_RESET,
	SYNC_MODE_FREQUENCY_ANALYZE,
	SYNC_MODE_FREQUENCY_APPLY,
	SYNC_NMODES
} sync_mode_t;

#define CODEC_FLAG_MODE_QUEUED    0x01
#define CODEC_FLAG_NF_COMPLETE    0x02
#define CODEC_FLAG_IR_COMPLETE    0x04
#define CODEC_FLAG_SWEEP_COMPLETE 0x08
#define CODEC_FLAG_PRINT          0X10

typedef enum {
	TX_MODE_SILENCE,
	TX_MODE_IMPULSE,
	TX_MODE_COSINE,
	TX_MODE_SYNC_ONLY,
	TX_NMODES
} tx_mode_t;

typedef enum {
	DRAW_MODE_TIME,
    DRAW_MODE_IMPULSE,
	DRAW_MODE_FREQUENCY_ABS_LINEAR,
	DRAW_MODE_FREQUENCY_ABS_LOG,
	DRAW_MODE_FREQUENCY_CONSTELLATION,
} draw_mode_t;

typedef enum {
	DRAW_SRC_RX,
    DRAW_SRC_TX,
    DRAW_SRC_IMPULSE
} draw_src_t;

#define DRAW_FLAG_NF 0x01
#define DRAW_FLAG_IR 0x02

class fftcodec {
public:
	unsigned int m_nfftsize;
    unsigned int m_nspecsize;
	double m_fsamplerate;
	//
	// codec
	//
	codec_mode_t    m_codec_mode;
	long            m_codec_flags;
	codec_mode_t    m_codec_mode_next;
	pthread_mutex_t m_codec_mutex;
	//
	// reciever
	//
	rx_mode_t     m_rx_mode;
	sync_mode_t   m_sync_mode;
	double        m_delta_fv;
	double        m_error;
    vsample*      m_pvsample;
	long          m_irx;
	long          m_irxframe;
	double*       m_rx_raw;
	std::complex<double>* m_rx_x;
	std::complex<double>* m_rx_X;
	std::complex<double>* m_rx_Xframes;
    std::complex<double>* m_rx_impulse;
    std::complex<double>* m_rx_Impulse;
    std::complex<double>* m_rx_Sinv;
    double*               m_rx_sinv;
    std::complex<double>* m_rx_Sinvd;
	fftw_plan     m_rx_plan;
    fftw_plan     m_rx_impulse_plan;
    fftw_plan     m_rx_sinv_plan;
    fftw_plan     m_rx_Sinvd_plan;
	double*       m_rx_nf_mean;
	double*       m_rx_nf_rms;
	double*       m_rx_ir_mean;
	double*       m_rx_ir_rms;
	double*       m_rx_Xarg;
	double        m_rx_Xargmean;
	//
	// transmitter
	//
	tx_mode_t     m_tx_mode;
	long          m_itx;
	long          m_itxframe;
	long          m_itxsweep;
	std::complex<double>* m_tx_x;
	std::complex<double>* m_tx_X;
	fftw_plan     m_tx_plan;
	//
	// renderer
	//
    FreeTypeFont* m_font;
    LGraph*       m_timegraph;
    LGraph*       m_specgraph;
    PGraph*       m_pgraph;
    int           m_text_height;
	draw_mode_t   m_draw_mode;
	draw_src_t    m_draw_src;
	long          m_draw_flags;
	long          m_iconstellation0;
	long          m_iconstellation1;
	double        m_draw_constellation_max;
	double        m_draw_db_min;
	double        m_draw_db_max;
	double        m_draw_time_max;
    float         m_draw_abs_max;
	std::complex<double>* m_draw_data_x;
	std::complex<double>* m_draw_data_X;
    std::complex<double>* m_draw_data_impulse;
    float* m_graphy;
    glm::vec4  m_timeColor1;
    glm::vec4  m_timeColor0;
    glm::vec4  m_freqColor1;
    glm::vec4  m_freqColor0;
    glm::vec4  m_nfColor1;
    glm::vec4  m_nfColor0;
    glm::vec4  m_irColor1;
    glm::vec4  m_irColor0;

public:
    fftcodec(int nfftsize, int fsamplerate);
    ~fftcodec();

    void store_sample( double p_val );
    double read_sample();
    void sync();
    void codec_mode( codec_mode_t p_mode );
private:
    void encode_spectrum();
    void decode_spectrum();
    void update_codec_mode();
    void init_codec_mode();
    void init_rx_mode( rx_mode_t p_rxmode );
    void init_tx_mode( tx_mode_t p_txmode );
    const char* draw_mode_string();
    const char* codec_mode_string();
    const char* rx_mode_string();
    const char* tx_mode_string();

    void spectrum_draw_abs( std::complex<double>* src );
    void spectrum_draw_log( std::complex<double>* src );
    void dspectrum_draw_abs( double* src );
    void dspectrum_draw_log( double* src );
    void dspectrum_draw_log_bounds( double* p_px, double* p_prms );

public:
    void draw_init( Display* p_pDisplay );
    void draw( int p_width, int p_height );

};

void color_set( GLdouble* p_v, double p_r, double p_g, double p_b );
void color_set_hsv( GLdouble* p_c, double p_h, double p_s, double p_v );
void spectrum_stats( std::complex<double>* src, double* p_pmean, double* p_prms, long p_nfftsize, long p_nframes );
void spectrum_copy( std::complex<double>* src, std::complex<double>* dst, long p_npoints );
void spectrum_set( std::complex<double>* dst, long p_npoints, std::complex<double> p_val );

