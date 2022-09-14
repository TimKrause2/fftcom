typedef struct
{
	int m_nss;
	int m_nwindow;
	double  m_fs;
	double  m_fv;
	double  m_tsample;
	double* m_pdata_in;
	int     m_idata_in;
	double* m_pdata_out;
	int     m_idata_outr;
	int     m_idata_outw;
	int     m_ndata_out;
	double* m_ptable;
} vsample_t;

vsample_t* vsample_create( int p_nwindow, int p_nss, int p_fsamplerate );
void vsample_store_virtual( vsample_t* p_pvsample, double p_value );
void vsample_sample_in( vsample_t* p_pvsample, double p_value );
int vsample_count_virtual( vsample_t* p_pvsample );
double vsample_read_virtual( vsample_t* p_pvsample );
