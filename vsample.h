class vsample
{
public:
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

    vsample(int p_nwindow, int p_nss, int p_fsamplerate);
    ~vsample();
    void store_virtual( double p_value );
    void sample_in( double p_value );
    double read_virtual( );

};

