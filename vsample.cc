#include <math.h>
#include <stdlib.h>
#include "vsample.h"

double sinc(double p_x)
{
	if(fabs(p_x)<1e-9)
		return 1.0;
	else
		return sin(p_x)/p_x;
}

vsample::vsample( int p_nwindow, int p_nss, int p_fsamplerate )
{
	if(p_nwindow&1)p_nwindow++;
    m_nwindow = p_nwindow;
    m_nss = p_nss;
    m_fs = p_fsamplerate;
    m_fv = p_fsamplerate;
    m_tsample = (double)(p_nwindow/2);
    m_pdata_in = new double[p_nwindow];
    m_idata_in = 0;
    m_pdata_out = new double[p_nwindow];
    m_idata_outr = 0;
    m_idata_outw = 0;
    m_ndata_out = 0;
    m_ptable = new double[p_nwindow*p_nss];
	int i;
	int s;
    double* l_pdata = m_ptable;
	for(i=0;i<p_nss;i++){
		double l_deltat = (double)i/p_nss;
		for(s=0;s<p_nwindow;s++){
			double l_n = s-(p_nwindow/2-1);
			*l_pdata = sinc( (l_n-l_deltat)*M_PI );
			l_pdata++;
		}
	}
}

void vsample::sample_in( double p_value )
{
	// write the data into the delay line buffer
    m_pdata_in[m_idata_in] = p_value;
    int l_ioldest = m_idata_in + 1;
    if(l_ioldest == m_nwindow)
		l_ioldest = 0;
	
	// render samples within the convolver range
    while( m_tsample >= 0.0 && m_tsample< 1.0 ){
		// determine the filter to use
        int l_nfilter = floor(m_tsample*m_nss);
		
		// correlate with impulse sample
        double* l_pdata = &m_pdata_in[l_ioldest];
        double* l_ptable = &m_ptable[l_nfilter*m_nwindow];
		double l_acc = 0.0;
        int l_nloop = m_nwindow - l_ioldest;
		int n;
		for(n=l_nloop;n;--n){
			l_acc += *l_pdata * *l_ptable;
			l_pdata++;
			l_ptable++;
		}
        l_nloop = m_nwindow - l_nloop;
        l_pdata = m_pdata_in;
		for(n=l_nloop;n;--n){
			l_acc += *l_pdata * *l_ptable;
			l_pdata++;
			l_ptable++;
		}
        store_virtual(l_acc);
		// set the virtual sample time to the next sample
        m_tsample += m_fs / m_fv;
	}

	// update the time for the virtual sample
    m_tsample -= 1.0;

	// update the write index
    m_idata_in = l_ioldest;
}

void vsample::store_virtual( double p_value )
{
    if(m_ndata_out==m_nwindow)
		return;
    m_pdata_out[m_idata_outw] = p_value;
    if(++m_idata_outw==m_nwindow)
        m_idata_outw = 0;
    m_ndata_out++;
}

double vsample::read_virtual( )
{
    if(m_ndata_out==0)
		return 0.0;
    double l_r=m_pdata_out[m_idata_outr];
    if(++m_idata_outr==m_nwindow)
        m_idata_outr = 0;
    m_ndata_out--;
	return l_r;
}

