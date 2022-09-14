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

vsample_t* vsample_create( int p_nwindow, int p_nss, int p_fsamplerate )
{
	vsample_t* l_pvsample;
	l_pvsample = (vsample_t*)malloc(sizeof(vsample_t));
	if(!l_pvsample){
		return NULL;
	}
	if(p_nwindow&1)p_nwindow++;
	l_pvsample->m_nwindow = p_nwindow;
	l_pvsample->m_nss = p_nss;
	l_pvsample->m_fs = p_fsamplerate;
	l_pvsample->m_fv = p_fsamplerate;
	l_pvsample->m_tsample = (double)(p_nwindow/2);
	l_pvsample->m_pdata_in = (double*)malloc(p_nwindow*sizeof(double));
	l_pvsample->m_idata_in = 0;
	l_pvsample->m_pdata_out = (double*)malloc(p_nwindow*sizeof(double));
	l_pvsample->m_idata_outr = 0;
	l_pvsample->m_idata_outw = 0;
	l_pvsample->m_ndata_out = 0;
	l_pvsample->m_ptable = (double*)malloc(p_nwindow*p_nss*sizeof(double));
	int i;
	int s;
	double* l_pdata = l_pvsample->m_ptable;
	for(i=0;i<p_nss;i++){
		double l_deltat = (double)i/p_nss;
		for(s=0;s<p_nwindow;s++){
			double l_n = s-(p_nwindow/2-1);
			*l_pdata = sinc( (l_n-l_deltat)*M_PI );
			l_pdata++;
		}
	}
	return l_pvsample;
}

void vsample_sample_in( vsample_t* p_pvsample, double p_value )
{
	// write the data into the delay line buffer
	p_pvsample->m_pdata_in[p_pvsample->m_idata_in] = p_value;
	int l_ioldest = p_pvsample->m_idata_in + 1;
	if(l_ioldest == p_pvsample->m_nwindow)
		l_ioldest = 0;
	
	// render samples within the convolver range
	while( p_pvsample->m_tsample >= 0.0 && p_pvsample->m_tsample< 1.0 ){
		// determine the filter to use
		int l_nfilter = floor(p_pvsample->m_tsample*p_pvsample->m_nss);
		
		// correlate with impulse sample
		double* l_pdata = &p_pvsample->m_pdata_in[l_ioldest];
		double* l_ptable = &p_pvsample->m_ptable[l_nfilter*p_pvsample->m_nwindow];
		double l_acc = 0.0;
		int l_nloop = p_pvsample->m_nwindow - l_ioldest;
		int n;
		for(n=l_nloop;n;--n){
			l_acc += *l_pdata * *l_ptable;
			l_pdata++;
			l_ptable++;
		}
		l_nloop = p_pvsample->m_nwindow - l_nloop;
		l_pdata = p_pvsample->m_pdata_in;
		for(n=l_nloop;n;--n){
			l_acc += *l_pdata * *l_ptable;
			l_pdata++;
			l_ptable++;
		}
		vsample_store_virtual(p_pvsample,l_acc);
		// set the virtual sample time to the next sample
		p_pvsample->m_tsample += p_pvsample->m_fs / p_pvsample->m_fv;
	}

	// update the time for the virtual sample
	p_pvsample->m_tsample -= 1.0;

	// update the write index
	p_pvsample->m_idata_in = l_ioldest;
}

void vsample_store_virtual( vsample_t* p_pvsample, double p_value )
{
	if(p_pvsample->m_ndata_out==p_pvsample->m_nwindow)
		return;
	p_pvsample->m_pdata_out[p_pvsample->m_idata_outw] = p_value;
	if(++p_pvsample->m_idata_outw==p_pvsample->m_nwindow)
		p_pvsample->m_idata_outw = 0;
	p_pvsample->m_ndata_out++;
}

double vsample_read_virtual( vsample_t* p_pvsample )
{
	if(p_pvsample->m_ndata_out==0)
		return 0.0;
	double l_r=p_pvsample->m_pdata_out[p_pvsample->m_idata_outr];
	if(++p_pvsample->m_idata_outr==p_pvsample->m_nwindow)
		p_pvsample->m_idata_outr = 0;
	p_pvsample->m_ndata_out--;
	return l_r;
}

