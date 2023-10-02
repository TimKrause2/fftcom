#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <pthread.h>
#include <signal.h>
#include <getopt.h>
#include <jack/jack.h>
#include <jack/ringbuffer.h>

#include <semaphore.h>
#include <math.h>

#include <X11/Xlib.h>
#include <X11/X.h>
#include <X11/keysym.h>
//#include <GL/gl.h>
//#include <GL/glext.h>
#include <GL/glew.h>
#include <GL/glx.h>

//#include <complex.h>
#include <fftw3.h>

#include "fftcodec.h"

#define FFTSIZE 2400

jack_client_t* g_pclient;

jack_port_t*       g_pportrx;
jack_ringbuffer_t* g_prbrx;
jack_port_t*       g_pporttx;
jack_ringbuffer_t* g_prbtx;
sem_t              g_sem_process;
pthread_t          g_pthread_process;
fftcodec*          g_pfftcodec;

Display* g_pDisplay;
Window   g_window;
GLXContext g_context;
int g_width=400;
int g_height=400;
Atom g_wmDeleteMessage;
Atom g_wmProtocols;

bool g_process = true;
int g_samplerate;
jack_nframes_t g_buffersize;

void synthesize_tx_data( void ){
	int l_samplesize = sizeof(jack_default_audio_sample_t);
	while( jack_ringbuffer_write_space(g_prbtx) >= l_samplesize ){
        jack_default_audio_sample_t l_x = g_pfftcodec->read_sample( );
		jack_ringbuffer_write(g_prbtx,(char*)&l_x,l_samplesize);
	}
}

void* process_thread(void* p_arg)
{
	printf("process_thread\n");
	pthread_setcanceltype( PTHREAD_CANCEL_ASYNCHRONOUS, NULL );
	int l_samplesize = sizeof(jack_default_audio_sample_t);
	while(g_process){
		// wait for signal
		sem_wait(&g_sem_process);
		// synthesize output data
		synthesize_tx_data();
		// process input data
		while( jack_ringbuffer_read_space(g_prbrx) >= l_samplesize ){
			jack_default_audio_sample_t l_x;
			jack_ringbuffer_read(g_prbrx,(char*)&l_x,l_samplesize);
            g_pfftcodec->store_sample( l_x );
		}
	}
	printf("process_thread exiting\n");
	return NULL;
}

void ConfigureEvaluate( XConfigureEvent* p_event )
{
	g_width = p_event->width;
	g_height = p_event->height;
}


void start_thread(void)
{
	sem_init(&g_sem_process,0,0);
	g_process=1;
	pthread_create(&g_pthread_process,NULL,process_thread,NULL);
}

int process( jack_nframes_t nframes, void* p_arg)
{
	jack_default_audio_sample_t* l_prxdata;
	jack_default_audio_sample_t* l_ptxdata;
	int nframebytes = nframes*sizeof(jack_default_audio_sample_t);
	
	// transfer input port data to rx ringbuffer
	l_prxdata = (jack_default_audio_sample_t*)jack_port_get_buffer(g_pportrx,nframes);
	if( jack_ringbuffer_write_space(g_prbrx) >= nframebytes ){
		jack_ringbuffer_write(g_prbrx,(const char*)l_prxdata,nframebytes);
	}
	// transfer tx ringbuffer data to output port
	l_ptxdata = (jack_default_audio_sample_t*)jack_port_get_buffer(g_pporttx,nframes);
	if( jack_ringbuffer_read_space(g_prbtx) >= nframebytes ){
		jack_ringbuffer_read(g_prbtx,(char*)l_ptxdata,nframebytes);
	}
	// signal the process thread
	int l_sval;
	sem_getvalue(&g_sem_process,&l_sval);
	if( l_sval == 0 ){
		sem_post(&g_sem_process);
	}
	
	return 0;
}

void signal_handler(int sig)
{
	g_process=0;
}

void glWindowCreate( int argc, char** argv )
{
	PFNGLXSWAPINTERVALSGIPROC swap_interval;
	Window l_root;
	int l_screen_number;
	Screen *l_screen;
	XVisualInfo *g_visinfo;
	l_screen_number = DefaultScreen( g_pDisplay );
	l_root = RootWindow( g_pDisplay, l_screen_number );
	l_screen = XScreenOfDisplay( g_pDisplay, l_screen_number );
    int l_screen_width = WidthOfScreen(l_screen);
    int l_screen_height = HeightOfScreen(l_screen);
    g_width = l_screen_width/2;
    g_height = l_screen_height/2;
	int l_attribs[]={
		GLX_RGBA,
		GLX_DOUBLEBUFFER,
		GLX_RED_SIZE,8,
		GLX_GREEN_SIZE,8,
		GLX_BLUE_SIZE,8,
		GLX_ALPHA_SIZE,8,
		GLX_DEPTH_SIZE,24,
		GLX_X_VISUAL_TYPE,GLX_DIRECT_COLOR,
		None};
	g_visinfo = glXChooseVisual( g_pDisplay, l_screen_number, l_attribs );
	if(!g_visinfo){
		printf( "Could not find a suitable visual.\n" );
		XCloseDisplay(g_pDisplay);
		exit(0);
	}else{
		printf( "Selected visual 0x%04lX\r\n", (unsigned long)g_visinfo->visualid );
// 		printf( "\tscreen:%d\n", g_visinfo->screen );
// 		printf( "\tdepth:%d\n", g_visinfo->depth );
// 		printf( "\tclass:%d\n", g_visinfo->class );
// 		printf( "\tred_mask  :0x%016lX\n", g_visinfo->red_mask );
// 		printf( "\tgreen_mask:0x%016lX\n", g_visinfo->green_mask );
// 		printf( "\tblue_mask :0x%016lX\n", g_visinfo->blue_mask );
// 		printf( "\tcolormap_size:%d\n", g_visinfo->colormap_size );
// 		printf( "\tbits_per_rgb:%d\n", g_visinfo->bits_per_rgb );
	}
	XSetWindowAttributes l_attr;
	unsigned long l_mask;
	l_attr.background_pixel = 0xFFFFFF;
	l_attr.border_pixel = 0;
	l_attr.colormap = XCreateColormap( g_pDisplay, l_root, g_visinfo->visual, AllocNone );
	l_mask = CWBorderPixel | CWColormap;
	g_window = XCreateWindow(
		g_pDisplay,
		l_root,
        l_screen_width/4,l_screen_height/4,
		g_width,g_height,
        4,
		g_visinfo->depth,
		InputOutput,
		g_visinfo->visual,
		l_mask, &l_attr );
	
	XSelectInput( g_pDisplay, g_window,
                  KeyPressMask |
                  KeyReleaseMask |
                  ButtonPressMask |
                  ButtonReleaseMask |
                  StructureNotifyMask |
                  EnterWindowMask |
                  LeaveWindowMask |
                  PointerMotionMask );
	g_context = glXCreateContext( g_pDisplay, g_visinfo, NULL, True );
	if(!g_context){
		printf("Could not create the GL context.\n");
		XCloseDisplay(g_pDisplay);
		exit(0);
	}
	glXMakeCurrent( g_pDisplay, g_window, g_context );
	swap_interval = (PFNGLXSWAPINTERVALSGIPROC) glXGetProcAddress((unsigned char *)"glXSwapIntervalSGI");
    if(swap_interval){
        swap_interval(1);
    }
	XFree( g_visinfo );
    g_wmDeleteMessage = XInternAtom( g_pDisplay, "WM_DELETE_WINDOW", False);
    g_wmProtocols = XInternAtom( g_pDisplay, "WM_PROTOCOLS", False);
	XSetWMProtocols( g_pDisplay, g_window, &g_wmDeleteMessage, 1);

	XIconSize *l_icon_size_list;
	int l_icon_size_count;
	int l_icon_width;
	int l_icon_height;
	if( XGetIconSizes( g_pDisplay, l_root, &l_icon_size_list, &l_icon_size_count ) == 0 ){
		printf("Window manager doesn't specify icon sizes.\n");
        l_icon_width = 64;
        l_icon_height = 64;
	}else{
		printf("Window manager specifies the following icon sizes.\n");
		int c;
		for(c=0;c<l_icon_size_count;c++){
			printf("min_width:%d min_height:%d max_width:%d max_height:%d width_inc:%d height_inc:%d\n",
						 l_icon_size_list[c].min_width,
					l_icon_size_list[c].min_height,
					l_icon_size_list[c].max_width,
					l_icon_size_list[c].max_height,
					l_icon_size_list[c].width_inc,
					l_icon_size_list[c].height_inc );
		}
		if( l_icon_size_count ){
			l_icon_width = l_icon_size_list[0].max_width;
			l_icon_height = l_icon_size_list[0].max_height;
		}else{
			l_icon_width = 48;
			l_icon_height = 48;
		}
	}
	printf("icon size: %dx%d\n",l_icon_width,l_icon_height);

	Pixmap l_icon_pixmap = XCreatePixmap( g_pDisplay, g_window,
                                        l_icon_width, l_icon_height,
                                        DefaultDepthOfScreen(l_screen) );
	XGCValues l_gcvalues;
	l_gcvalues.foreground = 0x0000FF00;
	l_gcvalues.background = 0x00000000;
	l_gcvalues.line_width = 2;
	GC l_gc = XCreateGC( g_pDisplay, l_icon_pixmap,
										GCForeground | GCBackground | GCLineWidth,
									 &l_gcvalues );
	
	XSetForeground( g_pDisplay, l_gc, 0x00000000 );

	XFillRectangle( g_pDisplay, l_icon_pixmap, l_gc,
									0, 0,
								 l_icon_width, l_icon_height );
	
	XSetForeground( g_pDisplay, l_gc, 0x0000FF00 );
	
	XDrawLine( g_pDisplay, l_icon_pixmap, l_gc,
						 0, l_icon_height/2, l_icon_width, l_icon_height/2 );


	XSizeHints    *l_size_hints;
	XWMHints      *l_wm_hints;
	XClassHint    *l_class_hint;
	XTextProperty  l_tp_window_name;
    char          *l_cp_window_name=(char*)"FFTCOM:";
	XTextProperty  l_tp_icon_name;
	char          *l_cp_icon_name=(char*)"fftcom";

	l_size_hints = XAllocSizeHints();
	if(!l_size_hints){
		printf("Couldn't allocate size hints.\n");
		exit(0);
	}
	l_wm_hints = XAllocWMHints();
	if(!l_wm_hints){
		printf("Couldn't allocate window manager hints.\n");
		exit(0);
	}
	l_class_hint = XAllocClassHint();
	if(!l_class_hint){
		printf("Couldn't allocate class hint.\n");
		exit(0);
	}
	if( XStringListToTextProperty( &l_cp_window_name, 1, &l_tp_window_name ) == 0 ){
		printf("Couldn't allocate text property for window name.\n");
		exit(0);
	}
	if( XStringListToTextProperty( &l_cp_icon_name, 1, &l_tp_icon_name ) == 0 ){
		printf("Couldn't allocate text property for icon name.\n");
		exit(0);
	}
	l_size_hints->flags = PPosition | PSize | PMinSize;
	l_size_hints->min_width = 200;
	l_size_hints->min_height = 200;
	l_size_hints->x = 0;
	l_size_hints->y = 0;
	l_size_hints->width = 300;
	l_size_hints->height = 300;
	l_wm_hints->flags = InputHint | StateHint | IconPixmapHint;
	l_wm_hints->input = True;
	l_wm_hints->initial_state = NormalState;
	l_wm_hints->icon_pixmap = l_icon_pixmap;
	l_class_hint->res_name = (char*)"fftcom";
	l_class_hint->res_class = (char*)"fftcom";
	XSetWMProperties( g_pDisplay, g_window, &l_tp_window_name, &l_tp_icon_name,
                      argv, argc, l_size_hints, l_wm_hints, l_class_hint );

    XMapWindow( g_pDisplay, g_window );
	XFree( l_size_hints );
	XFree( l_wm_hints );
	XFree( l_class_hint );
}

void glWindowDestroy( void )
{
	glXMakeCurrent( g_pDisplay, None, NULL );
	glXDestroyContext( g_pDisplay, g_context );
	XDestroyWindow( g_pDisplay, g_window );
}

int main(int p_narg, char** p_argv)
{
	jack_status_t  l_status;

	XInitThreads();
	//
	// fftw initialization
	//
	int l_result = fftw_init_threads();
	if( l_result ){
		printf("Setting fftw to 4 threads\n");
		fftw_plan_with_nthreads(4);
	}
	
	// open the display
    g_pDisplay = XOpenDisplay(":0");
	if(!g_pDisplay){
		printf("Couldn't open the display\n");
		exit(0);
	}

	glWindowCreate( p_narg, p_argv );
	
    GLenum err = glewInit();
    if(GLEW_OK != err){
        printf("glewInit Error:%s\n", glewGetErrorString(err));
        return 1;
    }

	// become a jack client
	g_pclient = jack_client_open(
		"fftcom",
		JackNoStartServer,
		&l_status);
	if(!g_pclient){
		printf("Couldn't become a client l_status:0x%08X\r\n",l_status);
		exit(0);
	}

    char window_name[1024];
    snprintf(window_name, 1024, "FFTCOM:%s", jack_get_client_name(g_pclient));
    XStoreName(g_pDisplay, g_window, window_name);

    signal(SIGQUIT, signal_handler );
	signal(SIGTERM, signal_handler );
	signal(SIGHUP, signal_handler );
	signal(SIGINT, signal_handler );

	g_samplerate = jack_get_sample_rate( g_pclient );
	g_buffersize = jack_get_buffer_size( g_pclient );
	
//	int l_buffer_size = FFTSIZE / g_buffersize;
//	if( l_buffer_size == 0 ){
//		l_buffer_size = 1;
//	}
    int l_buffer_size = 2*g_buffersize*sizeof(jack_default_audio_sample_t);
	
	g_prbrx = jack_ringbuffer_create(l_buffer_size);
	g_prbtx = jack_ringbuffer_create(l_buffer_size);

	// register the input port
	g_pportrx = jack_port_register( g_pclient, "rx", JACK_DEFAULT_AUDIO_TYPE, JackPortIsInput, 0);
	
	// register the output port
	g_pporttx = jack_port_register( g_pclient, "tx", JACK_DEFAULT_AUDIO_TYPE, JackPortIsOutput, 0);

    g_pfftcodec = new fftcodec( FFTSIZE, g_samplerate );

	if( !g_pfftcodec ){
		printf("FFTCodec creation failed\n");
		goto clean_up_and_exit;
	}
	
    g_pfftcodec->draw_init( g_pDisplay );
	
	start_thread();
	
	jack_set_process_callback(g_pclient, process, NULL);

	synthesize_tx_data();

	jack_activate(g_pclient);

    if(!l_status) jack_connect(g_pclient, "system:capture_1", "fftcom:rx");

    while(g_process){
		while(XPending(g_pDisplay)){
			XEvent l_event;
			XNextEvent(g_pDisplay,&l_event);
			switch(l_event.type){
            case KeyPress:
                g_process = g_pfftcodec->KeyPressEvaluate((XKeyEvent*)&l_event);
                break;
            case ConfigureNotify:
                ConfigureEvaluate((XConfigureEvent*)&l_event);
                break;
            case EnterNotify:
                g_pfftcodec->EnterNotifyEvaluate(&l_event);
                break;
            case LeaveNotify:
                g_pfftcodec->LeaveNotifyEvaluate(&l_event);
                break;
            case MotionNotify:
                g_pfftcodec->MotionNotifyEvaluate(&l_event);
                break;
            case ClientMessage:
                if( l_event.xclient.message_type == g_wmProtocols &&
                    l_event.xclient.data.l[0] == g_wmDeleteMessage ){
                    g_process = false;
                }
                break;
            default:
                break;
			}
		}
		if(!g_process)
			break;
        g_pfftcodec->draw( g_width, g_height );
		glXSwapBuffers( g_pDisplay, g_window );
	}
	printf("waiting for thread to exit\n");
	pthread_join(g_pthread_process,NULL);
    jack_deactivate(g_pclient);

clean_up_and_exit:
	
	jack_client_close(g_pclient);
	
	glWindowDestroy();
	XCloseDisplay( g_pDisplay );
	return 0;
}
