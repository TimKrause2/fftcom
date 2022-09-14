#include "font.h"
#include <ft2build.h>
#include FT_FREETYPE_H
#include FT_STROKER_H
#include <stdarg.h>
#include <string.h>

FreeTypeFont::FreeTypeFont( void ){
	textures = NULL;
}

FreeTypeFont::~FreeTypeFont( void ){
	Free();
}

void FreeTypeFont::Free( void )
{
	if(textures){
		glDeleteLists( list_base, 128 );
		glDeleteTextures( 128, textures );
		delete [] textures;
		textures = NULL;
	}
}


void make_dlist ( FT_Face ftface, unsigned char charcode, GLuint list_base, GLuint * textures ) {

	FT_UInt glyph_index;
	glyph_index = FT_Get_Char_Index( ftface, charcode );
	FT_Error error = FT_Load_Glyph( ftface, glyph_index, FT_LOAD_RENDER );
	if( error ){
		printf("FT_Load_Glyph:0x%X\n",error);
		return;
	}
// 	printf("make_dlist charcode=0x%hhX\n",charcode);
	if( !(ftface->glyph->format == FT_GLYPH_FORMAT_BITMAP) )
		return;
// 	printf("bitmap_left:%d\n",ftface->glyph->bitmap_left);
// 	printf("bitmap_top:%d\n",ftface->glyph->bitmap_top);
// 	printf("bitmap width:%d\n",ftface->glyph->bitmap.width);
// 	printf("bitmap rows:%d\n",ftface->glyph->bitmap.rows);
// 	printf("bitmap pitch:%d\n",ftface->glyph->bitmap.pitch);
	bool texture_present = (ftface->glyph->bitmap.rows!=0) && (ftface->glyph->bitmap.width!=0);
// 	if( texture_present ){	
// 		unsigned char *data = ftface->glyph->bitmap.buffer;
// 		for(int y=0;y<ftface->glyph->bitmap.rows;y++){
// 			printf("|");
// 			for(int x=0;x<ftface->glyph->bitmap.width;x++){
// 				const char *str;
// 				if(*data==0){
// 					str = " ";
// 				}else if(*data<128){
// 					str = "_";
// 				}else if(*data<192){
// 					str = "*";
// 				}else{
// 					str = "#";
// 				}
// 				printf("%s",str);
// 				data++;
// 			}
// 			printf("|\n");
// 		}
// 	}

	if( texture_present ){
		glBindTexture( GL_TEXTURE_2D, textures[charcode] );
		glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA8, 
									ftface->glyph->bitmap.width,ftface->glyph->bitmap.rows,0,
								GL_ALPHA, GL_UNSIGNED_BYTE,
								(const GLvoid*)ftface->glyph->bitmap.buffer );
		glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
		glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
		glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
		glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
	}

	glNewList(list_base+charcode,GL_COMPILE);
	
	if( texture_present ){
		glBindTexture( GL_TEXTURE_2D, textures[charcode] );
	
		glBegin( GL_QUADS );
			glTexCoord2d( 0.0, 0.0 );
			glVertex2i( ftface->glyph->bitmap_left, -ftface->glyph->bitmap_top);
			glTexCoord2d( 1.0, 0.0 );
			glVertex2i( ftface->glyph->bitmap_left+ftface->glyph->bitmap.width, -ftface->glyph->bitmap_top);
			glTexCoord2d( 1.0, 1.0 );
			glVertex2i( ftface->glyph->bitmap_left+ftface->glyph->bitmap.width, -ftface->glyph->bitmap_top+ftface->glyph->bitmap.rows);
			glTexCoord2d( 0.0, 1.0 );
			glVertex2i( ftface->glyph->bitmap_left, -ftface->glyph->bitmap_top+ftface->glyph->bitmap.rows);
		glEnd();
	}
	
	glTranslated( ftface->glyph->advance.x / 64, 0.0, 0.0 );
	
	glEndList();
}

void FreeTypeFont::Load( const char *filename, unsigned int pixelsize )
{
	if(textures)Free();
	FT_Library ftlibrary;
	FT_Face    ftface;

	FT_Error error = FT_Init_FreeType( &ftlibrary );
	if(error){
		printf("FT_Init_FreeType:0x%X\n",error);
		return;
	}

	error = FT_New_Face( ftlibrary, filename,
											 0, &ftface );
	if( error ){
		printf("FT_New_Face:0x%X\n",error);
		FT_Done_FreeType( ftlibrary );
		return;
	}
	
	if(!FT_IS_SCALABLE( ftface ) ){
		printf("Font is not scalable. Not using it.\n");
		FT_Done_Face( ftface );
		FT_Done_FreeType( ftlibrary );
		return;
	}

	error = FT_Set_Pixel_Sizes( ftface, pixelsize, pixelsize );
	if( error ){
		printf("FT_Set_Pixel_Sizes:0x%X\n",error);
		FT_Done_Face( ftface );
		FT_Done_FreeType( ftlibrary );
		return;
	}

	textures = new GLuint[128];
	
	list_base = glGenLists(128);
	glGenTextures( 128, textures );
	glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );
	for( unsigned char i=0;i<128;i++)
		make_dlist( ftface, i, list_base, textures );
	
	FT_Done_Face( ftface );
	FT_Done_FreeType( ftlibrary );
}

#define TEXT_BUF_LENGTH 1024

void FreeTypeFont::Printf( double x, double y, const char *format, ... )
{
	if(!textures)return;
	char text[TEXT_BUF_LENGTH];
	va_list ap;
	va_start( ap, format );
	vsnprintf( text, TEXT_BUF_LENGTH, format, ap );
	va_end( ap );
	GLsizei length = strlen( text );

	GLint	viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	GLdouble left = viewport[0];
	GLdouble right = viewport[0]+viewport[2];
	GLdouble top = viewport[1];
	GLdouble bottom = viewport[1]+viewport[3];
	GLdouble nearVal = -1.0;
	GLdouble farVal = 1.0;
	glOrtho( left, right, bottom, top, nearVal, farVal );
	
	glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
	glEnable( GL_TEXTURE_2D );
	glEnable( GL_BLEND );
	glDisable( GL_LIGHTING );
	glDisable( GL_DEPTH_TEST );
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	glTranslated( x, y, 0.0 );
	
	glListBase( list_base );
	
	glCallLists( length, GL_UNSIGNED_BYTE, (const GLvoid *)text );
	
	glDisable( GL_TEXTURE_2D );
	glDisable( GL_BLEND );
	
}

// Each time the renderer calls us back we just push another span entry on
// our list.

void
RasterCallback(const int y,
               const int count,
               const FT_Span * const spans,
               void * const user) 
{
  Spans *sptr = (Spans *)user;
  for (int i = 0; i < count; ++i) 
    sptr->push_back(Span(spans[i].x, y, spans[i].len, spans[i].coverage));
}


// Set up the raster parameters and render the outline.

void
RenderSpans(FT_Library &library,
            FT_Outline * const outline,
            Spans *spans) 
{
  FT_Raster_Params params;
  memset(&params, 0, sizeof(params));
  params.flags = FT_RASTER_FLAG_AA | FT_RASTER_FLAG_DIRECT;
  params.gray_spans = RasterCallback;
  params.user = spans;

  FT_Outline_Render(library, outline, &params);
}

void make_dlist_outline( FT_Library library,
												 FT_Face face,
												 unsigned char charcode,
												 GLuint list_base, GLuint * textures,
												 Pixel32 fontCol,
												 Pixel32 outlineCol,
												 double outlineWidth )
{
// 	printf("make_dlist_outline: charcode=%hhd\n",charcode);
	// Load the glyph we are looking for.
	FT_UInt gindex = FT_Get_Char_Index(face, charcode);
	if (FT_Load_Glyph(face, gindex, FT_LOAD_NO_BITMAP) == 0)
	{
		// Need an outline for this to work.
		if (face->glyph->format == FT_GLYPH_FORMAT_OUTLINE)
		{
			// Render the basic glyph to a span list.
			Spans spans;
			RenderSpans(library, &face->glyph->outline, &spans);

			// Next we need the spans for the outline.
			Spans outlineSpans;

			// Set up a stroker.
			FT_Stroker stroker;
			FT_Stroker_New(library, &stroker);
			FT_Stroker_Set(stroker,
											(int)(outlineWidth * 64),
											FT_STROKER_LINECAP_ROUND,
											FT_STROKER_LINEJOIN_ROUND,
											0);

			FT_Glyph glyph;
			if (FT_Get_Glyph(face->glyph, &glyph) == 0)
			{
				FT_Glyph_StrokeBorder(&glyph, stroker, 0, 1);
				// Again, this needs to be an outline to work.
				if (glyph->format == FT_GLYPH_FORMAT_OUTLINE)
				{
					// Render the outline spans to the span list
					FT_Outline *o =
						&reinterpret_cast<FT_OutlineGlyph>(glyph)->outline;
					RenderSpans(library, o, &outlineSpans);
				}

				// Clean up afterwards.
				FT_Stroker_Done(stroker);
				FT_Done_Glyph(glyph);

				// Now we need to put it all together.
				if (!spans.empty())
				{
					// Figure out what the bounding rect is for both the span lists.
					Rect rect(spans.front().x,
										spans.front().y,
										spans.front().x,
										spans.front().y);
					for (Spans::iterator s = spans.begin();
								s != spans.end(); ++s)
					{
						rect.Include(Vec2(s->x, s->y));
						rect.Include(Vec2(s->x + s->width - 1, s->y));
					}
					for (Spans::iterator s = outlineSpans.begin();
								s != outlineSpans.end(); ++s)
					{
						rect.Include(Vec2(s->x, s->y));
						rect.Include(Vec2(s->x + s->width - 1, s->y));
					}

					// This is unused in this test but you would need this to draw
					// more than one glyph.
// 					float bearingX = (float)face->glyph->metrics.horiBearingX/64;
// 					float bearingY = (float)face->glyph->metrics.horiBearingY/64;
// 					float advance = (float)face->glyph->advance.x/64;
// 					printf("bearingX:%f bearingY:%f advance:%f\n",bearingX,bearingY,advance);
// 					printf("rect: xmin:%f, xmax:%f, ymin:%f ymax:%f\n",
// 								 rect.xmin, rect.xmax, rect.ymin, rect.ymax );

					// Get some metrics of our image.
					int imgWidth = rect.Width(),
							imgHeight = rect.Height(),
							imgSize = imgWidth * imgHeight;

					// Allocate data for our image and clear it out to transparent.
					Pixel32 *pxl = new Pixel32[imgSize];
					memset(pxl, 0, sizeof(Pixel32) * imgSize);

					// Loop over the outline spans and just draw them into the
					// image.
					for (Spans::iterator s = outlineSpans.begin();
								s != outlineSpans.end(); ++s)
						for (int w = 0; w < s->width; ++w)
							pxl[(int)((imgHeight - 1 - (s->y - rect.ymin)) * imgWidth
												+ s->x - rect.xmin + w)] =
								Pixel32(outlineCol.r, outlineCol.g, outlineCol.b,
												s->coverage);

					// Then loop over the regular glyph spans and blend them into
					// the image.
					for (Spans::iterator s = spans.begin();
								s != spans.end(); ++s)
						for (int w = 0; w < s->width; ++w)
						{
							Pixel32 &dst =
								pxl[(int)((imgHeight - 1 - (s->y - rect.ymin)) * imgWidth
										+ s->x - rect.xmin + w)];
							Pixel32 src = Pixel32(fontCol.r, fontCol.g, fontCol.b,
																		s->coverage);
							uint8 one_minus_src_alpha = 255 - src.a;
							dst.r = (int)( (dst.r*one_minus_src_alpha + src.r*src.a)/255.0f );
							dst.g = (int)( (dst.g*one_minus_src_alpha + src.g*src.a)/255.0f );
							dst.b = (int)( (dst.b*one_minus_src_alpha + src.b*src.a)/255.0f );
							dst.a = MIN( 255, dst.a + src.a );
// 							dst.r = (int)(dst.r + ((src.r - dst.r) * src.a) / 255.0f);
// 							dst.g = (int)(dst.g + ((src.g - dst.g) * src.a) / 255.0f);
// 							dst.b = (int)(dst.b + ((src.b - dst.b) * src.a) / 255.0f);
// 							dst.a = MIN(255, dst.a + src.a);
						}

					// write the image to a texture
					glBindTexture( GL_TEXTURE_2D, textures[charcode] );
					glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA,
												imgWidth, imgHeight, 0,
												GL_RGBA, GL_UNSIGNED_BYTE, pxl );
					glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
					glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
					glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
					glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
					
					glNewList(list_base+charcode,GL_COMPILE);
					
					glBindTexture( GL_TEXTURE_2D, textures[charcode] );
					
					glBegin( GL_QUADS );
						glTexCoord2d( 0.0, 0.0 );
						glVertex2f( rect.xmin, -rect.ymax-1 );
						glTexCoord2d( 1.0, 0.0 );
						glVertex2f( rect.xmax+1, -rect.ymax-1 );
						glTexCoord2d( 1.0, 1.0 );
						glVertex2f( rect.xmax+1, -rect.ymin );
						glTexCoord2d( 0.0, 1.0 );
						glVertex2f( rect.xmin, -rect.ymin );
					glEnd();
	
					glTranslated( (double)face->glyph->advance.x / 64, 0.0, 0.0 );
					
					glEndList();

					delete [] pxl;
				}else{
					glNewList(list_base+charcode,GL_COMPILE);

					glTranslated( (double)face->glyph->advance.x / 64, 0.0, 0.0 );
					
					glEndList();
				}
				
			}
		}
	}
}

void FreeTypeFont::LoadOutline( const char *filename, unsigned int fontsize,
										Pixel32 fontColor, Pixel32 outlineColor,
									 double outlineWidth)
{
	if(textures)Free();
	FT_Library ftlibrary;
	FT_Face    ftface;

	FT_Error error = FT_Init_FreeType( &ftlibrary );
	if(error){
		printf("FT_Init_FreeType:0x%X\n",error);
		return;
	}

	error = FT_New_Face( ftlibrary, filename,
											 0, &ftface );
	if( error ){
		printf("FT_New_Face:0x%X\n",error);
		FT_Done_FreeType( ftlibrary );
		return;
	}
	
	if(!FT_IS_SCALABLE( ftface ) ){
		printf("Font is not scalable. Not using it.\n");
		FT_Done_Face( ftface );
		FT_Done_FreeType( ftlibrary );
		return;
	}

	error = FT_Set_Pixel_Sizes( ftface, fontsize, fontsize );
	if( error ){
		printf("FT_Set_Pixel_Sizes:0x%X\n",error);
		FT_Done_Face( ftface );
		FT_Done_FreeType( ftlibrary );
		return;
	}

	textures = new GLuint[128];
	
	list_base = glGenLists(128);
	glGenTextures( 128, textures );
	glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );
	for( unsigned char i=0;i<128;i++)
		make_dlist_outline( ftlibrary, ftface, i, list_base, textures, fontColor,
			outlineColor, outlineWidth );
	
	FT_Done_Face( ftface );
	FT_Done_FreeType( ftlibrary );
}



