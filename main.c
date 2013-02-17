/* sierpinsky-hearts.c - Creates a sierpinsky hearts type of thing.
 */

#include<stdlib.h>
#include<stdio.h>
#include<stdint.h>
#include<png.h>
#include<unistd.h>
#include<math.h>
#include<pthread.h>

typedef struct {
  uint8_t r;
  uint8_t g;
  uint8_t b;
  uint8_t a;
} pixel_t;

typedef struct {
  pixel_t **data;
  int w,h;
} image_t;

typedef struct heart_t {
  int x;
  int y;
  double size;
} heart_t;

image_t *drawHeart(image_t *img, heart_t *heart);
image_t *newImage(int w,int h);

image_t *imagePNGLoadFromFile( const char *filename )
{
  FILE *f = fopen( filename, "r" );
  if (!f) {
    printf("File %s does not exist.\n", filename);
    return NULL;
    //exit(0);
  };

  /* Check if it is actually PNG... */
  char sig[8];
  fread( sig, 1, 8, f );
  if (png_sig_cmp( (unsigned char*)sig, 0, 8 )) {
    fclose( f );
    printf("File is not a PNG file.\n");
    return NULL;
  };

  /* Start over from the beginning of the file... */
  fseek( f, 0, SEEK_SET );

  image_t *image = malloc(sizeof(image_t));

  png_structp png_ptr = png_create_read_struct( PNG_LIBPNG_VER_STRING,
						NULL, NULL, NULL );
  if (!png_ptr) {
    fclose( f );
    return NULL;
  };

  png_infop info_ptr = png_create_info_struct( png_ptr );
  if (!info_ptr) {
    png_destroy_read_struct( &png_ptr, (png_infopp)NULL, (png_infopp)NULL );
    fclose( f );
    return NULL;
  };

  png_infop end_info = png_create_info_struct( png_ptr );
  if (!end_info) {
    png_destroy_read_struct( &png_ptr, &info_ptr, (png_infopp)NULL );
    fclose( f );
    return NULL;
  };

  png_init_io( png_ptr, f );

  /* Lets go the low level route, shall we? */

  png_read_info( png_ptr, info_ptr );

  png_uint_32 w, h;
  int depth, color_type, filter, compression, interlace;
  png_get_IHDR( png_ptr, info_ptr, &w, &h, &depth, &color_type, &interlace, &compression, &filter );

  /* Setup the input transformations... */

  if (color_type == PNG_COLOR_TYPE_PALETTE) {
    png_set_palette_to_rgb( png_ptr );
  };
  if ((color_type == PNG_COLOR_TYPE_GRAY) && (depth < 8)) {
#if PNG_LIBPNG_VER_MINOR >= 4
    png_set_expand_gray_1_2_4_to_8( png_ptr );
#else
    png_set_gray_1_2_4_to_8( png_ptr );
#endif
  };
  if (depth == 16) {
    png_set_strip_16( png_ptr );
  };
  if (depth < 8) {
    png_set_packing( png_ptr );
  };
  if (color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
    png_set_gray_to_rgb( png_ptr );
    if (color_type == PNG_COLOR_TYPE_GRAY)
      png_set_filler( png_ptr, 0xffff, PNG_FILLER_AFTER );
  };
  if (color_type == PNG_COLOR_TYPE_RGB) {
    png_set_filler( png_ptr, 0xffff, PNG_FILLER_AFTER );
  };
  if (png_get_valid( png_ptr, info_ptr, PNG_INFO_tRNS)) {
    png_set_tRNS_to_alpha( png_ptr );
  };

  /* Read out the data... */
  unsigned char *row_pointers[h];
  int i, j;
  for (i=0; i<h; i++) {
    row_pointers[i] = malloc(w*4);
  };

  png_read_image( png_ptr, row_pointers );

  /* Yay! Data! */
  image->data = malloc(sizeof(pixel_t*)*w);//malloc(w * h * 4);
  for (i=0; i<w; i++) {
    image->data[i] = malloc(sizeof(pixel_t)*h);
    for (j=0; j<h; j++) {
      image->data[i][j].r = row_pointers[h-j-1][i*4];
      image->data[i][j].g = row_pointers[h-j-1][i*4+1];
      image->data[i][j].b = row_pointers[h-j-1][i*4+2];
      image->data[i][j].a = row_pointers[h-j-1][i*4+3];
    }
  }

  image->w = w;
  image->h = h;

  fclose( f );

  return image;
  return NULL;
}

void imagePNGWriteStatusCallback( png_structp png_ptr, png_uint_32 row,
				  int pass )
{
}

int imagePNGWriteToFile( const char *filename, image_t *image )
{
  /* Open the file... */
  FILE *f = fopen( filename, "wb" );
  if (!f) {
    return 1;
  };

  /* Start setting up libpng to write... */
  png_structp png_ptr = png_create_write_struct( PNG_LIBPNG_VER_STRING,
						 NULL, NULL, NULL );
  if (!png_ptr) {
    /* Error? */
  };

  png_infop info_ptr = png_create_info_struct( png_ptr );
  if (!info_ptr) {
    png_destroy_write_struct( &png_ptr, NULL );
    /* Error? */
  };

  /* Init the io... */
  png_init_io( png_ptr, f );

  /* Tell it the write status function we want to use... */
  png_set_write_status_fn( png_ptr, imagePNGWriteStatusCallback );

  /* Setup the header... */
  png_set_IHDR( png_ptr, info_ptr, image->w, image->h, 8,
		PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT );

  /* Setup the pointers to the rows in the image... */
  unsigned char *row_pointers[image->h];
  int i, j=0;
  for (i=0; i<image->h; i++) {
    row_pointers[image->h-i-1] = malloc(image->w*4);
    for (j=0; j<image->w; j++) {
      int y = i;//image->h - i - 1;
      row_pointers[image->h-i-1][j*4+0] = image->data[j][y].r;
      row_pointers[image->h-i-1][j*4+1] = image->data[j][y].g;
      row_pointers[image->h-i-1][j*4+2] = image->data[j][y].b;
      row_pointers[image->h-i-1][j*4+3] = image->data[j][y].a;
    }
    //row_pointers[image->h-i-1] = &image->data[j];
    //j+=4*image->w;
  }

  png_set_rows( png_ptr, info_ptr, row_pointers );

  /* Write out the image... */
  png_write_png( png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL );

  /* Done! Let's clean up. */
  fclose( f );

  return 0;
}

void imageAddBorder(image_t *img)
{
  int x;
  for (x=0; x<img->w; x++) {
    img->data[x][0].r = 0;
    img->data[x][img->h-1].r = 0;
  }
  int y;
  for (y=0; y<img->h; y++) {
    img->data[0][y].r = 0;
    img->data[img->w-1][y].r = 0;
  }
}

#define HEART_FN_CACHE_SIZE 10000

double heart_Function_Cache[HEART_FN_CACHE_SIZE+1];

// Gives the r in the direction of dx,dy
double heartFunction(double dx, double dy, double size)
{
  double tmp = dx;
  dx = dy;
  dy = tmp;
  double theta = atan2(dy,dx);
  //theta = M_PI/2.0 - theta;
  //double r = (2.0 - 2.0 * sin(theta) + sin(theta) * sqrt(fabs(cos(theta))) / (sin(theta) + 1.4)) / 4.0 * size;
  //printf("%f => %f\n", theta, r);
  int index = ((double)HEART_FN_CACHE_SIZE*(theta+M_PI)/(M_PI*2));
  //printf("%f => %f => %i\n", theta, (theta+M_PI)/(M_PI*2), index);
  double r = heart_Function_Cache[index] * size;
  //double r = fabs(theta) / 3.14159265 * size;
  return r;
}

void buildHeartFunctionCache(void)
{
  int i;
  for (i=0; i<HEART_FN_CACHE_SIZE; i++) {
    double theta = (M_PI*2) * ((double)i / (double)HEART_FN_CACHE_SIZE) - M_PI;
#if 1
    theta = M_PI/2.0 - theta;
    double r = (2.0 - 2.0 * sin(theta) + sin(theta) * sqrt(fabs(cos(theta))) / (sin(theta) + 1.4)) / 4.0;
#else
    double r = fabs(theta) / 3.14159265;
#endif
    heart_Function_Cache[i] = r;
  }
  heart_Function_Cache[i] = heart_Function_Cache[0];
}

// Returns if the point is inside the heart function
int heartFn(double dx, double dy, double size)
{
  double r = heartFunction(dx, dy, size);
  double d = sqrt(dx*dx+dy*dy);
  if (dy > 0) {
    //printf("Hello? %f %f %f %f %f\n", dx, dy, theta, r, d);
  }
  //return (tmp>0)?1:0;
  return (d<r)?1:0;
}

void findHeartCentroid(int size, double *c_x, double *c_y)
{
  // First render the heart
  image_t *img = newImage(size+2,size*2);
  heart_t h;
  h.x = size/2;
  h.y = size;
  h.size = size;
  drawHeart(img, &h);

  // Now, find the centroid
  *c_x = *c_y = 0.0;
  int x,y;
  int npoints = 0;
  for (x=0; x<size+2; x++) {
    for (y=0; y<size*2; y++) {
      if (img->data[x][y].r > 127) {
	*c_x += x;
	*c_y += y;
	npoints++;
      }
    }
    free(img->data[x]);
  }
  free(img->data);
  free(img);

  *c_x = *c_x / (double)npoints;
  *c_y = *c_y / (double)npoints;
}

int doHeartsOverlap(heart_t *a, heart_t *b)
{
  double dx = b->x-a->x;
  double dy = b->y-a->y;
  double d = sqrt(dx*dx+dy*dy);
  if (d<a->size && d<b->size) {
    return 1;
  }
  if (d > a->size+b->size) {
    return 0;
  }
  double d1 = heartFunction(b->x-a->x, b->y-a->y, a->size);
  double d2 = heartFunction(a->x-b->x, a->y-b->y, b->size);
  //printf("%f %f %f\n", d1,d2,d);
  //return (d1+d2>d)?1:0;
  if (d1+d2>d) return 1;

  // Alright, now we brute force it.
  // Render them both on the same canvas.
  // See if they intersect.
  int x,y;
  int minx = a->x - a->size;
  int miny = a->y - a->size;
  int maxx = a->x + a->size;
  int maxy = a->y + a->size;
  if (b->x-b->size > minx) minx = b->x-b->size;
  if (b->x+b->size > maxx) maxx = b->x+b->size;
  if (b->y-b->size > miny) miny = b->y-b->size;
  if (b->y+b->size > maxy) maxy = b->y+b->size;

  for (x=minx; x<maxx; x++) {
    for (y=miny; y<maxy; y++) {
      if (heartFn(x-a->x, y-a->y, a->size) && heartFn(x-b->x, y-b->y, b->size)) {
	return 1;
      }
    }
  }
  return 0;
}

double dist(double x1, double y1, double x2, double y2)
{
  double dx = x2-x1;
  double dy = y2-y1;
  return sqrt(dx*dx + dy*dy);
}

int willHeartFit(image_t *img, int x, int y, int size)
{
  int dx,dy;
  int min_dx = (x-size<0)?-x:-size;
  int max_dx = (x+size>=img->w)?img->w-x-1:size;
  int min_dy = (y-size<0)?-y:-size;
  int max_dy = (y+size>=img->h)?img->h-y-1:size;
  //printf("%i,%i,%i => %i %i %i %i\n", x,y,size,min_dx, max_dx, min_dy, max_dy);
  for (dx=min_dx; dx<max_dx; dx++) {
    for (dy=min_dy; dy<max_dy; dy++) {
      // Are we within the <size> circle?
      //double d = dist(0.0,0.0, (double)dx, (double)dy);
      //if (d < size) {
      if (img->data[x+dx][y+dy].g < 127) {
	if (heartFn(dx,dy, size+3)) {
	  return 0;
	}
      }
    }
  }
  return 1;
}

// This returns the largest heart that will fit at x,y
int findDistanceToPixel(image_t *img, int x, int y)
{
  // Walk up a counter of sizes
  int size = 1;
  int incr = 10;
  while (1) {
    if (!willHeartFit(img, x,y, size)) {
      if (incr == 1) {
	break;
      }
      size -= incr;
      if (size < 1) size = 1;
      incr = incr / 2;
      //printf("%i %i\n", size, incr);
    }
    size += incr;
  }
  return size;
}

typedef struct distToPixelThread_t {
  image_t *img;
  image_t *dst;
  int min_x;
  int max_x;
} distToPixelThread_t;

void *findDistanceToPixel_thread(void *arg)
{
  distToPixelThread_t *data = arg;

  int x,y;
  for (x=data->min_x; x<data->max_x; x++) {
    for (y=0; y<data->img->h; y++) {
      if (data->img->data[x][y].r > 127) {
	data->dst->data[x][y].r = findDistanceToPixel(data->img, x, y);
	data->dst->data[x][y].a = 255;
      }
    }
    //printf("x=%i on %i-%i\n", x, data->min_x, data->max_x);
  }

  return NULL;
}

image_t *newImage(int w,int h)
{
  image_t *dst = malloc(sizeof(image_t));
  dst->w = w;
  dst->h = h;
  dst->data = malloc(sizeof(pixel_t*)*w);//malloc(w * h * 4);
  int i;
  for (i=0; i<w; i++) {
    dst->data[i] = malloc(sizeof(pixel_t)*h);
  }
  return dst;
}

// Find the distances of all white areas to black areas.
image_t *findDistances(image_t *src_img)
{
  int i;
  image_t *dst = malloc(sizeof(image_t));
  int w = src_img->w;
  int h = src_img->h;
  dst->w = w;
  dst->h = h;
  dst->data = malloc(sizeof(pixel_t*)*w);//malloc(w * h * 4);
  for (i=0; i<w; i++) {
    dst->data[i] = malloc(sizeof(pixel_t)*h);
  }

  // Create a few threads
  int num_threads = 16;
  pthread_t threads[num_threads];
  distToPixelThread_t t_args[num_threads];
  for (i=0; i<num_threads; i++) {
    t_args[i].img = src_img;
    t_args[i].dst = dst;
    t_args[i].min_x = i*(w/num_threads);
    t_args[i].max_x = (i+1)*(w/num_threads);
    pthread_create(&threads[i], NULL, findDistanceToPixel_thread, (void*)&t_args[i]);
  }

  for (i=0; i<num_threads; i++) {
    pthread_join(threads[i], NULL);
  }

  // For each pixel...
#if 0
  int x;
  for (x=0; x<w; x++) {
    int y;
    for (y=0; y<h; y++) {
      if (src_img->data[x][y].r > 127) {
	dst->data[x][y].r = findDistanceToPixel(src_img, x,y);
	dst->data[x][y].a = 255;
      }
      //printf("y=%i\n", y);
    }
    printf("x=%i\n", x);
  }
#endif
  return dst;
}

image_t *drawHeart(image_t *img, heart_t *heart)
{
  int x = heart->x;
  int y = heart->y;
  double size = heart->size;
  int min_dx = (x-size<0)?-x:-size;
  int max_dx = (x+size>=img->w)?img->w-1:size;
  int min_dy = (y-size<0)?-y:-size;
  int max_dy = (y+size>=img->h)?img->h-1:size;
  //printf("%i %i %f => %i %i %i %i\n", x, y, size, min_dx, max_dx, min_dy, max_dy);
  int dx,dy;
  for (dx=-x/*min_dx*/; dx<img->w-x-1/*max_dx*2*/; dx++) {
    for (dy=-y/*min_dy*/; dy<img->h-y-1/*max_dy*/; dy++) {
      if (dx < 0) {
	//printf("Hello?\n");
      }
      if (heartFn(dx,dy, size)) {
	//printf("Printed!\n");
	img->data[x+dx][y+dy].r = 255;
	img->data[x+dx][y+dy].g = 0;
	img->data[x+dx][y+dy].b = 0;
	img->data[x+dx][y+dy].a = 255;
      } else {
	//img->data[x+dx][y+dy].r = 0;
      }
      //img->data[x+dx][y+dy].g = 255;//dy*128;
      //img->data[x+dx][y+dy].b = 0;
      //img->data[x+dx][y+dy].a = 255;
    }
  }
  return img;
}

// Determines if the given point is on a "ridge" in the image.
int imageIsOnRidge(image_t *img, int x, int y)
{
  // Check the five-neighborhood. If they're all lower, then we're on a ridge/peak.
  int dx,dy;
  int size = 3;
  if (x<size || img->w-size<x || y<size || img->h-size<y) {
    return 0;
  }
  int min_dx = -size;
  int max_dx = size;
  int min_dy = -size;
  int max_dy = size;
  int comp = img->data[x][y].r;
  //printf("%i,%i => %i,%i,%i,%i\n", x,y, min_dx,max_dx, min_dy, max_dy);
  for (dx=min_dx; dx<max_dx; dx++) {
    for (dy=min_dy; dy<max_dy; dy++) {
      if (comp+3 < img->data[x+dx][y+dy].r && img->data[x+dx][y+dy].a > 127) {
	//printf("%i,%i %i %i\n", x,y, comp, img->data[x+dx][y+dy].r);
	return 0;
      }
    }
  }
  return 1;
}

int main(int argc, char **argv)
{
  buildHeartFunctionCache();
  printf("Built heart function cache.\n");

  // Render a test heart
  image_t *test_heart_img = newImage(4096,4096);
  heart_t test_heart;
  test_heart.x = 2048;
  test_heart.y = 2048;
  test_heart.size = 2040;
  drawHeart(test_heart_img, &test_heart);
  imagePNGWriteToFile("heart.png", test_heart_img);
  //return 0;

  time_t start_time = time(NULL);
  int depth = 1;

  // Load the image to sierpinskify
  //image_t *src_img = imagePNGLoadFromFile("input2.png");
  image_t *src_img = imagePNGLoadFromFile(argv[1]);

  // Invert the colors
  int x,y;
  for (x=0; x<src_img->w; x++) {
    for (y=0; y<src_img->h; y++) {
      src_img->data[x][y].r = 255-src_img->data[x][y].r;
      src_img->data[x][y].g = 255-src_img->data[x][y].g;
      src_img->data[x][y].b = 255-src_img->data[x][y].b;
      src_img->data[x][y].a = 255;
    }
  }

  imageAddBorder(src_img);

  // Create the final destination image
  image_t *final_img = newImage(src_img->w,src_img->h);
  for (x=0; x<src_img->w; x++) {
    for (y=0; y<src_img->h; y++) {
      final_img->data[x][y] = src_img->data[x][y];
    }
  }

  // Figure out the distances

 repeat:;
  image_t *distances = findDistances(final_img); // Note that this is now updated when we add hearts to the image, so it doesn't have to iterate.
  char buf[1024];
  snprintf(buf, 1024, "distances_%i.png", depth);
  imagePNGWriteToFile(buf, distances);

  // Turn peaks in the distance data into hearts
  printf("%ld: Converting image peaks.\n", time(NULL)-start_time);
  heart_t *hearts = malloc(sizeof(heart_t)*1024);
  int nhearts = 0;
  for (x=1; x<distances->w-1; x++) {
    for (y=1; y<distances->h-1; y++) {
      int comp = distances->data[x][y].r;
      if (distances->data[x+1][y].r < comp &&
	  distances->data[x-1][y].r < comp &&
	  distances->data[x][y+1].r < comp &&
	  distances->data[x][y-1].r < comp &&
	  distances->data[x+1][y+1].r < comp &&
	  distances->data[x-1][y-1].r < comp &&
	  distances->data[x-1][y+1].r < comp &&
	  distances->data[x+1][y-1].r < comp) {
	//drawHeart(final_img, x, y, comp);
	if (comp > 15) {
	  hearts[nhearts].x = x;
	  hearts[nhearts].y = y;
	  hearts[nhearts].size = comp;
	  nhearts++;
	}
      }
    }
  }

  // Detect ridges, and turn them into hearts
  printf("%ld: Converting image ridges.\n", time(NULL)-start_time);
  for (x=1; x<distances->w-1; x++) {
    for (y=1; y<distances->h-1; y++) {
      if (imageIsOnRidge(distances, x, y)) {
	if (distances->data[x][y].r > 5) {
	  hearts = realloc(hearts, sizeof(heart_t)*(nhearts+1));
	  hearts[nhearts].x = x;
	  hearts[nhearts].y = y;
	  hearts[nhearts].size = distances->data[x][y].r;
	  nhearts++;
	}
      }
    }
  }

  // Disqualify hearts if they're too close together.
  printf("%ld: We have %i potential heart candidates.\n", time(NULL)-start_time, nhearts);
  char *enabled_hearts = malloc(sizeof(char)*nhearts);
  memset(enabled_hearts, 1, nhearts);
  int valid_hearts = nhearts;
  int i, j;
  for (i=0; i<nhearts; i++) {
    for (j=i+1; j<nhearts; j++) {
      if (enabled_hearts[i] && enabled_hearts[j]) {
	if (doHeartsOverlap(&hearts[i], &hearts[j])) {
	  //printf("Disabling %i (%f) or %i (%f)\n", i, hearts[i].size, j, hearts[j].size);
	  if (hearts[i].size > hearts[j].size) {
	    enabled_hearts[j] = 0;
	    valid_hearts--;
	  } else {
	    enabled_hearts[i] = 0;
	    valid_hearts--;
	    //break;
	  }
	}
      }
    }
    if (i%100==0)
      printf("Working on heart %i (There are still %i valid hearts)\n", i, valid_hearts);
  }

  // Render the hearts!
  printf("%ld: Rendering %i hearts.\n", time(NULL)-start_time, valid_hearts);
  for (i=0; i<nhearts; i++) {
    if (enabled_hearts[i]) {
      drawHeart(final_img, &hearts[i]);
      drawHeart(src_img, &hearts[i]);
    }
  }

#if 0
  // Let's add the hearts to the distance map, so we don't have to regen the
  // entire thing from scratch next time.
  printf("%ld: Rendering hearts to distance map.\n",time(NULL)-start_time);
  for (i=0; i<nhearts; i++) {
    if (!enabled_hearts[i]) continue;
    heart_t *h = &hearts[i];

    // Render the clear alpha channel through
    for (x=h->x-h->size; x<h->x+h->size; x++) {
      for (y=h->y-h->size; y<h->y+h->size; y++) {
	if (x>0 && y>0 && x<distances->w && y<distances->h) {
	  if (heartFn(x-h->size, y-h->size, h->size)) {
	    distances->data[x][y].a = 0;
	  }
	}
      }
    }

    // Now, expand the heart, filling in new distance values as we go.
    int size = h->size+1;
    while (1) {
      int did_reset_val = 0;

      for (x=h->x-size; x<h->x+size; x++) {
	for (y=h->y-size; y<h->y+size; y++) {
	  if (x>0 && y>0 && x<distances->w && y<distances->h) {
	    if (distances->data[x][y].a > 127) {
	      if (distances->data[x][y].r > size) {
		if (heartFn(h->x, h->y, size)) {
		  did_reset_val = 1;
		  distances->data[x][y].r = size;
		}
	      }
	    }
	  }
	}
      }

      if (!did_reset_val) {
	printf("Stopped resetting values at size=%i for heart %i\n", size, i);
	break;
      }
      size++;
    }

#if 0
    for (x=0; x<distances->w; x++) {
      for (y=0; y<distances->h; y++) {
	if (distances->data[x][y].a < 127) continue;
	double dx = x - hearts[i].x;
	double dy = y - hearts[i].y;
	double d_to_center = sqrt(dx*dx+dy*dy);
	int needs_long_processing = 1;
	if (d_to_center < hearts[i].size) {
	  if (heartFn(x-hearts[i].x, y-hearts[i].y, hearts[i].size)) {
	    distances->data[x][y].a = 0;
	    needs_long_processing = 0;
	  }
	}
	if (needs_long_processing) {
	  // If the largest heart that can fit in this pixel is already smaller
	  // than the distance to heart i - size of heart i, then we can save
	  // the processing required to determine the largest heart.
	  if (distances->data[x][y].r > d_to_center-hearts[i].size) {
	    // Figure out the largest heart we can fit, as limited by heart i.
	    heart_t trial;
	    trial.x = x;
	    trial.y = y;
	    trial.size = 1;
	    while (trial.size<distances->data[x][y].r) {
	      if (doHeartsOverlap(&hearts[i], &trial)) {
		trial.size -= 1;
		distances->data[x][y].r = trial.size;
		break;
	      }
	      trial.size += 1;
	    }
	  }
	}
      }
    }
#endif
    printf("%ld: Now on heart %i\n", time(NULL)-start_time, i);
  }
#endif

  // Iterate! Sierpinsky is iterative!
  depth++;
  printf("Acheived depth %i in %ld seconds.\n", depth, time(NULL)-start_time);
  snprintf(buf, 1024, "intermediate_%i.png", depth);
  imagePNGWriteToFile(buf, final_img);
  if (depth < 5) {
    goto repeat;
  }

  imagePNGWriteToFile("final.png", final_img);

  return 0;
}
