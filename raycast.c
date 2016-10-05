#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


typedef struct {
  char *type; // 0 = cylinder, 1 = sphere, 2 = teapot
  double color[3];
  union {
    struct {
      double width;
      double height;
    } camera;
    struct {
      double position[3];
      double radius;
    } sphere;
    struct {
      double position[3];
      double normal[3];
    } plane;
  };
} Object;

typedef struct {
  double origin[3];
  double dir[3];
}ray;



typedef struct ObjectInfo{
  Object *objectArray;
  size_t objectNumber;
} ObjectInfo;

typedef struct ImageData{
  double width;
  double height;
  char* color;
}ImageData;

struct ObjectInfo read_scene(FILE* json);

int line = 1;

// next_c() wraps the getc() function and provides error checking and line
// number maintenance
int next_c(FILE* json) {
  int c = fgetc(json);
#ifdef DEBUG
  printf("next_c: '%c'\n", c);
#endif
  if (c == '\n') {
    line += 1;
  }
  if (c == EOF) {
    fprintf(stderr, "Error: Unexpected end of file on line number %d.\n", line);
    exit(1);
  }
  return c;
}


// expect_c() checks that the next character is d.  If it is not it emits
// an error.
void expect_c(FILE* json, int d) {
  int c = next_c(json);
  if (c == d) return;
  fprintf(stderr, "Error: Expected '%c' on line %d.\n", d, line);
  exit(1);
}


// skip_ws() skips white space in the file.
void skip_ws(FILE* json) {
  int c = next_c(json);
  while (isspace(c)) {
    c = next_c(json);
  }
  ungetc(c, json);
}


// next_string() gets the next string from the file handle and emits an error
// if a string can not be obtained.
char* next_string(FILE* json) {
  char buffer[129];
  int c = next_c(json);
  if (c != '"') {
    fprintf(stderr, "Error: Expected string on line %d.\n", line);
    exit(1);
  }
  c = next_c(json);
  int i = 0;
  while (c != '"') {
    if (i >= 128) {
      fprintf(stderr, "Error: Strings longer than 128 characters in length are not supported.\n");
      exit(1);
    }
    if (c == '\\') {
      fprintf(stderr, "Error: Strings with escape codes are not supported.\n");
      exit(1);
    }
    if (c < 32 || c > 126) {
      fprintf(stderr, "Error: Strings may contain only ascii characters.\n");
      exit(1);
    }
    buffer[i] = c;
    i += 1;
    c = next_c(json);
  }
  buffer[i] = 0;
  return strdup(buffer);
}

double next_number(FILE* json) {
  double value;
  fscanf(json, "%lf", &value);
  // Error check this..
  return value;
}

double* next_vector(FILE* json) {
  double* v = malloc(3*sizeof(double));
  expect_c(json, '[');
  skip_ws(json);
  v[0] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[1] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[2] = next_number(json);
  skip_ws(json);
  expect_c(json, ']');
  return v;
}


struct ObjectInfo read_scene(FILE* json) {
  int c;
  ObjectInfo object;
  object.objectArray = NULL;
  object.objectNumber = 0;

  skip_ws(json);

  // Find the beginning of the list
  expect_c(json, '[');

  skip_ws(json);

  // Find the objects

  while (1) {
    c = fgetc(json);
    if (c == ']') {
      fprintf(stderr, "Error: This is the worst scene file EVER.\n");
      fclose(json);
      return object;
    }
    if (c == '{') {
      skip_ws(json);

      // Parse the object
      char* key = next_string(json);
      if (strcmp(key, "type") != 0) {
	       fprintf(stderr, "Error: Expected \"type\" key on line number %d.\n", line);
	        exit(1);
      }

      skip_ws(json);

      expect_c(json, ':');

      skip_ws(json);

      char* value = next_string(json);
      object.objectNumber += 1;
      object.objectArray =realloc(object.objectArray, sizeof(Object)*object.objectNumber);
      object.objectArray[object.objectNumber-1].type = value;

      if (strcmp(value, "camera") == 0) {
      } else if (strcmp(value, "sphere") == 0) {
      } else if (strcmp(value, "plane") == 0) {
      } else {
	       fprintf(stderr, "Error: Unknown type, \"%s\", on line number %d.\n", value, line);
	        exit(1);
      }

      skip_ws(json);

      while (1) {
        c = next_c(json);
        if (c == '}') {
	  // stop parsing this object
          break;
  	     }
         else if (c == ',') {
  	  // read another field
          skip_ws(json);
  	      char* key = next_string(json);
  	      skip_ws(json);
  	      expect_c(json, ':');
  	      skip_ws(json);
  	  if ((strcmp(key, "width") == 0) ||
  	      (strcmp(key, "height") == 0) ||
  	      (strcmp(key, "radius") == 0)) {
  	    double value = next_number(json);
        if (value <= 0){
          fprintf(stderr, "Error: Width, Height, Radius must be greater than 0");
          exit(1);
        }
        else {
          if (strcmp(key, "width") == 0){
            object.objectArray[object.objectNumber-1].camera.width = value;
          }
          else if (strcmp(key, "height") == 0){
              object.objectArray[object.objectNumber-1].camera.height = value;
          }
 		     else if (strcmp(key, "radius") == 0){
              object.objectArray[object.objectNumber-1].sphere.radius = value;
          }
        }
  	   }

       else if ((strcmp(key, "color")==0)||(strcmp(key, "position")==0)||(strcmp(key, "normal")==0)){
          double* value = next_vector(json);

          if((strcmp(key, "color")==0)){
            object.objectArray[object.objectNumber-1].color[0]=value[0];
            object.objectArray[object.objectNumber-1].color[1]=value[1];
            object.objectArray[object.objectNumber-1].color[2]=value[2];
          }
          else if((strcmp(key, "position")==0)){
            object.objectArray[object.objectNumber-1].sphere.position[0]=value[0];
            object.objectArray[object.objectNumber-1].sphere.position[1]=value[1];
            object.objectArray[object.objectNumber-1].sphere.position[2]=value[2];
          }
          else if((strcmp(key, "normal")==0)){
            object.objectArray[object.objectNumber-1].plane.normal[0]=value[0];
            object.objectArray[object.objectNumber-1].plane.normal[1]=value[1];
            object.objectArray[object.objectNumber-1].plane.normal[2]=value[2];
          }
       }
      else {
  	    fprintf(stderr, "Error: Unknown property, \"%s\", on line %d.\n",
  		    key, line);
  	    //char* value = next_string(json);
  	  }
  	  skip_ws(json);
  	} else {
  	  fprintf(stderr, "Error: Unexpected value on line %d\n", line);
  	  exit(1);
  	}
        }
        skip_ws(json);
        c = next_c(json);
        if (c == ',') {
  	// noop
  	skip_ws(json);
        } else if (c == ']') {
  	fclose(json);
  	return object;
        } else {
  	fprintf(stderr, "Error: Expecting ',' or ']' on line %d.\n", line);
  	exit(1);
      }
    }
  }
}

static inline double sqr(double v) {
  return v*v;
}



static void normalize(double* v) {
  double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

void shade_pixel(double *color, int row, int col, ImageData *image){
  image->color[(int)(row*image->width*3 + col*3)] = (char)(color[0]*255);
  image->color[(int)(row*image->width*3 + col*3+1)] = (char)(color[1]*255);
  image->color[(int)(row*image->width*3 + col*3+2)] = (char)(color[2]*255);
}

double planeintersection(double* Ro, double* Rd, double* position, double* normal){
  double val = (normal[0]* Rd[0])+(normal[1]* Rd[1])+(normal[2]*Rd[2]);
  if(fabs(val)< .0001){
    return -1;
  }
  double secval[3];
  int i;
  for (i=0; i<=2; i++){
    secval[i] = position[i]-Ro[i];
  }

  double thirdval = (secval[0]* normal[0])+(secval[1]* normal[1])+(secval[2]* normal[2]);
  double final = thirdval/val;

  if (final < 0){
    return -1;
  }
  return final;
}


void ppmprint(ImageData *image, FILE* outputfile, int ppmmagicnumber){
  size_t number;
  
  int imagesize = image->width*image->height*3;
  if (ppmmagicnumber == 6){
    fprintf(outputfile, "P%d\n%lf %lf\n%d\n", ppmmagicnumber, image->width, image->height, 255);
    int i;
    for(i=0; i<imagesize; i++){
      char c = image->color[i];
      if(i%3 !=0){
        fwrite(&c, 1, 1, outputfile);
      }
    }
  }
  fclose(outputfile);
}

double sphereintersection(double* Ro, double* Rd, double* C, double radius){
  double val = (sqr(Rd[0]) + sqr(Rd[1]) + sqr(Rd[2]));
  double secval = 2 * (Rd[0]*(Ro[0]-C[0]) + Rd[1]*(Ro[1]-C[1]) + Rd[2]*(Ro[2]-C[2]));
  double thirdval = sqr(Ro[0]-C[0]) + sqr(Ro[1]-C[1]) + sqr(Ro[2]-C[2]) - sqr(radius);

  double det = sqr(secval) - 4 * val * thirdval;

  if (det<0){
    return -1;
  }
  det = sqrt(det);

  double t1 = (-secval - det)/(2*val);
  if (t1 > 0){
    return t1;
  }
  double t2 = (-secval + det)/(2*val);
  if (t1 >0){
    return t2;
  }
  return -1;
}


/*
double cylinder_intersection(double* Ro, double* Rd,
			     double* C, double r) {
  // Step 1. Find the equation for the object you are
  // interested in..  (e.g., cylinder)
  //
  // x^2 + z^2 = r^2
  //
  // Step 2. Parameterize the equation with a center point
  // if needed
  //
  // (x-Cx)^2 + (z-Cz)^2 = r^2
  //
  // Step 3. Substitute the eq for a ray into our object
  // equation.
  //
  // (Rox + t*Rdx - Cx)^2 + (Roz + t*Rdz - Cz)^2 - r^2 = 0
  //
  // Step 4. Solve for t.
  //
  // Step 4a. Rewrite the equation (flatten).
  //
  // -r^2 +
  // t^2 * Rdx^2 +
  // t^2 * Rdz^2 +
  // 2*t * Rox * Rdx -
  // 2*t * Rdx * Cx +
  // 2*t * Roz * Rdz -
  // 2*t * Rdz * Cz +
  // Rox^2 -
  // 2*Rox*Cx +
  // Cx^2 +
  // Roz^2 -
  // 2*Roz*Cz +
  // Cz^2 = 0
  //
  // Steb 4b. Rewrite the equation in terms of t.
  //
  // t^2 * (Rdx^2 + Rdz^2) +
  // t * (2 * (Rox * Rdx - Rdx * Cx + Roz * Rdz - Rdz * Cz)) +
  // Rox^2 - 2*Rox*Cx + Cx^2 + Roz^2 - 2*Roz*Cz + Cz^2 - r^2 = 0
  //
  // Use the quadratic equation to solve for t..
  double a = (sqr(Rd[0]) + sqr(Rd[2]));
  double b = (2 * (Ro[0] * Rd[0] - Rd[0] * C[0] + Ro[2] * Rd[2] - Rd[2] * C[2]));
  double c = sqr(Ro[0]) - 2*Ro[0]*C[0] + sqr(C[0]) + sqr(Ro[2]) - 2*Ro[2]*C[2] + sqr(C[2]) - sqr(r);

  double det = sqr(b) - 4 * a * c;
  if (det < 0) return -1;

  det = sqrt(det);

  double t0 = (-b - det) / (2*a);
  if (t0 > 0) return t0;

  double t1 = (-b + det) / (2*a);
  if (t1 > 0) return t1;

  return -1;
}ObjectInfo
*/

int main(int argc, char *argv[]) {
  if (argc != 5){
    fprintf(stderr, "Usage: raycast width height input.json output.ppm");
    return(1);
  }

  if (atoi(argv[1]) <= 0){
    fprintf(stderr, "error: width must be greater than 0");
  }
  if (atoi(argv[2]) <= 0){
    fprintf(stderr, "error: height must be greater than 0");
  }

  FILE* json = fopen(argv[3], "rb");

  if (json == NULL) {
    fprintf(stderr, "Error: Could not open file \"%s\"\n", argv[3]);
    exit(1);
  }

  //create object?

  ObjectInfo object;

  int objectsnumber;

  double cameraheight;
  double camerawidth;



  double radius;



  int M = atoi(argv[1]);
  int N = atoi(argv[2]);


  object = read_scene(json);

  objectsnumber = object.objectNumber;

  int i;

  int ppmmagicnumber=6;

  for (i = 0; i < objectsnumber; i++){
    if (object.objectArray[i].camera.width && object.objectArray[1].camera.height) {
      camerawidth = object.objectArray[i].camera.width;
      cameraheight = object.objectArray[i].camera.height;
    }
  }

  if(!camerawidth || !cameraheight){
    fprintf(stderr, "Error: invalid camera width/height");
    exit(1);
  }

  double pixelheight = cameraheight / M;
  double pixelwidth = camerawidth / N;

  ImageData *image = (ImageData *)malloc(sizeof(ImageData));
  image->height = N;
  image->width = M;
  image->color = (char *)malloc(sizeof(char) * image->height * image->width * 4);

  int y;
  int x;
  for (y=0; y < M; y += 1) {
    for (x=0; x < N; x += 1) {
      double Ro[3] = {0, 0, 0};
      // Rd = normalize(P - Ro)
      double Rd[3] = {object.objectArray[i].sphere.position[0] -
        (camerawidth/2) + pixelwidth * (x + 0.5),
        object.objectArray[i].sphere.position[1] -
        (cameraheight/2) + pixelheight * (y+0.5), 1};
      normalize(Rd);

      int helper = 0;
      double best_t = INFINITY;
      for (i=0; i<objectsnumber; i++) {
	       double t = 0;

	        if(object.objectArray[i].sphere.position && object.objectArray[i].sphere.radius){
            t = sphereintersection(Ro, Rd, object.objectArray[i].sphere.position, object.objectArray[i].sphere.radius);

          }
          else if(object.objectArray[i].plane.position && object.objectArray[i].plane.normal){
            t = planeintersection(Ro, Rd, object.objectArray[i].plane.position, object.objectArray[i].plane.normal);

          }
          if (t>0 && t < best_t){
            best_t = t;
            helper = i;
          }
        }
        if (best_t > 0 && best_t != INFINITY) {
          if (object.objectArray[helper].sphere.position && object.objectArray[helper].sphere.radius){
            shade_pixel(object.objectArray[helper].color, y, x, image);
          }
          if (object.objectArray[helper].plane.position && object.objectArray[helper].plane.normal){
            shade_pixel(object.objectArray[helper].color, y, x, image);
          }
  	       printf("#");
         } else {
  	        printf(".");
          }
        }
        printf("\n");
      }

      FILE* output = fopen(argv[4], "w+");

      ppmprint(image, output, ppmmagicnumber);
}
