#include "FPToolkit.c"
#include "M3d_matrix_tools.c"
#include <stdio.h>
#include <stdlib.h>


#define MAXOBJS 50
#define MAXPTS 59000
#define MAXPOLYS 57500

#define MaxDiffuse .5
#define Ambient .2
#define SpecPower 50


double light_in_eyespace[3];
int numobjects;
int numpoints[MAXOBJS];
int numpolys[MAXOBJS];
double x[MAXOBJS][MAXPTS];
double y[MAXOBJS][MAXPTS];
double z[MAXOBJS][MAXPTS];
double cx[MAXOBJS],cy[MAXOBJS],cz[MAXOBJS]; //tracks the center of the 3d object
int psize[MAXOBJS][MAXPOLYS];
int con[MAXOBJS][MAXPOLYS][20];
int rotation[MAXOBJS]; //changes the orientation of the backface elimination
double red[MAXOBJS], green[MAXOBJS], blue[MAXOBJS];
int s=1;

typedef
struct {
  int objnum;
  int polynum;
  double dist;
}
Dt;
Dt p[MAXPOLYS];

double r2()
{
  return (double)rand() / (double)RAND_MAX;
}


void translate(int onum, double dx, double dy, double dz){
  double temp[4][4];

  M3d_make_translation(temp, dx, dy, dz);
  M3d_mat_mult_points(x[onum],y[onum],z[onum],temp,x[onum],y[onum],z[onum], numpoints[onum]+1);

  cx[onum]=cx[onum]+dx;
  cy[onum]=cy[onum]+dy;
  cz[onum]=cz[onum]+dz;
}


int read_object(FILE *f, int onum)
{
  int i,j ;
    // point info
    fscanf(f,"%d",&numpoints[onum]);
    if (numpoints[onum] >= MAXPTS) {
      // need an extra for object centering
      printf("MAXPTS = %d :  exceeded.\n",MAXPTS);
      exit(1);
    }
    for (i = 0 ; i < numpoints[onum] ; i++) {
      fscanf(f,"%lf %lf %lf",&x[onum][i], &y[onum][i], &z[onum][i]) ;
    }
    // connectivity info
    fscanf(f,"%d",&numpolys[onum]);
    if (numpolys[onum] > MAXPOLYS) {
      printf("MAXPOLYS = %d :  exceeded.\n",MAXPOLYS) ;
      exit(1) ;
    }
    for (i = 0 ; i < numpolys[onum] ; i++) {
      fscanf(f,"%d",&psize[onum][i]) ;
      for (j = 0 ; j < psize[onum][i] ; j++) {
        fscanf(f,"%d",&con[onum][i][j]) ;
      } // end for j
    } // end for i
      
}

int print_object (FILE *fout, int onum)
{
  int i,j ;
  fprintf(fout, "%d\n",numpoints[onum]) ;
  for (i = 0 ; i < numpoints[onum] ; i++) {
    fprintf(fout, "%12.6lf %12.6lf %12.6lf\n",x[onum][i],y[onum][i],z[onum][i]) ;
  }
  for (i = 0 ; i < numpolys[onum] ; i++) {
    fprintf(fout, "%3d    ",psize[onum][i]) ;
    for (j = 0 ; j < psize[onum][i] ; j++) {
      fprintf(fout, "%2d ", con[onum][i][j]) ;
    }
    fprintf(fout, "\n") ;
  }    
}

void poly_center(int onum, int poly, double center[]){
	double avx = 0,avy = 0, avz = 0;
	for(int i = 0; i < psize[onum][poly]; i++){
		avz += z[onum][con[onum][poly][i]];
		avy += y[onum][con[onum][poly][i]];
		avx += x[onum][con[onum][poly][i]];
	}
	center[0] = avx/psize[onum][poly];
	center[1] = avy/psize[onum][poly];
	center[2] = avz/psize[onum][poly];
}

void normal(double vector[]){
	double Length = sqrt(pow(vector[0],2)+pow(vector[1],2)+pow(vector[2],2));
	vector[0] /= Length;
	vector[1] /= Length;
	vector[2] /= Length;
}

void colorchange(double rgb[],double intensity, double r, double g, double b){
	rgb[0] = r;
	rgb[1] = g;
	rgb[2] = b;
	double tot = Ambient + MaxDiffuse;
	if(intensity == tot){
		return;
	}else if(intensity < tot){
		rgb[0] = rgb[0]*(intensity/tot);
		rgb[1] = rgb[1]*(intensity/tot);
		rgb[2] = rgb[2]*(intensity/tot);
	}else if(intensity > tot){
		rgb[0] = rgb[0] + (1-rgb[0])*((intensity - tot)/(1-tot));
		rgb[1] = rgb[1] + (1-rgb[1])*((intensity - tot)/(1-tot));
		rgb[2] = rgb[2] + (1-rgb[2])*((intensity - tot)/(1-tot));
	}
}

double light_model(int onum, int poly){
  double Lu[3], Ru[3], center[3];
  double a[3], b[3], Nu[3];

  a[0] = x[onum][con[onum][poly][1]] - x[onum][con[onum][poly][0]];
  a[1] = y[onum][con[onum][poly][1]] - y[onum][con[onum][poly][0]];
  a[2] = z[onum][con[onum][poly][1]] - z[onum][con[onum][poly][0]];
	
  b[0] = x[onum][con[onum][poly][2]] - x[onum][con[onum][poly][0]];
  b[1] = y[onum][con[onum][poly][2]] - y[onum][con[onum][poly][0]];
  b[2] = z[onum][con[onum][poly][2]] - z[onum][con[onum][poly][0]];
  M3d_x_product(Nu,a,b);

  poly_center(onum,poly,center);
  Lu[0] = light_in_eyespace[0] - center[0];
  Lu[1] = light_in_eyespace[1] - center[1];
  Lu[2] = light_in_eyespace[2] - center[2];
  
  normal(Lu);
  normal(Nu);

  Nu[0]*=rotation[onum]; Nu[1]*=rotation[onum]; Nu[2]*=rotation[onum];
  
  double dotln=(Lu[0]*Nu[0] + Lu[1]*Nu[1] + Lu[2]*Nu[2]);
  double diffuse = MaxDiffuse*dotln;
  
  if(dotln < 0){diffuse=0;}
  Ru[0]=((2*dotln)*Nu[0])-Lu[0];
  Ru[1]=((2*dotln)*Nu[1])-Lu[1];
  Ru[2]=((2*dotln)*Nu[2])-Lu[2];
  double Eu[3] = {-center[0],-center[1],-center[2]};
  normal(Eu);
  
  double doter=(Eu[0]*Ru[0] + Eu[1]*Ru[1] + Eu[2]*Ru[2]);
  double spec = (1-Ambient-MaxDiffuse)*pow(doter,SpecPower);
  if(doter<0){spec=0;}
  return spec + diffuse;
}

void perspective_polygon(double xp[100], double yp[100], double zp[100], int np, int onum, int poly, double color[3])
{
  double xbb[100], ybb[100], colorun[3];
  double H=tan(40*(M_PI/180));
  double intensity;
  double t1 = 0.2, t2 = 0.5, t3 = 0.8;
  
  for(int i =0; i<np; i++){
    xbb[i]=(400/H)*(xp[i]/zp[i])+400;
    ybb[i]=(400/H)*(yp[i]/zp[i])+400;
  }
  intensity = Ambient + light_model(onum, poly);
  /*
  if (intensity<=t1) {G_rgb(0, 0, 0);}
  else if (intensity>t1 && intensity<t2) {G_rgb(0.4, 0.4, 0.4);}
  else {G_rgb(0.8, 0.8, 0.8);} */
  
  if(onum == 1){
    if (intensity<=t1) {G_rgb(color[0]-color[0],color[1]-color[1],color[2]-color[2]);}
    else if (intensity>t1 && intensity<t2) {G_rgb(color[0]-(color[0]/2),color[1]-(color[1]/2),color[2]-(color[2]/2));}
    else {G_rgb(color[0],color[1],color[2]);}
  } 
  
  else if(onum == 0){
    if (intensity<=t1) {G_rgb(color[0]-color[0],color[1]-color[1],color[2]-color[2]);}
    else if (intensity>t1 && intensity<t2) {G_rgb(color[0]-(color[0]/2),color[1]-(color[1]/2),color[2]-(color[2]/2));}
    else {G_rgb(color[0],color[1],color[2]);}
    //G_rgb(color[0],color[1],color[2]);
  } 
  
  else {
    if (intensity<=t1) {G_rgb(color[0]-color[0],color[1]-color[1],color[2]-color[2]);}
    else if (intensity>t1 && intensity<t2) {G_rgb(color[0]-(color[0]/2),color[1]-(color[1]/2),color[2]-(color[2]/2));}
    else {G_rgb(color[0],color[1],color[2]);}
    //colorchange(color,intensity,.4,.3,.2);
    //G_rgb(color[0],color[1],color[2]);
  }
 
  //G_rgb(red[onum], green[onum], blue[onum]);
  G_fill_polygon(xbb, ybb, np);
  G_rgb(0,0,0);
  //G_polygon(xbb,ybb,np);
   
}

void perspective_polygon_bf(double xp[100], double yp[100], double zp[100], int np, int onum)
{
  double xbb[100], ybb[100];
  double H=tan(40*(M_PI/180));
  G_rgb(0, 0, 0);
  for (int j=-2; j<3; j++) { 
    for (int k=-2; k<3; k++) { 
      for(int i =0; i<np; i++){
        xbb[i]=(400/H)*(xp[i]/zp[i])+400+j;
        ybb[i]=(400/H)*(yp[i]/zp[i])+400+k;
      }
    G_fill_polygon(xbb, ybb, np);
    }
  }
}
  
int compare (const void *p, const void *q)
{
  Dt *a, *b ;

  a = (Dt*)p ;
  b = (Dt*)q ;

  if  (((*a).dist) < ((*b).dist)) return -1 ;
  else if (((*a).dist) > ((*b).dist)) return 1 ;
  else return 0 ;
}

int draw_all (int onum, double r, double g, double b)
{
  int h, np, nnp, totalpolys=0;
  double A[3], B[3], P[3], E[3], dot;
  double xp[100], yp[100], zp[100], zdist=-1000;
  int k = 0 ;
  double color[3];
  color[0]=r; color[1]=g; color[2]=b;
  
  for(onum=0; onum < numobjects; onum++){
    totalpolys+=numpolys[onum];
    for (int i = 0 ; i < numpolys[onum] ; i++) {
      np = psize[onum][i] ;

      
      zdist = 0 ;
      for (int j = 0 ; j < np ; j++) {
	h = con[onum][i][j];
        zdist += z[onum][h] ;
      }
      zdist /= np ;

      
      p[k].objnum=onum;
      p[k].polynum=i;
      p[k].dist=zdist ;
      k++ ;
    }
  }

  //printf("totalpolys = %d\n", totalpolys);  
  qsort(p, totalpolys, sizeof(Dt), compare);
  // for(int i=0; i<totalpolys; i++){printf("%lf\n", zdist);}
  //sorted from closest to farthest
  //need to draw in reverse (totalpolys-1 to 0)
  //by drawing all the polys the objects should still be complete if done correctly 
  
  for(int i=totalpolys-1; i>=0; i--){
    nnp = psize[p[i].objnum][p[i].polynum];
    for (int j = 0 ; j < nnp ; j++) {
       h=con[p[i].objnum][p[i].polynum][j];
       xp[j] = x[p[i].objnum][h];
       yp[j] = y[p[i].objnum][h];
       zp[j] = z[p[i].objnum][h];
      }
    E[0]=-xp[0]; E[1]=-yp[0]; E[2]=-zp[0];
    A[0]=xp[1]-xp[0]; A[1]=yp[1]-yp[0]; A[2]=zp[1]-zp[0];
    B[0]=xp[2]-xp[0]; B[1]=yp[2]-yp[0]; B[2]=zp[2]-zp[0];
    M3d_x_product(P,A,B);
    dot=(P[0]*E[0])+(P[1]*E[1])+(P[2]*E[2]);

    if(s==1){
      if (dot<0) {perspective_polygon_bf(xp, yp, zp, nnp, p[i].objnum);}
    }
    if (s==-1){
      if (dot>0) {perspective_polygon_bf(xp, yp, zp, nnp, p[i].objnum);}
    }
    
    perspective_polygon(xp, yp, zp, nnp, p[i].objnum, p[i].polynum, color);
  }
}

void center (int onum) {
  double ax, ay, az;
  double lx, hx, ly, hy, lz, hz;
  double s, temp[4][4];

  lx = hx = x[onum][0];
  ly = hy = y[onum][0];
  lz = hz = z[onum][0];
  for(int i =0; i <numpoints[onum]; i++)
    {
    if (lx>x[onum][i]) {lx=x[onum][i];}
    if (hx<x[onum][i]) {hx=x[onum][i];}
    if (ly>y[onum][i]) {ly=y[onum][i];}
    if (hy<y[onum][i]) {hy=y[onum][i];}
    if (lz>z[onum][i]) {lz=z[onum][i];}
    if (hz<z[onum][i]) {hz=z[onum][i];}
    }
  ax=hx-lx;
  ay=hy-ly;
  az=hz-lz;
  cx[onum]=(ax/2)+lx;
  cy[onum]=(ay/2)+ly;
  cz[onum]=(az/2)+lz;

  translate(onum,-cx[onum],-cy[onum],-cz[onum]);
  //x[onum][numpoints[onum]]=0;
  //y[onum][numpoints[onum]]=0;
  //z[onum][numpoints[onum]]=0;
}

int main(int argc, char **argv)
{
  FILE *fln ;
  int key,w ;
  char fname[100] ;
  int onum;
  int count=0;
  int ke;
  
  
  int sign=1;
  int action='t';
  double temp[4][4] ;
  
  //s=1;

  //  scanf("%d", &numobjects);
  numobjects = argc - 1 ;
  
  if (numobjects>MAXOBJS) {
    printf("MAXOBJS= %d : exceeded.\n", MAXOBJS);
    exit(1);
  }
  for (onum=0; onum < numobjects; onum++) { 
    //    printf("enter name of xy file ") ;
    //    scanf("%s",fname) ;
    //    fln = fopen(fname,"r") ;
    fln = fopen(argv[onum+1],"r") ;    
    if (fln == NULL) {
      printf("can't read file, %s\n",fname) ;
      exit(1) ;
    }
    read_object(fln, onum) ;
    center(onum);
    rotation[onum]=1;
    //red[onum]=r2();
    //green[onum]=r2();
    //blue[onum]=r2();
    }

    printf("enter your light coords:\n");
    printf("x = ");
    scanf("%lf",&light_in_eyespace[0]);
    printf("\ny = ");
    scanf("%lf", &light_in_eyespace[1]);
    printf("\nz = ");
    scanf("%lf",&light_in_eyespace[2]);
  
  
  M3d_make_identity(temp);
 
  onum = 0 ;
  G_init_graphics(800,800) ;
  while (0==0) {
    G_rgb(1, 1, 1);
    G_clear();
    G_rgb(1,0,1); 
    if (key == 'q') {
      exit(0) ;

    } else if (key == 'c') {
      sign = -sign ;

    } else if (key == 't') {
      action = key ;

    } else if (key == 'r') {
      action = key ;

    } else if (key == 's') {
      rotation[onum]= -rotation[onum];
      s = -s;
      
    } else if (('0' <= key) && (key <= '9')) {
      w = key - '0' ;  
      if (w < numobjects) { onum = w ; }

    } else if ((key == 'x') && (action == 't')) {
      translate(onum, sign*2, 0, 0);
    } else if ((key == 'y') && (action == 't')) {
      translate(onum, 0, sign*2, 0);
    } else if ((key == 'z') && (action == 't')) {
      translate(onum, 0, 0, sign*2);
    } else if ((key == 'x') && (action == 'r')) {
      M3d_make_translation(temp, -cx[onum],-cy[onum],-cz[onum]);
      M3d_mat_mult_points(x[onum],y[onum],z[onum],temp,x[onum],y[onum],z[onum], numpoints[onum]+1);
      M3d_make_x_rotation_cs(temp, cos(sign*(2*M_PI/180)), sin(sign*(2*M_PI/180)));
      M3d_mat_mult_points(x[onum],y[onum],z[onum],temp,x[onum],y[onum],z[onum], numpoints[onum]+1);
      M3d_make_translation(temp,cx[onum],cy[onum],cz[onum]);
      M3d_mat_mult_points(x[onum],y[onum],z[onum],temp,x[onum],y[onum],z[onum], numpoints[onum]+1);
      
    } else if ((key == 'y') && (action == 'r')) {
      M3d_make_translation(temp, -cx[onum],-cy[onum],-cz[onum]);
      M3d_mat_mult_points(x[onum],y[onum],z[onum],temp,x[onum],y[onum],z[onum], numpoints[onum]+1);
      M3d_make_y_rotation_cs(temp, cos(sign*(2*M_PI/180)), sin(sign*(2*M_PI/180)));
      M3d_mat_mult_points(x[onum],y[onum],z[onum],temp,x[onum],y[onum],z[onum], numpoints[onum]+1);
      M3d_make_translation(temp,cx[onum],cy[onum],cz[onum]);
      M3d_mat_mult_points(x[onum],y[onum],z[onum],temp,x[onum],y[onum],z[onum], numpoints[onum]+1);
      
    } else if ((key == 'z') && (action == 'r')) {
      M3d_make_translation(temp, -cx[onum],-cy[onum],-cz[onum]);
      M3d_mat_mult_points(x[onum],y[onum],z[onum],temp,x[onum],y[onum],z[onum], numpoints[onum]+1);
      M3d_make_z_rotation_cs(temp, cos(sign*(2*M_PI/180)), sin(sign*(2*M_PI/180)));
      M3d_mat_mult_points(x[onum],y[onum],z[onum],temp,x[onum],y[onum],z[onum], numpoints[onum]+1);
      M3d_make_translation(temp,cx[onum],cy[onum],cz[onum]);
      M3d_mat_mult_points(x[onum],y[onum],z[onum],temp,x[onum],y[onum],z[onum], numpoints[onum]+1);
      
    } else {
      printf("no action\n") ;
    }
    
    draw_all(onum, 0.75, 0.75, 0.75);
    char fname[100];
    sprintf(fname, "hh%04d.xwd", count);
    G_save_image_to_file(fname) ;
    count++;
    printf("%d ", count);
    key=G_wait_key(); 
  }  
}
