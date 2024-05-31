#include "FPToolkit.c"
#include "M3d_matrix_tools.c"
#include <stdio.h>
#include <stdlib.h>


#define MAXOBJS 50
#define MAXPTS 59000
#define MAXPOLYS 57500

int numobjects;
int numpoints[MAXOBJS];
int numpolys[MAXOBJS];
double x[MAXOBJS][MAXPTS];
double y[MAXOBJS][MAXPTS];
double z[MAXOBJS][MAXPTS];
double cx[MAXOBJS],cy[MAXOBJS],cz[MAXOBJS]; //tracks the center of the 3d object
int psize[MAXOBJS][MAXPOLYS];
int con[MAXOBJS][MAXPOLYS][20];
int rotation; //changes the orientation of the backface elimination
double red[MAXOBJS], green[MAXOBJS], blue[MAXOBJS];

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


void perspective_polygon(double xp[100], double yp[100], double zp[100], int np, int onum)
{
  double xbb[100], ybb[100];
  double H=tan(40*(M_PI/180));
  
  for(int i =0; i<np; i++){
    xbb[i]=(400/H)*(xp[i]/zp[i])+400;
    ybb[i]=(400/H)*(yp[i]/zp[i])+400;
  }
  G_rgb(red[onum], green[onum], blue[onum]);
  G_fill_polygon(xbb, ybb, np);
  //G_rgb(0, 0, 0);
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

int draw_all (int onum)
{
  int h, np, nnp, totalpolys=0;
  double A[3], B[3], P[3], E[3], dot;
  double xp[100], yp[100], zp[100], zdist=-1000;
  int k = 0 ;
  
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

  //printf("k = %d  totalpolys = %d\n",k,totalpolys) ;
  

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

    if(rotation==1){
      if(dot>0){perspective_polygon(xp, yp, zp, nnp, p[i].objnum);}
      else if (dot<0) {perspective_polygon_bf(xp, yp, zp, nnp, p[i].objnum);;}
    }
    if (rotation==-1){
      if(dot<0){perspective_polygon(xp, yp, zp, nnp, p[i].objnum);}
      else if (dot>0) {perspective_polygon_bf(xp, yp, zp, nnp, p[i].objnum);;}
    }
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
  int ke;
  
  int sign=1;
  int action='t';
  double temp[4][4] ;
  rotation=1;

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
    red[onum]=r2();
    green[onum]=r2();
    blue[onum]=r2();
    }
  //print_object(stdout, 0) ;
  
  M3d_make_identity(temp);
 
  onum = 0 ;
  G_init_graphics(800,800) ;
  while (0==0) {
    G_rgb(1, 1, 1);
    G_clear();
    G_rgb(1,0,1);
    // print_object(stdout, 0); 
    if (key == 'q') {
      exit(0) ;

    } else if (key == 'c') {
      sign = -sign ;

    } else if (key == 't') {
      action = key ;

    } else if (key == 'r') {
      action = key ;

    } else if (key == 's') {
      rotation= -rotation;
      
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
    
    draw_all(onum);
    key=G_wait_key(); 
  }  
}

