#include <stdio.h>
#include <math.h>


/*

 ( x')          (x)
 ( y')  =   M * (y)  
 ( z')          (z)
 ( 1 )          (1)

instead of (x',y',z',1) = (x,y,z,1) * M  

*/




int M3d_print_mat (double a[4][4])
{
  int r,c ;
  for (r = 0 ; r < 4 ; r++ ) {
      for (c = 0 ; c < 4 ; c++ ) {
           printf(" %12.4lf ",a[r][c]) ;
      }
      printf("\n") ;
  }

  return 1 ;
} 





int M3d_copy_mat (double a[4][4], double b[4][4])
// a = b
{
  int r,c ;
  for (r = 0 ; r < 4 ; r++ ) {
      for (c = 0 ; c < 4 ; c++ ) {
           a[r][c] = b[r][c] ;
      }
  }

  return 1 ;
} 





int M3d_make_identity (double a[4][4])
// a = I
{
  int r,c ;
  for (r = 0 ; r < 4 ; r++ ) {
      for (c = 0 ; c < 4 ; c++ ) {
           if (r == c) a[r][c] = 1.0 ;
               else    a[r][c] = 0.0 ;
      }
  }

  return 1 ;
} 





int M3d_make_translation (double a[4][4], double dx, double dy, double dz)
{
  M3d_make_identity(a) ;
  a[0][3] =  dx ;  a[1][3] = dy ;  a[2][3] = dz ;
  return 1 ;
}





int M3d_make_scaling (double a[4][4], double sx, double sy, double sz)
{
  M3d_make_identity(a) ;
  a[0][0] =  sx ;  a[1][1] = sy ;  a[2][2] = sz ;
  return 1 ;
}












int M3d_make_x_rotation_cs (double a[4][4], double cs, double sn)
// this one assumes cosine and sine are already known
{
  M3d_make_identity(a) ;

  a[1][1] =   cs ;  a[1][2] = -sn ;
  a[2][1] =   sn ;  a[2][2] =  cs ;

  return 1 ;
}



int M3d_make_y_rotation_cs (double a[4][4], double cs, double sn)
// this one assumes cosine and sine are already known
{
  M3d_make_identity(a) ;

  a[0][0] =   cs ;  a[0][2] =  sn ;
  a[2][0] =  -sn ;  a[2][2] =  cs ;

  return 1 ;
}


int M3d_make_z_rotation_cs (double a[4][4], double cs, double sn)
// this one assumes cosine and sine are already known
{
  M3d_make_identity(a) ;

  a[0][0] =   cs ;  a[0][1] = -sn ;
  a[1][0] =   sn ;  a[1][1] =  cs ;

  return 1 ;
}





int M3d_mat_mult (double res[4][4], double a[4][4], double b[4][4])
// res = a * b
// this is SAFE, i.e. the user can make a call such as 
// M3d_mat_mult(p,  p,q) or M3d_mat_mult(p,  q,p) or  M3d_mat_mult(p, p,p)
{
  double sum ;
  int k ;
  int r,c ;
  double tmp[4][4] ;

  for (r = 0 ; r < 4 ; r++ ) {
      for (c = 0 ; c < 4 ; c++ ) {
           sum = 0.0 ;
           for (k = 0 ; k < 4 ; k++) {
                 sum = sum + a[r][k]*b[k][c] ;
           }
           tmp[r][c] = sum ;
      }
  }


  M3d_copy_mat (res,tmp) ;

  return 1 ;
}





int M3d_mat_mult_pt (double P[3],   double m[4][4], double Q[3])
// P = m*Q
// SAFE, user may make a call like M3d_mat_mult_pt (W, m,W) ;
{
  double u,v,t ;

  u = m[0][0]*Q[0] + m[0][1]*Q[1] + m[0][2]*Q[2] + m[0][3] ;
  v = m[1][0]*Q[0] + m[1][1]*Q[1] + m[1][2]*Q[2] + m[1][3] ;
  t = m[2][0]*Q[0] + m[2][1]*Q[1] + m[2][2]*Q[2] + m[2][3] ;  

  P[0] = u ;
  P[1] = v ;
  P[2] = t ;
  
  return 1 ;
}





int M3d_mat_mult_points (double X[], double Y[], double Z[],
                         double m[4][4],
                         double x[], double y[], double z[], int numpoints)
// |X0 X1 X2 ...|       |x0 x1 x2 ...|
// |Y0 Y1 Y2 ...| = m * |y0 y1 y2 ...|
// |Z0 Z1 Z2 ...|       |z0 z1 z2 ...|  
// | 1  1  1 ...|       | 1  1  1 ...|

// SAFE, user may make a call like M3d_mat_mult_points (x,y,z,  m, x,y,z,  n) ;
{
  double u,v,t ;
  int i ;

  for (i = 0 ; i < numpoints ; i++) {
    u = m[0][0]*x[i] + m[0][1]*y[i] + m[0][2]*z[i] + m[0][3] ;
    v = m[1][0]*x[i] + m[1][1]*y[i] + m[1][2]*z[i] + m[1][3] ;
    t = m[2][0]*x[i] + m[2][1]*y[i] + m[2][2]*z[i] + m[2][3] ;    

    X[i] = u ;
    Y[i] = v ;
    Z[i] = t ;
  }

  return 1 ;
}






int M3d_x_product (double res[3], double a[3], double b[3])
// res = a x b  , cross product of two vectors
// SAFE: it is ok to make a call such as
// D3d_x_product (a,  a,b) or
// D3d_x_product (b,  a,b) or
// D3d_x_product (a,  a,a) 
{
    double r[3] ;
    int v ;
    
    r[0] = a[1]*b[2] - b[1]*a[2] ;
    r[1] = b[0]*a[2] - a[0]*b[2] ;
    r[2] = a[0]*b[1] - b[0]*a[1] ;

    res[0] = r[0] ;
    res[1] = r[1] ;
    res[2] = r[2] ;

    if ((res[0] == 0) && (res[1] == 0) && (res[2] == 0)) {
	v = 0 ;
    } else {
	v = 1 ;
    }

    return v ;
}






//===========================================================================
// For Advanced Graphics :
//===========================================================================

#define SX 0
#define SY 1
#define SZ 2

#define RX 3
#define RY 4
#define RZ 5

#define TX 6
#define TY 7
#define TZ 8

#define NX 9
#define NY 10
#define NZ 11

int M3d_make_movement_sequence_matrix(double v[4][4], double vi[4][4], int n, int mtype[], double mparam[]) {
  double matrix[4][4] ;
  double iMatrix[4][4] ;
  double radians, radians2 = 0 ;
  int j = n - 1 ;
  M3d_make_identity(v) ;
  M3d_make_identity(vi) ;
  M3d_make_identity(matrix) ;
  
  for(int i = 0; i < n; i++) {
    radians = (M_PI/180)*mparam[i] ;
    radians2 = (M_PI/180)*mparam[j] ;
    printf("mtype[%d] = %d\n", i, mtype[i]) ;

    // Scaling
    if(mtype[i] == 0) { M3d_make_scaling(matrix, mparam[i], 1, 1) ; }
    else if(mtype[i] == 1) { M3d_make_scaling(matrix, 1, mparam[i], 1) ; }
    else if(mtype[i] == 2) { M3d_make_scaling(matrix, 1, 1, mparam[i]) ; }

    // Rotation
    else if(mtype[i] == 3) {M3d_make_x_rotation_cs(matrix, cos(radians), sin(radians)) ; }
    else if(mtype[i] == 4) { M3d_make_y_rotation_cs(matrix, cos(radians), sin(radians)) ; }
    else if(mtype[i] == 5) { M3d_make_z_rotation_cs(matrix, cos(radians), sin(radians)) ; }

    // Translation
    else if(mtype[i] == 6) { M3d_make_translation(matrix, mparam[i], 0, 0) ; }
    else if(mtype[i] == 7) { M3d_make_translation(matrix, 0, mparam[i], 0) ; }
    else if(mtype[i] == 8) { M3d_make_translation(matrix, 0, 0, mparam[i]) ; }

    // Negation
    else if(mtype[i] == 9) { M3d_make_scaling(matrix, -1, 1, 1) ; }
    else if(mtype[i] == 10) { M3d_make_scaling(matrix, 1, -1, 1) ; }
    else if(mtype[i] == 11) { M3d_make_scaling(matrix, 1, 1, -1) ; }

    // Scaling Inverse
    if(mtype[j] == 0) { M3d_make_scaling(iMatrix, 1/mparam[j], 1, 1) ; }
    else if(mtype[j] == 1) { M3d_make_scaling(iMatrix, 1, 1/mparam[j], 1) ; }
    else if(mtype[j] == 2) { M3d_make_scaling(iMatrix, 1, 1, 1/mparam[j]) ; }

    // Rotation inverse
    else if(mtype[j] == 3) {M3d_make_x_rotation_cs(iMatrix, cos(radians2), -sin(radians2)) ; }
    else if(mtype[j] == 4) { M3d_make_y_rotation_cs(iMatrix, cos(radians2), -sin(radians2)) ; }
    else if(mtype[j] == 5) { M3d_make_z_rotation_cs(iMatrix, cos(radians2), -sin(radians2)) ; }

    // Translation Inverse
    else if(mtype[j] == 6) { M3d_make_translation(iMatrix, -mparam[j], 0, 0) ; }
    else if(mtype[j] == 7) { M3d_make_translation(iMatrix, 0, -mparam[j], 0) ; }
    else if(mtype[j] == 8) { M3d_make_translation(iMatrix, 0, 0, -mparam[j]) ; }

    // Negation Inverse
    else if(mtype[j] == 9) { M3d_make_scaling(iMatrix, -1, 1, 1) ; }
    else if(mtype[j] == 10) { M3d_make_scaling(iMatrix, 1, -1, 1) ; }
    else if(mtype[j] == 11) { M3d_make_scaling(iMatrix, 1, 1, -1) ; }

    M3d_mat_mult(v, matrix, v) ;
    M3d_mat_mult(vi, iMatrix, vi) ;
    j-- ;
  }
  
  return 1 ;
}

int M3d_view(double v[4][4], double vi[4][4], double eye[3], double coi[3], double up[3]) {
  /*int mtype[6];
  double mparam[6], temp[4][4];
  double hyp, angle;
  
  //translation to center
  mtype[0]=TX; mparam[0]=-eye[0];
  mtype[1]=TY; mparam[1]=-eye[1];
  mtype[2]=TZ; mparam[2]=-eye[2];

  M3d_make_translation(temp, mparam[0], mparam[1], mparam[2]);
  M3d_mat_mult_pt(coi, temp, coi); M3d_mat_mult_pt(up, temp, up);

  //rotate y
  hyp=sqrt(coi[0]*coi[0]+coi[2]*coi[2]);
  angle=coi[0]/hyp;
  mtype[3]=RY; mparam[3]=-asin(angle);

  M3d_make_y_rotation_cs(temp, cos(mparam[3]), -sin(mparam[3]));
  M3d_mat_mult_pt(coi, temp, coi); M3d_mat_mult_pt(up, temp, up);

  //rotate x
  hyp=sqrt(coi[1]*coi[1]+coi[2]*coi[2]);
  angle=coi[1]/hyp;
  mtype[4]=RX; mparam[4]=asin(angle);

  M3d_make_x_rotation_cs(temp, cos(mparam[4]), sin(mparam[4]));
  M3d_mat_mult_pt(coi, temp, coi); M3d_mat_mult_pt(up, temp, up);

  //rotate z
  hyp=sqrt(up[0]*up[0]+up[1]*up[1]);
  angle=up[0]/hyp;
  mtype[5]=RZ; mparam[5]=asin(angle);
  
  M3d_make_z_rotation_cs(temp, cos(mparam[5]), sin(mparam[5]));
  M3d_mat_mult_pt(up, temp, up);

  M3d_make_movement_sequence_matrix(v, vi, 6, mtype, mparam);*/
  
  int mtype[6];
  double mparam[6];
  double temp1[4][4], temp2[4][4], temp3[4][4], temp4[4][4];
  double itemp1[4][4], itemp2[4][4], itemp3[4][4], itemp4[4][4];
  double hyp;

  M3d_make_translation(temp1, -eye[0], -eye[1], -eye[2]); M3d_make_translation(itemp1, eye[0], eye[1], eye[2]);
  M3d_mat_mult_pt(coi, temp1, coi); M3d_mat_mult_pt(up, temp1, up);
  
  hyp = sqrt(coi[0]*coi[0] + coi[2]*coi[2]);

  M3d_make_y_rotation_cs(temp2, coi[2]/hyp, -coi[0]/hyp); M3d_make_y_rotation_cs(itemp2, coi[2]/hyp, -(-coi[0]/hyp));
  M3d_mat_mult_pt(coi, temp2, coi); M3d_mat_mult_pt(up, temp2, up);

  hyp = sqrt(coi[1]*coi[1] + coi[2]*coi[2]);

  M3d_make_x_rotation_cs(temp3, coi[2]/hyp, coi[1]/hyp); M3d_make_x_rotation_cs(itemp3, coi[2]/hyp, -(coi[1]/hyp));
  M3d_mat_mult_pt(coi, temp3, coi); M3d_mat_mult_pt(up, temp3, up);

  hyp = sqrt(up[0]*up[0] + up[1]*up[1]);

  M3d_make_z_rotation_cs(temp4, up[1]/hyp, up[0]/hyp); M3d_make_z_rotation_cs(itemp4, up[1]/hyp, -(up[0]/hyp));
  M3d_mat_mult_pt(up, temp4, up);

  M3d_make_identity(v);
  M3d_make_identity(vi);
  
  M3d_mat_mult(v, temp1, v); M3d_mat_mult(v, temp2, v); M3d_mat_mult(v, temp3, v); M3d_mat_mult(v, temp4, v);
  M3d_mat_mult(vi, itemp4, vi); M3d_mat_mult(vi, itemp3, vi); M3d_mat_mult(vi, itemp2, vi); M3d_mat_mult(vi, itemp1, vi);
  
  return 1;
}

  
