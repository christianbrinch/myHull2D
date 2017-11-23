#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

int npoints = 200;

typedef struct {
  int id;
  double x[3];
} point;

typedef struct {
  double r[3];
} edge;

typedef struct {
  point *pointID[3];
  edge edge[3];
  double normal[3];
  int isLowerMFace;
} face;

typedef struct {
  face *facet;
  size_t used;
  size_t size;
} Container;





double
dot(double a[3], double b[3]){
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void
cross(double a[3], double b[3], double *normal){
  normal[0] =  (a[1]*b[2] - b[1]*a[2]);
  normal[1] = -(a[0]*b[2] - b[0]*a[2]);
  normal[2] =  (a[0]*b[1] - b[0]*a[1]);
  return;
}

void
greatCircle(face facet, double *c){
  int i;
  double qxy[3], qx[3], qy[3], det[4];
  for(i = 0; i < 3; i++){
    qxy[i] = pow(facet.pointID[i]->x[0], 2) + pow(facet.pointID[i]->x[1], 2);
    qx[i] = facet.pointID[i]->x[0];
    qy[i] = facet.pointID[i]->x[1];
  }

  det[0] = qx[0]*qy[1] + qy[0]*qx[2] + qx[1]*qy[2] - qy[1]*qx[2] - \
           qy[0]*qx[1] - qx[0]*qy[2];
  det[1] = qxy[0]*qy[1] + qy[0]*qxy[2] + qxy[1]*qy[2] - qy[1]*qxy[2] - \
           qy[0]*qxy[1] - qxy[0]*qy[2];
  det[2] = qxy[0]*qx[1] + qx[0]*qxy[2] + qxy[1]*qx[2] - qx[1]*qxy[2] - \
           qx[0]*qxy[1] - qxy[0]*qx[2];
  det[3] = qxy[0]*qx[1]*qy[2] + qx[0]*qy[1]*qxy[2] + qy[0]*qxy[1]*qx[2] \
          -qy[0]*qx[1]*qxy[2] - qx[0]*qxy[1]*qy[2] - qxy[0]*qy[1]*qx[2];

  c[0] = 0.5 * det[1]/det[0];
  c[1] = -0.5 * det[2]/det[0];
  c[2] = sqrt(pow(qx[0] - c[0], 2) + pow(qy[0] - c[1], 2));
}


int
checkFacet(face facet, point *P){
  int j, i, sign;
  double proj, u[3];
  // Check if new facet is an M-face
  sign = 0;
  for(j = 0; j < npoints; j++){
    if(j != facet.pointID[0]->id && j != facet.pointID[1]->id && j != facet.pointID[2]->id){
      for(i = 0; i < 3; i++){
        u[i] = P[j].x[i] - facet.pointID[0]->x[i];
      }
      proj = dot(u, facet.normal);
      if(fabs(proj) > 0.){
        proj = proj/fabs(proj);
      }
      if(sign == 0){
        sign = (int) proj;
      } else {
        if((int) proj != sign){
          return 0;
        }
      }
    }
  }
  return 1;
}

int
checkLowerHull(face facet, point *P){
  int j;
  double circ[3];
  // Check if facet belongs to lower hull
  greatCircle(facet, circ);
  for(j = 0; j < npoints; j++){
    if(j != facet.pointID[0]->id && j != facet.pointID[1]->id && j != facet.pointID[2]->id){
      if(sqrt(pow(P[j].x[0] - circ[0], 2) + pow(P[j].x[1] - circ[1], 2)) < circ[2]){
        return -1;
      }
    }
  }
  return 1;
}

int
checkForDuplicateFacet(Container *c, int f){
  int perm[npoints], i, j, t;
  // Check if facet already exists
  for (i = 0; i < f; i++){
    for(j = 0; j < npoints; j++){
      perm[j] = 0;
    }
    for(j = 0; j < 3; j++){
      perm[c->facet[i].pointID[j]->id] += 1;
    }
    for(j= 0; j < 3; j++){
      perm[c->facet[f].pointID[j]->id] -= 1;
    }
    t = 0;
    for(j = 0; j < npoints; j++){
      t += abs(perm[j]);
    }
    if(t == 0) return 1;
  }
  return 0;
}

void
defineFacet(Container *c, int f, point *p1, point *p2, point *p3){
  int k, i;
  double norm;
  // Create new facet
  if(c->used == c->size - 1){
    c->size *= 2;
    c->facet = (face *)realloc(c->facet, c->size * sizeof(face));
  }
  c->facet[f].pointID[0] = p1;
  c->facet[f].pointID[1] = p2;
  c->facet[f].pointID[2] = p3;
  c->used = f;

  for(k = 0; k < 3; k++){
    for(i = 0; i < 3; i++){
      c->facet[f].edge[k].r[i] = c->facet[f].pointID[(k+1)%3]->x[i] - c->facet[f].pointID[k]->x[i];
    }

  }
  cross(c->facet[f].edge[0].r, c->facet[f].edge[1].r, c->facet[f].normal);
  norm = sqrt(dot(c->facet[f].edge[0].r, c->facet[f].edge[0].r)*dot(c->facet[f].edge[1].r, c->facet[f].edge[1].r) - pow(dot(c->facet[f].edge[0].r, c->facet[f].edge[1].r), 2) );
  for(i = 0; i < 3; i++){
    c->facet[f].normal[i] = c->facet[f].normal[i]/norm;
  }
}




int
traverse(int f, int root, int side, Container *c, point *P){
  int j,best, tmp;
  double best_angle=0.,angle;
  tmp=f;

  for(j=0;j<npoints;j++){
    if(j!=c->facet[root].pointID[0]->id && j!=c->facet[root].pointID[1]->id && j!=c->facet[root].pointID[2]->id){
      defineFacet(c, f, &P[c->facet[root].pointID[(side+0)%3]->id], &P[c->facet[root].pointID[(side+1)%3]->id], &P[j]);
      angle = acos(dot(c->facet[root].normal,c->facet[f].normal)/(sqrt(dot(c->facet[root].normal,c->facet[root].normal))*sqrt(dot(c->facet[f].normal,c->facet[f].normal))) );
      if(angle>best_angle){
        best_angle = angle;
        best = j;
      }
    }
  }

  // Define new facet
  defineFacet(c, f, &P[c->facet[root].pointID[(side+0)%3]->id], &P[c->facet[root].pointID[(side+1)%3]->id], &P[best]);

  if(!checkFacet(c->facet[f], P)){
    printf("New facet is not M face\n");
    exit(0);
  }

  if(checkForDuplicateFacet(c, f) == 1){
    return f-1;
  }

  c->facet[f].isLowerMFace = checkLowerHull(c->facet[f],P);


  f = traverse(f+1, tmp, (side+1)%3, c, P);
  f = traverse(f+1, tmp, (side+2)%3, c, P);
  f = traverse(f+1, tmp, (side+3)%3, c, P);
  return f;
}

void
triangulate(Container *c, point *P){
  int m, n, p, f=0;

  // Lift points to higher dimensional paraboloid
  for(p = 0; p < npoints; p++){
    P[p].x[2] = pow(P[p].x[0], 2) + pow(P[p].x[1], 2);
  }
  // Construct convex hull in 2+1 dimensions
  // First we find an initial hull facet
  m = 1;
  n = 2;
  do{
    defineFacet(c, f, &P[0], &P[m], &P[n]);
    n += 1;
    if(n == npoints){
      m += 1;
      n = m + 1;
    }
  } while(!checkFacet(c->facet[f], P) && m < npoints - 1);

  if(m == npoints - 1){
    printf("Can't find a facet on the hull\n");
    exit(0);
  }
  c->facet[f].isLowerMFace = checkLowerHull(c->facet[0], P);
  // initial M-face found

  f=traverse(f+1, 0, 0, c, P);
}


void initContainer(Container *c, size_t np) {
  c->facet = (face *)malloc(np * sizeof(face));
  c->used = 0;
  c->size = np;
}

void freeContainer(Container *c){
  free(c->facet);
  c->facet = NULL;
  c->used = c->size = 0;
}


int
main(){
  gsl_rng *ran = gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(ran,time(0));
  point P[npoints];
  Container c;
  double dist[3];
  int i,f,j,k,min;

  // Build set P={p1,p2,...pn}
  for(i=0;i<npoints;i++){
    P[i].id = i;
    for(j=0;j<2;j++){
      P[i].x[j]=gsl_ran_gaussian(ran,0.4);
      //P[i].x[j]=(2*gsl_rng_uniform(ran))-1.;
    }
  }

  initContainer(&c,npoints);
  triangulate(&c, P);

  // Lloyd smoothing
  for(i=0;i<25;i++){
    for(f=0;f<c.used;f++){
      if(c.facet[f].isLowerMFace==1){
        for(k=0;k<3;k++){
           dist[k]= sqrt(dot(c.facet[f].edge[k].r,c.facet[f].edge[k].r));
        }
        if(dist[0] < dist[1] && dist[0] < dist[2]){
          min = 0;
        } else if (dist[1] < dist[0] && dist[0] < dist[2]){
          min = 1;
        } else {
          min = 2;
        }
        for(j=0;j<3;j++){
          P[c.facet[f].pointID[(min+1)%3]->id].x[j] += c.facet[f].edge[min].r[j]/dist[min]*0.10*((dist[(min+1)%3]/2.+dist[(min+2)%3]/2.)-dist[min]);
          P[c.facet[f].pointID[min]->id].x[j] -= c.facet[f].edge[min].r[j]/dist[min]*0.10*((dist[(min+1)%3]/2.+dist[(min+2)%3]/2.)-dist[min]);
        }

      }
    }
    freeContainer(&c);
    initContainer(&c,npoints);
    triangulate(&c, P);
  }


  for(j=0;j<npoints;j++){
    printf("%f %f %f\n", P[j].x[0],P[j].x[1], P[j].x[2]);
  }
  for(f=0;f<c.used;f++){
    if(c.facet[f].isLowerMFace==1){
      printf("%d %d %d\n", c.facet[f].pointID[0]->id, c.facet[f].pointID[1]->id, c.facet[f].pointID[2]->id);
    }
  }
  freeContainer(&c);
  gsl_rng_free(ran);
}
