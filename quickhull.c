#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>

#define npoints 2000

typedef struct {
  int id;
  double x[3];
  int used;
} point;

typedef struct {
  point p[2];
} edgeID;

typedef struct {
  double r[3];
} edge;

typedef struct {
  point *pointID[3];
  edge edge[3];
  double normal[3];
  int neigh[3];
  point *subP[npoints];
  int np;
  int isLowerHull;
  int isMFace;
  int id;
} face;

typedef struct {
  face *facet;
  size_t used;
  size_t size;
} Container;


#define STACK_MAX 100000
struct Stack {
  int     data[STACK_MAX];
  int     size;
};
typedef struct Stack Stack;

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


void Stack_Init(Stack *S)
{
  S->size = 0;
}

int Stack_Top(Stack *S)
{
  if (S->size == 0) {
    fprintf(stderr, "Error: stack empty\n");
    return -1;
  } 

  return S->data[S->size-1];
}

void Stack_Push(Stack *S, int d)
{
  if (S->size < STACK_MAX)
    S->data[S->size++] = d;
  else {
    fprintf(stderr, "Error: stack full\n");
    exit(0);
  }
}

void Stack_Pop(Stack *S)
{
  if (S->size == 0)
    fprintf(stderr, "Error: stack empty\n");
  else
    S->size--;
}



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
checkFacet(Container *c, int f, int np, point **P){
  int j, i, sign;
  double proj, u[3];
  // Check if new facet is an M-face
  sign = 0;
  for(j = 0; j < np; j++){
    if(P[j]->id != c->facet[f].pointID[0]->id && P[j]->id != c->facet[f].pointID[1]->id && P[j]->id != c->facet[f].pointID[2]->id){
      for(i = 0; i < 3; i++){
	u[i] = P[j]->x[i] - c->facet[f].pointID[0]->x[i];
      }
      proj = dot(u, c->facet[f].normal);
      if(fabs(proj) > 0){
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
checkLowerHull(face facet, int np, point **P, int f){
  int j;
  double circ[3];
  // Check if facet belongs to lower hull
  greatCircle(facet, circ);
  for(j = 0; j < np; j++){
    if(P[j]->id != facet.pointID[0]->id && P[j]->id != facet.pointID[1]->id && P[j]->id != facet.pointID[2]->id){
      if(sqrt(pow(P[j]->x[0] - circ[0], 2) + pow(P[j]->x[1] - circ[1], 2)) < circ[2]){
	return -1;
      }
    }
  }
  return 1;
}





void
defineFacet(Container *c, int f, point *P1, point *P2, point *P3){
  int k, i;
  // Create new facet
  if(c->used == c->size){
    c->size *= 2;
    c->facet = (face *)realloc(c->facet, c->size * sizeof(face));
  }
  c->facet[f].pointID[0] = P1;
  c->facet[f].pointID[1] = P2;
  c->facet[f].pointID[2] = P3;
  c->facet[f].np = 0;
  c->facet[f].id=f;
  c->used = f+1;
  for(i=0;i<3;i++){
    c->facet[f].neigh[i]=-1;
  }

  P1->used = 1;
  P2->used = 1;
  P3->used = 1;

  for(k = 0; k < 3; k++){
    for(i = 0; i < 3; i++){
      c->facet[f].edge[k].r[i] = c->facet[f].pointID[(k+1)%3]->x[i] - c->facet[f].pointID[k]->x[i];
    }

  }
  cross(c->facet[f].edge[0].r, c->facet[f].edge[1].r, c->facet[f].normal);
}

void
adjustNormal(face *facet, point *p){
  int i;
  double proj, u[3];
  for(i = 0; i < 3; i++){
    u[i] = p->x[i] - facet->pointID[0]->x[i];
  }
  proj = dot(u, facet->normal);
  if(proj<0){
    for(i=0;i<3;i++){
      facet->normal[i]=-facet->normal[i];
    }
  }
} 



void
triangulate(Container *c, point **P){
  int m, n, p, k, i, j,f=0, np,  count,ct,*freq, v, visited[20000], tmp, c1, c2;
  int index[3], best=-1, used[npoints];
  double dist, u[3], proj, proj_best,initT[3];
  Container heap;
  Stack S,visible;

  Stack_Init(&S);

  np = npoints;
  for(i=0;i<np;i++){
    used[i] = 0;
  }

  // Lift points to higher dimensional paraboloid
  for(p = 0; p < np; p++){
    P[p]->x[2] = pow(P[p]->x[0], 2) + pow(P[p]->x[1], 2);
  }

  // Construct convex hull in 2+1 dimensions
  // First we create an initial internal simplex
  initT[0]=1e30; initT[1]=0.; initT[2]=0;
  for(i = 0; i < np; i++){
    dist=sqrt(P[i]->x[0]*P[i]->x[0]+P[i]->x[1]*P[i]->x[1]+P[i]->x[2]*P[i]->x[2]);
    if(dist>initT[2]){
      initT[1]=initT[2];
      index[1]=index[2];
      initT[2]=dist;
      index[2]=P[i]->id;
    } else if (dist>initT[1] && dist<initT[2]){
      initT[1]=dist;
      index[1]=P[i]->id;
    }
    if(dist < initT[0]){
      initT[0]=dist;
      index[0]=P[i]->id;
    }
  }

  defineFacet(c, 1, P[index[0]], P[index[1]], P[index[2]]);
  proj_best = 0.;
  for(i = 0; i < np; i++){
    if(i != index[0] && i != index[1] && i != index[2]){
      for(j = 0; j < 3; j++){
	u[j] = P[i]->x[j] - c->facet[1].pointID[0]->x[j];
      }
      proj = dot(u, c->facet[1].normal);
      if(fabs(proj) > proj_best){
	proj_best = fabs(proj);
	best = i;
      }
    }
  }
  adjustNormal(&c->facet[1],P[best]); 
  defineFacet(c, 0, P[index[1]], P[index[2]], P[best]);
  adjustNormal(&c->facet[0],P[index[0]]); 
  defineFacet(c, 2, P[index[0]], P[index[1]], P[best]);
  adjustNormal(&c->facet[2],P[index[2]]); 
  defineFacet(c, 3, P[index[0]], P[index[2]], P[best]);
  adjustNormal(&c->facet[3],P[index[1]]); 
  // Here loop over three initial facets; define neighbors and outside set
  for(p=0;p<4;p++){
    for(i=0;i<3;i++){
      c->facet[p].neigh[i]=(p+i+1)%4;
    }
    c->facet[p].np=0;
    for(j=0;j<np;j++){
      if(P[j]->used==0 && used[j]==0){
	for(i = 0; i < 3; i++){
	  u[i] = P[j]->x[i] - c->facet[p].pointID[0]->x[i];
	}
	proj = dot(u, c->facet[p].normal);
	if(proj < 0 ){
	  c->facet[p].subP[c->facet[p].np++]=P[j];
	  used[j]=1;
	}
      }
    }
    if(c->facet[p].np > 0){
      Stack_Push(&S, p);
    }
  }



  // Main iteration loop
  f=3;
  while(S.size>0){
    p=Stack_Top(&S); 
    Stack_Pop(&S);

    // Find most distant point from facet
    best = -1;
    proj_best = 0; 
    for(i=0;i<c->facet[p].np; i++){
      if(c->facet[p].subP[i]->used==0){
	for(j = 0; j < 3; j++){
	  u[j] = c->facet[p].subP[i]->x[j] - c->facet[p].pointID[0]->x[j];
	}
	proj = dot(u, c->facet[p].normal);

	if(proj < proj_best){
	  proj_best = proj;
	  best = i;
	}
      }
    }
    if(best != -1){
      // best is now the most distant point


      // Now find all faces that can be seen from best
      initContainer(&heap,npoints);
      Stack_Init(&visible);
      Stack_Push(&visible, p);
      for(i=0;i<20000;i++){
	visited[i]=0;
      }
      count=0;
      heap.facet[count++] = c->facet[p];
      while(visible.size>0){
	v=Stack_Top(&visible);
	visited[c->facet[v].id]=1;
	Stack_Pop(&visible);
	for(k=0;k<3;k++){
	  if(visited[c->facet[v].neigh[k]] == 0){
	    visited[c->facet[v].neigh[k]]=1;
	    for(i = 0; i < 3; i++){
	      u[i] = c->facet[p].subP[best]->x[i] - c->facet[c->facet[v].neigh[k]].pointID[0]->x[i];
	    }
	    proj = dot(u, c->facet[c->facet[v].neigh[k]].normal);
	    if(proj<0){
	      heap.facet[count++] = c->facet[c->facet[v].neigh[k]];
	      Stack_Push(&visible, c->facet[v].neigh[k]);
	    }
	  }
	}
      }

      // Extract unique edges from heap facets.
      freq=(int *)malloc(3*count*sizeof(int));
      for(m=0;m<3*count;m++){
	freq[m]=-1;
      } 
      for(m=0;m<3*count;m++){
	ct=1;
	for(n=m+1;n<3*count;n++){

	  if( (heap.facet[n/3].pointID[(n%3)]==heap.facet[m/3].pointID[(m%3)] &&
		heap.facet[n/3].pointID[((n%3)+1)%3]==heap.facet[m/3].pointID[((m%3)+1)%3]) ||
	      (heap.facet[n/3].pointID[(n%3)]==heap.facet[m/3].pointID[((m%3)+1)%3] &&\
	       heap.facet[n/3].pointID[((n%3)+1)%3]==heap.facet[m/3].pointID[(m%3)])    ) {
	    ct++;
	    freq[n] = 0;
	  }

	}
	if(freq[m]!=0) freq[m]=ct;
      }      

      tmp=f;
      // Define new facets, determine new outside sets and link new facets to hull
      for(i=0;i<np;i++){
	used[i] = 0;
      }
      for(m=0;m<3*count;m++){
	if(freq[m]==1){
	  defineFacet(c, ++f, heap.facet[m/3].pointID[m%3], heap.facet[m/3].pointID[((m%3)+1)%3], c->facet[p].subP[best]);
	  adjustNormal(&c->facet[f],heap.facet[m/3].pointID[((m%3)+2)%3]); 
	  c1=heap.facet[m/3].pointID[m%3]->id;
	  c2=heap.facet[m/3].pointID[((m%3)+1)%3]->id;

	  for(j=0;j<3;j++){
	    for(k=0;k<3;k++){
	      if( (c->facet[heap.facet[m/3].neigh[j]].pointID[k]->id == c1 && c->facet[heap.facet[m/3].neigh[j]].pointID[(k+1)%3]->id == c2) || 
		  (c->facet[heap.facet[m/3].neigh[j]].pointID[k]->id == c2 && c->facet[heap.facet[m/3].neigh[j]].pointID[(k+1)%3]->id == c1) ){
		for(i=0;i<3;i++){
		  if(c->facet[heap.facet[m/3].neigh[j]].neigh[i] == heap.facet[m/3].id){
		    c->facet[heap.facet[m/3].neigh[j]].neigh[i] = f;
		    c->facet[f].neigh[0] = c->facet[heap.facet[m/3].neigh[j]].id;
		  }
		}
	      }
	    }
	  }

	  for(j=0;j<np;j++){
	    if(P[j]->used==0 && used[j]==0){
	      for(i = 0; i < 3; i++){
		u[i] = P[j]->x[i] - c->facet[f].pointID[0]->x[i];
	      }
	      proj = dot(u, c->facet[f].normal);
	      if(proj < 0.){
		c->facet[f].subP[c->facet[f].np++]=P[j];
		used[j]=1;
	      }
	    }
	  }
	  if(c->facet[f].np>0){
	    Stack_Push(&S, f);
	  }
	}  
      }
      free(freq);

      for(m=tmp+1;m<f+1;m++){
	ct=1;
	for(n=tmp+1;n<f+1;n++){
	  if(m!=n){
	    if(c->facet[m].pointID[0]->id==c->facet[n].pointID[0]->id || 
		c->facet[m].pointID[0]->id==c->facet[n].pointID[1]->id ||
		c->facet[m].pointID[1]->id==c->facet[n].pointID[0]->id ||
		c->facet[m].pointID[1]->id==c->facet[n].pointID[1]->id){
	      c->facet[m].neigh[ct++]=c->facet[n].id;
	    }
	  }
	}	
      }
      freeContainer(&heap);
    }
  }
  return;
}










int
main(){
  gsl_rng *ran = gsl_rng_alloc(gsl_rng_ranlxs2);
  gsl_rng_set(ran, time(0));
  point P[npoints], (*Pptr[npoints]);
  Container c;
  //double dist[3];
  int i,f,j;



  // Build set P={p1,p2,...pn}
  for(i=0;i<npoints;i++){
    Pptr[i] = &P[i];
    P[i].id = i;
    P[i].used = 0;
    for(j=0;j<2;j++){
      P[i].x[j]=gsl_ran_gaussian(ran,0.4);
      //P[i].x[j]=(2*gsl_rng_uniform(ran))-1.;
    }
  }


  initContainer(&c,npoints);
  triangulate(&c, Pptr);


/*  
  // Lloyd smoothing
  for(i=0;i<25;i++){
    for(f=0;f<c.used;f++){
      if(c.facet[f].isLowerHull==1){
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
    triangulate(&c, Pptr);
  }
*/

  for(j=0;j<npoints;j++){
  //  printf("%f %f %f\n", P[j].x[0],P[j].x[1], P[j].x[2]);
  }
  i=0;
  for(f=0;f<c.used;f++){
    c.facet[f].isMFace = checkFacet(&c, f, npoints, Pptr);
    c.facet[f].isLowerHull = checkLowerHull(c.facet[f], npoints, Pptr, f);

    if(c.facet[f].isLowerHull == 1){
      //        printf("%d %d %d is Mface: %d is Lower hull: %d\n", c.facet[f].pointID[0]->id, c.facet[f].pointID[1]->id, c.facet[f].pointID[2]->id, c.facet[f].isMFace, c.facet[f].isLowerHull);
  //    printf("%d %d %d\n", c.facet[f].pointID[0]->id, c.facet[f].pointID[1]->id, c.facet[f].pointID[2]->id);
    i++;
    }
  }
  printf("Number of hyperplanes: %d\nNumber of Facets: %d\n", c.used, i);
  
  freeContainer(&c);
  gsl_rng_free(ran);
}
