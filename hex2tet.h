#include <vector>
#include <math.h>
#include <unordered_set>
#include <set>
#include <map>
#include <string>
#include <utility>
#include <math.h>
#include <ctime>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string.h>
#include <assert.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>

#include <sstream>
#include <iomanip>
#include "adapt_io.h"
#include "adapt_compute.h"
#include "adapt_operations.h"
#include "adapt_math.h"
#include "adapt_output.h"
#include "adapt_partition.h"
#include <sstream>
#include <iomanip>
#include "mmg/mmgs/libmmgs.h"
#include "mmg/mmg3d/libmmg3d.h"
#include "parmmg/libparmmg.h"

#define H2T_STRONGFAILURE 2

#define H2T_MAXTET_ERROR_MESSAGE(func,line,nemax,ncut,nhex) do          \
{                                                                     \
  fprintf(stdout,"%s:%d: max number of tet reached (%d). %d/%d hexa treated.\n", \
          (func),(line),(nemax),(ncut),(nhex));                       \
} while(0)

#define H2T_MINVAL_AND_LOC(v0,v1,v2,v3,minval,minloc) do \
{                                                      \
                                                       \
  minval = v0; minloc = 0;                             \
  if ( v1 < minval ) {                                 \
    minval = v1; minloc = 1;                           \
  }                                                    \
  if ( v2 < minval ) {                                 \
    minval = v2; minloc = 2;                           \
  }                                                    \
  if ( v3 < minval ) {                                 \
    minval = v3; minloc = 3;                           \
  }                                                    \
} while (0)

/** Constant for hash table key computation */
#define H2T_KA 31
/** Constant for hash table key computation */
#define H2T_KB 57
/** Constant for hash table key computation */
#define H2T_KC 79

typedef struct {
  int      min,max,iel,nxt;
} hedge;

typedef struct {
  int      size,nhmax,hnxt;
  hedge   *item;
} Hedge;
typedef Hedge * pHedge;

/** \brief hidir[i]: vertices of the face i */
static const unsigned char H2T_hidir[6][4] = { {0,3,2,1},{0,4,7,3},{0,1,5,4},{4,5,6,7},{1,2,6,5},{2,3,7,6} };
/** \brief hidiropp[i]: vertices of the face opposite to the face i */
static const unsigned char  H2T_hidirop[6][4] = { {7,4,5,6}, {5,1,2,6}, {7,3,2,6}, {1,0,3,2}, {3,0,4,7}, {0,1,5,4} };
/** \brief hopp[i][j]: vertex of face i opposite to the vertex j  */
static const unsigned char H2T_hopp [6][4] = { {2,1,0,3},{7,3,0,4},{5,4,0,1},{6,7,4,5},{6,5,1,2},{7,6,2,3} };
/** \brief hied[i][.]: diagonal starting from the vertex i of one of the tetra face */
static const unsigned char H2T_hied[8][3] = { {2,5,7}, {3,4,6}, {0,5,7}, {1,4,6}, {1,3,6}, {0,2,7}, {1,3,4}, {0,2,5} };
/** \brief hop[i]: vertex opposite to the vertex i */
static const unsigned char H2T_hop[8] = { 6,7,4,5,2,3,0,1 };

static
int H2T_Add_tetra(MMG5_pMesh mesh,int v0,int v1,int v2,int v3,int ref,int ncut,int nhex);

int H2T_edgePut(pHedge hash,int a,int b,int np);

double H2T_quickvol(double *c1,double *c2,double *c3,double *c4);

int H2T_chkorient(MMG5_pMesh mmgMesh,int* hexa,int nhex);

int H2T_chkAdja(MMG5_pMesh mesh,int* listhexa,int* adjahex,int nhex);

static
int H2T_decouphex(MMG5_pMesh mesh, pHedge hed,int* p,int ref,int *ncut,int nhex);

int H2T_edgePoint(pHedge hash,int a,int b);

static int H2T_checkcase(int ph[8],int nu1,int nu2,pHedge hed);

static int H2T_checkcaseopp(int ph[8],int nu1,int nu2,pHedge hed);

int H2T_cuthex(MMG5_pMesh mesh,pHedge hed,int* listhexa,int* adjahex,int nhex);

int H2T_hashHexa(int* listhexa,int* adjahex, int nhex);
