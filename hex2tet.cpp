#include "hex2tet.h"

static
int H2T_Add_tetra(MMG5_pMesh mesh,int v0,int v1,int v2,int v3,int ref,int ncut,int nhex) {
  int iel;

  iel =  MMG3D_Add_tetrahedron ( mesh,v0,v1,v2,v3,ref);

  if ( !iel ) {
    H2T_MAXTET_ERROR_MESSAGE(__func__,__LINE__,mesh->nemax,ncut,nhex);
    return 0;
  }
  if ( iel < 0 ) {
    printf("  ## Error:%s: tetra %d %d %d %d: wrong orientation\n",
           __func__,v0,v1,v2,v3);
    return 0;
  }
  //mesh->tetra[iel].mark=0;

  return 1;
}

int H2T_edgePut(pHedge hash,int a,int b,int np) {
  int        key,mins,maxs;
  hedge     *ha;

  mins = a;
  maxs = b;
  if ( a > b ) {
    mins = b;
    maxs = a;
  }
  key = H2T_KA*mins + H2T_KB*maxs;
  key = key % hash->size;
  ha  = &hash->item[key];

  if ( ha->min ) {
    /* Same edge */
    if ( ha->min == mins && ha->max == maxs ) {
      return ha->iel;
    }
    else {
      while ( ha->nxt && ha->nxt < hash->nhmax ) {
        ha = &hash->item[ha->nxt];
        if ( ha->min == mins && ha->max == maxs ) {
          return ha->iel;
        }
      }
    }
    ha->nxt = hash->hnxt;
    ha      = &hash->item[hash->hnxt];
    ++hash->hnxt;
    if ( hash->hnxt >= hash->nhmax ) {
      fprintf(stdout,"  ## Memory alloc problem (edge): %d\n",hash->nhmax);
      return 0;
    }
  }

  /* insert */
  ha->min = mins;
  ha->max = maxs;
  ha->iel = np;
  ha->nxt = 0;

  return 1;
}




double H2T_quickvol(double *c1,double *c2,double *c3,double *c4) {
  double   ax,ay,az,bx,by,bz,vol;

  ax = c3[0] - c1[0];
  ay = c3[1] - c1[1];
  az = c3[2] - c1[2];

  bx = c4[0] - c1[0];
  by = c4[1] - c1[1];
  bz = c4[2] - c1[2];

  vol =   (c2[0]-c1[0]) * (ay*bz - az*by) \
        + (c2[1]-c1[1]) * (az*bx - ax*bz) \
        + (c2[2]-c1[2]) * (ax*by - ay*bx);

  return vol;
}

int H2T_chkorient(MMG5_pMesh mmgMesh,int* hexa,int nhex) {
  int     nbado,k,i,iadr,ph[8];
  double  volref,volhex;

  nbado = 0;
  volref = 1;
  double* c1 = new double[3];
  double* c2 = new double[3];
  double* c3 = new double[3];
  double* c4 = new double[3];
  for (k=1; k<=nhex; k++)
  {
    iadr = 9*k;
    for(i=0 ; i<8 ; i++)
      ph[i] = hexa[iadr+i];
      
//    c1[0] = xcn->getVal(ph[0],0); c1[1] = xcn->getVal(ph[0],1); c1[2] = xcn->getVal(ph[0],2);
//    c2[0] = xcn->getVal(ph[1],0); c2[1] = xcn->getVal(ph[1],1); c2[2] = xcn->getVal(ph[1],2);
//    c3[0] = xcn->getVal(ph[3],0); c3[1] = xcn->getVal(ph[3],1); c3[2] = xcn->getVal(ph[3],2);
//    c4[0] = xcn->getVal(ph[4],0); c4[1] = xcn->getVal(ph[4],1); c4[2] = xcn->getVal(ph[4],2);
      
    c1[0] = mmgMesh->point[ph[0]].c[0]; c1[1] = mmgMesh->point[ph[0]].c[1]; c1[2] = mmgMesh->point[ph[0]].c[2];
      
    c2[0] = mmgMesh->point[ph[1]].c[0]; c2[1] = mmgMesh->point[ph[1]].c[1]; c2[2] = mmgMesh->point[ph[1]].c[2];
      
    c3[0] = mmgMesh->point[ph[3]].c[0]; c3[1] = mmgMesh->point[ph[3]].c[1]; c3[2] = mmgMesh->point[ph[3]].c[2];
        
    c4[0] = mmgMesh->point[ph[4]].c[0]; c4[1] = mmgMesh->point[ph[4]].c[1]; c4[2] = mmgMesh->point[ph[4]].c[2];
      
//    c2[0] = xcn->getVal(ph[1],0); c2[1] = xcn->getVal(ph[1],1); c2[2] = xcn->getVal(ph[1],2);
//    c3[0] = xcn->getVal(ph[3],0); c3[1] = xcn->getVal(ph[3],1); c3[2] = xcn->getVal(ph[3],2);
//    c4[0] = xcn->getVal(ph[4],0); c4[1] = xcn->getVal(ph[4],1); c4[2] = xcn->getVal(ph[4],2);
    /** check the orientability of the hexahedra : vol of tet p0 p1 p3 p4 */
    volhex = H2T_quickvol(c1,c2,c3,c4);

    if ( volref*volhex < 0 )
    {
      nbado++;
      hexa[iadr + 3] = ph[1];
      hexa[iadr + 1] = ph[3];
      hexa[iadr + 5] = ph[7];
      hexa[iadr + 7] = ph[5];
    }
  }
  return nbado;
}

int H2T_chkAdja(MMG5_pMesh mesh,int* listhexa,int* adjahex,int nhex) {
  int *list,*mark;
  int icurc,ipil,iadr,adj,count,k,i;

  /** Allocs */
  list = mark = NULL;
  list = (int*) calloc(7*nhex+1,sizeof(int));
  assert(list);
  mark = (int*) calloc(nhex+1,sizeof(int));
  assert(mark);


  /** Stack initialization: mark the frist hexa as seen and append its
   * adjacents to the stack */
  icurc   = 0;
  ipil    = 0;
  mark[1] = 1;
  iadr    = 1;

  for ( i=0; i<6; i++ ) {
    adj  = adjahex[iadr + i];

    if ( !adj ) continue;

    list[ipil++] = adj;
  }

  /** Process the next hexa of the stack: mark it as seen and append its adjacents */
  while ( icurc++ < ipil ) {
    k = list[icurc-1]/6;

    if ( !k ) continue;

    mark[k] = 1;

    for ( i=0; i<6; i++ ) {
      iadr = 6*(k-1)+1;
      adj = adjahex[iadr + i];

      if ( !adj ) continue;

      /* Add only once each tetra */
      if ( mark[adj/6] )
        continue;
      else {
        mark[adj/6] = -1;
      }

      list[ipil++] = adj;
    }
  }

  /** Count the number of marked tetra (already seen) */
  count = 0;
  for ( i=1; i<=nhex; ++i ) {
    if ( mark[i] > 0 ) ++count;
  }

  /** Free the memory */
  free(list); list = NULL;
  free(mark); mark = NULL;

  return count;
}

static
int H2T_decouphex(MMG5_pMesh mesh, pHedge hed,int* p,int ref,int *ncut,int nhex) {

  /** Creation of the first tetra (0,1,3,7) */
  if ( !H2T_Add_tetra(mesh,p[0],p[1],p[3],p[7],ref,*ncut,nhex) ) return 0;

  /** Creation of the second tetra (7,2,6,1) */
  if ( !H2T_Add_tetra(mesh,p[7],p[2],p[6],p[1],ref,*ncut,nhex) ) return 0;

  /** Creation of the third tetra (1,4,5,7) */
  if ( !H2T_Add_tetra(mesh,p[1],p[4],p[5],p[7],ref,*ncut,nhex) ) return 0;

  /** Creation of the fourth tetra (7,4,0,1) */
  if ( !H2T_Add_tetra(mesh,p[7],p[4],p[0],p[1],ref,*ncut,nhex) ) return 0;

  /** Creation of the fourth tetra (1,6,7,5) */
  if ( !H2T_Add_tetra(mesh,p[1],p[6],p[7],p[5],ref,*ncut,nhex) ) return 0;

  /** Creation of the sixth tetra (1,3,7,2) */
  if ( !H2T_Add_tetra(mesh,p[1],p[3],p[7],p[2],ref,*ncut,nhex) ) return 0;

  /** Add edges to the hashtable */
  H2T_edgePut(hed,p[0],p[7],2);
  H2T_edgePut(hed,p[1],p[3],2);
  H2T_edgePut(hed,p[2],p[7],2);
  H2T_edgePut(hed,p[1],p[6],2);
  H2T_edgePut(hed,p[1],p[4],2);
  H2T_edgePut(hed,p[5],p[7],2);

  ++(*ncut);

  return 1;
}

int H2T_edgePoint(pHedge hash,int a,int b) {
  int        key,mins,maxs;
  hedge     *ha;

  /* compute key */
  mins = a;
  maxs = b;
  if ( a > b ) {
    mins = b;
    maxs = a;
  }
  key = H2T_KA*mins + H2T_KB*maxs;
  key = key % hash->size;
  ha  = &hash->item[key];

  if ( !ha->min )  return 0;
  else if ( ha->min == mins && ha->max == maxs ) {
    return ha->iel;
  }
  else if ( ha->nxt ) {
    do {
      ha = &hash->item[ha->nxt];
      if ( ha->min == mins && ha->max == maxs ) {
        return ha->iel;
      }
    }
    while ( ha->nxt && ha->nxt < hash->nhmax );
  }
  return 0;
}

static int H2T_checkcase(int ph[8],int nu1,int nu2,pHedge hed) {
  int i,nu3;

  for ( i=0; i<3; i++ ) {
    /** Check the existence of edge nu1-nu3 with nu3 the point
     * opposite to nu1 inside the hexa faces */
    nu3 = H2T_hied[nu1][i];

    /** Do not check the edge from which we come */
    if ( nu3==nu2 ) continue;

    /** Test the edge existence */
    if ( H2T_edgePoint(hed,ph[nu1],ph[nu3]) ) break;
  }

  if ( i<3 ) return 1;
  else return 0;
}

/**
 * \param ph hexa vertices
 * \param nu1 First extremity of a non created inside the hexa
 * \param nu2 Second extremity of a non created diagonale inside the hexa
 * \param hed edge hashtable
 *
 * \return 1 if at least 2 edges of the opposite case to the default case exist.
 *
 * Check if we have already hashed at least 1 diagonal issue from nu1
 * (of the hexa faces) without taking into account the nu1-nu2
 * diagonal.
 *
 */
static int H2T_checkcaseopp(int ph[8],int nu1,int nu2,pHedge hed) {
  int i,nu3,nu4;

  /** Point opposite to nu1 in the hexa */
  nu4 = H2T_hop[nu1];

  for ( i=0; i<3; i++ ) {
    /** Check the existence of edge nu4-nu3 with nu3 the point
     * opposite to nu1 inside the hexa faces */
    nu3 = H2T_hied[nu4][i];

    /** Do not check the edge hop[nu1]-hop[nu2] */
    if ( nu3==H2T_hop[nu2] ) continue;

    /** Test the edge existence */
    if ( H2T_edgePoint(hed,ph[nu4],ph[nu3]) ) break;
  }

  if ( i<3 ) return 1;
  else return 0;
}

int H2T_cuthex(MMG5_pMesh mesh,pHedge hed,int* listhexa,int* adjahex,int nhex) {
  MMG5_pPoint    ppt;
  int            i,ih,k,nu1,nu2,nu3,nu4,adj,icas0,icasopp,nncut;
  int            *list,*mark,p[8],ipil,icurc,iface,iadr;
  int            iel,ip,ph[8];
  double         c[3];
  int            ddebug,ncut;

  if ( mesh->info.ddebug ) {
    int count = H2T_chkAdja(mesh,listhexa,adjahex,nhex);
    printf("Number of hexa reached by adjacency: %d/%d\n",count,nhex);
  }
  /* Alloc */
  list = (int*) calloc(7*nhex+1,sizeof(int));
  assert(list);
  mark = (int*) calloc(nhex+1,sizeof(int));
  assert(mark);

  /* Stack initialization */
  mark[1] = -1;

    for ( ih=0; ih<8; ih++ ) p[ih] = listhexa[9+ih];

  if ( mesh->info.ddebug )
    printf(" First hexa %d %d %d %d %d %d %d %d\n",p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]);

  ncut = 0;
  if ( !H2T_decouphex(mesh,hed,p,listhexa[9+8],&ncut,nhex) ) return 0;

  icurc = 0;
  ipil  = 0;
  for ( i=0; i<6; i++ ) {
    iadr = 1;
    adj = adjahex[iadr + i];
    if ( !adj) continue;
    list[ipil++] = adj;
  }

  while(icurc++ < ipil) {
    k = list[icurc-1]/6;

    if ( !k ) continue;

    /** Store the hexa */
    for(ih=0;ih<8;ih++) ph[ih] = listhexa[9*k+ih];
    if ( mesh->info.ddebug ) {
      printf("ph %d %d %d %d %d %d %d %d\n",
             ph[0],ph[1],ph[2],ph[3],ph[4],ph[5],ph[6],ph[7]);
    }

    mark[k] = -1;

    iface = list[icurc-1]%6;
    if ( iface < 0 ) {
      printf(" ## Error: %s: wrong adjacency (%d)\n",__func__,iface);
      return 0;
    }
    ddebug=0;

    nu1 = H2T_hidir[iface][0];
    nu2 = H2T_hidir[iface][2];

    /** Test if the diagonal 0-2 of iface exist */
    if ( H2T_edgePoint(hed,ph[nu1],ph[nu2]) ) {

      /** Test if the opposite diagonal on the opposite face exist */
      nu3 = H2T_hidirop[iface][0];
      nu4 = H2T_hidirop[iface][2];
      if ( H2T_edgePoint(hed,ph[nu3],ph[nu4]) ) {
        /** If yes, mark the hexa */
        mark[k] = -10;
        continue;
      }
      if ( iface==1 || iface==5 ) {
        //find if other edge with ph->v[MMG_hidir[iface][0]], if yes->renum
        icas0 = H2T_checkcase(ph,nu1,nu2,hed);
        icasopp = H2T_checkcaseopp(ph,nu1,nu2,hed);
        if (  icas0 || icasopp ) {
          icas0 = H2T_checkcase(ph,nu2,nu1,hed);
          icasopp = H2T_checkcaseopp(ph,nu2,nu1,hed);
          if ( icas0 || icasopp ) {
            mark[k] = -10;
            continue;
          }
          icas0 = 1;
        }
        if ( icas0 ) {
          //debug check
          for ( i=0; i<3; i++ ) {
            nu3 = H2T_hied[nu2][i];
            if ( nu3==nu1 ) continue;
            if ( H2T_edgePoint(hed,ph[nu2],ph[nu3]) ) break;
          }
          assert(i==3);
          //printf("iface %d: another edge founded---> renum\n",iface);
          if ( iface==1 ) {
            p[0] = ph[3]; p[1] = ph[0]; p[2] = ph[1]; p[3] = ph[2];
            p[4] = ph[7]; p[5] = ph[4]; p[6] = ph[5]; p[7] = ph[6];
          } else {
            p[0] = ph[1]; p[1] = ph[2]; p[2] = ph[3]; p[3] = ph[0];
            p[4] = ph[5]; p[5] = ph[6]; p[6] = ph[7]; p[7] = ph[4];
          }
          H2T_decouphex(mesh,hed,p,listhexa[9*k+8],&ncut,nhex);
        } else {
          H2T_decouphex(mesh,hed,ph,listhexa[9*k+8],&ncut,nhex);
        }
      } else if ( iface==4 ) {
        icas0 = H2T_checkcase(ph,nu2,nu1,hed);
        icasopp = H2T_checkcaseopp(ph,nu2,nu1,hed);
        if ( icas0 || icasopp ) {
          icas0 = H2T_checkcase(ph,nu1,nu2,hed);
          icasopp = H2T_checkcaseopp(ph,nu1,nu2,hed);
          if ( icas0 || icasopp ) {
            mark[k] = -10;
            continue;
          }
          icas0 = 1;
        }
        if ( icas0 ) {
          //check debug
          for ( i=0; i<3; i++ ) {
            nu3 = H2T_hied[nu1][i];
            if ( nu3==nu2 ) continue;
            if ( H2T_edgePoint(hed,ph[nu1],ph[nu3]) ) break;
          }
          assert(i==3);
          //printf("iface 4 another edge founded ---> renum\n");
          p[0] = ph[3]; p[1] = ph[0]; p[2] = ph[1]; p[3] = ph[2];
          p[4] = ph[7]; p[5] = ph[4]; p[6] = ph[5]; p[7] = ph[6];
          H2T_decouphex(mesh,hed,p,listhexa[9*k+8],&ncut,nhex);
        } else {
          H2T_decouphex(mesh,hed,ph,listhexa[9*k+8],&ncut,nhex);
        }
      } else {
        if ( ddebug)  printf("face %d renumbering\n",iface);//iface 0,2,3
        icas0 = H2T_checkcase(ph,nu1,nu2,hed);
        icasopp = H2T_checkcaseopp(ph,nu1,nu2,hed);
        if ( icas0 || icasopp ) {
          icas0 = H2T_checkcase(ph,nu2,nu1,hed);
          icasopp = H2T_checkcaseopp(ph,nu2,nu1,hed);
          if ( icas0 || icasopp ) {
            mark[k] = -10;
            continue;
          }
          icas0 = 1;
        }
        if ( icas0 ) {
          //check debug
          for ( i=0; i<3; i++ ) {
            nu3 = H2T_hied[nu2][i];
            if ( nu3==nu1 ) continue;
            if ( H2T_edgePoint(hed,ph[nu2],ph[nu3]) ) break;
          }
          assert(i==3);
          icas0=1;
        }
        if ( ddebug) printf("icas %d\n",icas0);
        switch(iface) {
        case(0):
          if ( icas0 ) {
            p[0] = ph[3]; p[1] = ph[0]; p[2] = ph[1]; p[3] = ph[2];
            p[4] = ph[7]; p[5] = ph[4]; p[6] = ph[5]; p[7] = ph[6];
          } else {
            p[0] = ph[1]; p[1] = ph[2]; p[2] = ph[3]; p[3] = ph[0];
            p[4] = ph[5]; p[5] = ph[6]; p[6] = ph[7]; p[7] = ph[4];
          }
          break;
        case(2):
          if ( icas0 ) {
            p[0] = ph[3]; p[1] = ph[0]; p[2] = ph[1]; p[3] = ph[2];
            p[4] = ph[7]; p[5] = ph[4]; p[6] = ph[5]; p[7] = ph[6];
          } else {
            p[0] = ph[2]; p[1] = ph[3]; p[2] = ph[0]; p[3] = ph[1];
            p[4] = ph[6]; p[5] = ph[7]; p[6] = ph[4]; p[7] = ph[5];
          }
          break;
        case(3):
          if ( icas0 ) {
            p[0] = ph[1]; p[1] = ph[2]; p[2] = ph[3]; p[3] = ph[0];
            p[4] = ph[5]; p[5] = ph[6]; p[6] = ph[7]; p[7] = ph[4];
          } else {
            p[0] = ph[3]; p[1] = ph[0]; p[2] = ph[1]; p[3] = ph[2];
            p[4] = ph[7]; p[5] = ph[4]; p[6] = ph[5]; p[7] = ph[6];
          }
          break;
        }

        H2T_decouphex(mesh,hed,p,listhexa[9*k+8],&ncut,nhex);
      }
    }  else if ( H2T_edgePoint(hed,ph[H2T_hidir[iface][1]],ph[H2T_hidir[iface][3]]) ) {
      /** The edge 1-3 (opposite to 0-2) of iface exist */
      nu1 = H2T_hidir[iface][1];
      nu2 = H2T_hidir[iface][3];

      nu3 = H2T_hidirop[iface][1];
      nu4 = H2T_hidirop[iface][3];
      /** Test if the opposite edge on the opposite face exist */
      if ( H2T_edgePoint(hed,ph[nu3],ph[nu4]) ) {
        /** nothing to do */
        mark[k] = -10;
        continue;
      }
      if ( iface==0 || iface==3 ) {
        icas0 = H2T_checkcase(ph,nu1,nu2,hed);
        icasopp = H2T_checkcaseopp(ph,nu1,nu2,hed);
        if ( icas0 || icasopp ) {
          icas0 = H2T_checkcase(ph,nu2,nu1,hed);
          icasopp = H2T_checkcaseopp(ph,nu2,nu1,hed);
          if ( icas0 || icasopp ) {
            mark[k] = -10;
            continue;
          }
          icas0 = 1;
        }
        if ( ddebug) printf("icas0 %d\n",icas0);
        if ( icas0 ) {
          //check debug
          for ( i=0; i<3; i++ ) {
            nu3 = H2T_hied[nu2][i];
            if ( nu3==nu1 ) continue;
            if ( H2T_edgePoint(hed,ph[nu2],ph[nu3]) ) {
              break;
            }
          }
          assert(i==3);
          if ( iface==0 ) {
            p[0] = ph[2]; p[1] = ph[3]; p[2] = ph[0]; p[3] = ph[1];
            p[4] = ph[6]; p[5] = ph[7]; p[6] = ph[4]; p[7] = ph[5];
          } else {
            p[0] = ph[2]; p[1] = ph[3]; p[2] = ph[0]; p[3] = ph[1];
            p[4] = ph[6]; p[5] = ph[7]; p[6] = ph[4]; p[7] = ph[5];
          }
          H2T_decouphex(mesh,hed,p,listhexa[9*k+8],&ncut,nhex);
        } else {
          H2T_decouphex(mesh,hed,ph,listhexa[9*k+8],&ncut,nhex);
        }
      } else if ( iface==2 ) {
        icas0 = H2T_checkcase(ph,nu2,nu1,hed);
        icasopp = H2T_checkcaseopp(ph,nu2,nu1,hed);
        if ( icas0 || icasopp ) {
          icas0 = H2T_checkcase(ph,nu1,nu2,hed);
          icasopp = H2T_checkcaseopp(ph,nu1,nu2,hed);
          if ( icas0 || icasopp ) {
            mark[k] = -10;
            continue;
          }
          icas0 = 1;
        }
        if ( icas0 ) {
          //check debug
          for ( i=0; i<3; i++ ) {
            nu3 = H2T_hied[nu1][i];
            if ( nu3==nu2 ) continue;
            if ( H2T_edgePoint(hed,ph[nu1],ph[nu3]) ) break;
          }
          assert(i==3);
          p[0] = ph[1]; p[1] = ph[2]; p[2] = ph[3]; p[3] = ph[0];
          p[4] = ph[5]; p[5] = ph[6]; p[6] = ph[7]; p[7] = ph[4];
          H2T_decouphex(mesh,hed,p,listhexa[9*k+8],&ncut,nhex);
        } else {
          H2T_decouphex(mesh,hed,ph,listhexa[9*k+8],&ncut,nhex);
        }
      }
      else {
        if ( ddebug)  printf("face %d renumbering\n",iface);//iface 1,4,5
        icas0 = H2T_checkcase(ph,nu1,nu2,hed);
        icasopp = H2T_checkcaseopp(ph,nu1,nu2,hed);
        if ( icas0 || icasopp ) {
          icas0 = H2T_checkcase(ph,nu2,nu1,hed);
          icasopp = H2T_checkcaseopp(ph,nu2,nu1,hed);
          if ( icas0 || icasopp ) {
            mark[k] = -10;
            continue;
          }
          icas0 = 1;
        }
        switch(iface) {
        case(1):
          if ( icas0 ) {
            p[0] = ph[1]; p[1] = ph[2]; p[2] = ph[3]; p[3] = ph[0];
            p[4] = ph[5]; p[5] = ph[6]; p[6] = ph[7]; p[7] = ph[4];
          } else {
            p[0] = ph[2]; p[1] = ph[3]; p[2] = ph[0]; p[3] = ph[1];
            p[4] = ph[6]; p[5] = ph[7]; p[6] = ph[4]; p[7] = ph[5];
          }
          break;
        case(4):
          if ( ddebug ) {
            printf("at the beginning %d : %d %d %d %d %d %d %d %d\n",
                   k,ph[0],ph[1],ph[2],ph[3],ph[4],ph[5],ph[6],ph[7]);
          }

          if ( icas0 ) {
            p[0] = ph[1]; p[1] = ph[2]; p[2] = ph[3]; p[3] = ph[0];
            p[4] = ph[5]; p[5] = ph[6]; p[6] = ph[7]; p[7] = ph[4];
          } else {
            p[0] = ph[2]; p[1] = ph[3]; p[2] = ph[0]; p[3] = ph[1];
            p[4] = ph[6]; p[5] = ph[7]; p[6] = ph[4]; p[7] = ph[5];
          }
          if ( ddebug ) {
            printf("at the end %d : %d %d %d %d %d %d %d %d\n",
                   k,p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]);
          }
          break;
        case(5):
          if ( icas0 ) {
            p[0] = ph[2]; p[1] = ph[3]; p[2] = ph[0]; p[3] = ph[1];
            p[4] = ph[6]; p[5] = ph[7]; p[6] = ph[4]; p[7] = ph[5];
          } else {
            p[0] = ph[3]; p[1] = ph[0]; p[2] = ph[1]; p[3] = ph[2];
            p[4] = ph[7]; p[5] = ph[4]; p[6] = ph[5]; p[7] = ph[6];
          }
          break;
        }
        H2T_decouphex(mesh,hed,p,listhexa[9*k+8],&ncut,nhex);
      }
    } else {
      /** We have no edge on iface: not normal because the adjacent
       * through iface has been cutted */
      printf("Error: %s: no created diagonal on face %d of hexa: %d %d and %d %d tested\n",
             __func__,iface,ph[H2T_hidir[iface][1]],ph[H2T_hidir[iface][3]],
             ph[H2T_hidir[iface][0]],ph[H2T_hidir[iface][2]]);
      printf("hexa: %d %d %d %d %d %d %d %d\n",ph[0],ph[1],ph[2],ph[3],ph[4],ph[5],ph[6],ph[7]);
      return 0;
    }

    /** Append the adjacents of the current tetra to the stack */
    for ( i=0; i<6; i++ ) {
      iadr = 6*(k-1)+1;
      adj = adjahex[iadr + i];

      if ( !adj ) continue;

      if ( mark[adj/6] )
        continue;
      else
        mark[adj/6] = 10;

      list[ipil++] = adj;

      if ( mesh->info.ddebug )
        printf("k=%d: stack append: hexa %d (iface %d) in %d -- through face %d\n",
               k,adj/6,adj%6,ipil-1,i);
    }
  }

  if ( ncut < nhex ) {
    printf("     %d/%d Hexa succesfully cutted\n",ncut,nhex);
    printf("--> Try to cut the remaining hexa using point insertion\n");
  }

  /** Treat the remaining hexa */
  nncut = 0;
  for ( k=1; k<=nhex; k++ ) {
    if ( mark[k]==-1 ) continue;
    for(i=0;i<8;i++) ph[i] = listhexa[9*k+i];
    nncut++;
    /** create new vertex */
    c[0] = c[1] = c[2] = 0.;
    for ( i=0; i<8; i++ ) {
      ppt = &mesh->point[ph[i]];
      c[0] += ppt->c[0];        c[1] += ppt->c[1]; c[2] += ppt->c[2];
    }
    c[0] /= 8.; c[1] /= 8.; c[2] /= 8.;

    ip = MMG3D_Add_vertex(mesh,c[0],c[1],c[2],0);
    if ( !ip ) return 0;

    /** create 2 tets per faces */
    int ref = listhexa[9*k+8];

    for ( i=0 ;i<6; i++ ) {
      nu1 = H2T_hidir[i][0];
      nu2 = H2T_hidir[i][2];
      if ( H2T_edgePoint(hed,ph[nu1],ph[nu2]) ) {
        if ( mesh->ne+1 >= mesh->nemax ) {
          H2T_MAXTET_ERROR_MESSAGE(__func__,__LINE__,mesh->nemax,ncut+nncut,nhex);
          fprintf(stdout,"%d new points.\n",nncut);
          return 0;
        }

        /* Tetra ip,ph[nu1],ph[nu2],ph[H2T_hidir[i][1]] */
        iel =  MMG3D_Add_tetrahedron ( mesh,ip,ph[nu1],ph[nu2],ph[H2T_hidir[i][1]],ref );
        if ( !iel ) {
          H2T_MAXTET_ERROR_MESSAGE(__func__,__LINE__,mesh->nemax,ncut+nncut,nhex);
          fprintf(stdout,"%d new points.\n",nncut);
          return 0;
        }

        /* Tetra ip,ph[nu1],ph[H2T_hidir[i][3]],ph[nu2] */
        iel =  MMG3D_Add_tetrahedron ( mesh,ip,ph[nu1],ph[H2T_hidir[i][3]],ph[nu2],ref );
        if ( !iel ) {
          H2T_MAXTET_ERROR_MESSAGE(__func__,__LINE__,mesh->nemax,ncut+nncut,nhex);
          fprintf(stdout,"%d new points.\n",nncut);
          return 0;
        }
      } else {
        nu1 = H2T_hidir[i][1];
        nu2 = H2T_hidir[i][3];
        if ( !H2T_edgePoint(hed,ph[nu1],ph[nu2])) H2T_edgePut(hed,ph[nu1],ph[nu2],2);
        if ( mesh->ne+1 >= mesh->nemax ) {
          H2T_MAXTET_ERROR_MESSAGE(__func__,__LINE__,mesh->nemax,ncut+nncut,nhex);
          fprintf(stdout,"%d new points.\n",nncut);
          return 0;
        }

        /* Tetra ip,ph[nu1],ph[H2T_hidir[i][0]],ph[nu2] */
        iel =  MMG3D_Add_tetrahedron ( mesh,ip,ph[nu1],ph[H2T_hidir[i][0]],ph[nu2],ref );
        if ( !iel ) {
          H2T_MAXTET_ERROR_MESSAGE(__func__,__LINE__,mesh->nemax,ncut+nncut,nhex);
          fprintf(stdout,"%d new points.\n",nncut);
          return 0;
        }

        /* Tetra ip,ph[nu1],ph[nu2],ph[H2T_hidir[i][2]] */
        iel =  MMG3D_Add_tetrahedron ( mesh,ip,ph[nu1],ph[nu2],ph[H2T_hidir[i][2]],ref );
        if ( !iel ) {
          H2T_MAXTET_ERROR_MESSAGE(__func__,__LINE__,mesh->nemax,ncut+nncut,nhex);
          fprintf(stdout,"%d new points.\n",nncut);
          return 0;
        }
      }
    }
  }
  if ( nncut) fprintf(stdout,"  $$ %8d ADDED VERTEX\n",nncut);
  return 1;
}



int H2T_hashHexa(int* listhexa,int* adjahex, int nhex)
{
  int              k,kk,pp,l,ll,mins,mins1,opps,opps1,sum,sum1,iadr;
  int              *hcode,*link,hsize,imins,imins1;
  int              ph[8],ph1[8];
  long int         inival;
  unsigned char    i,ii,iii,i1,i2,i3,i4;
  unsigned int     key;

  /* default */
  fprintf(stdout,"  ** SETTING HEXA ADJACENCIES\n\n");
  fflush(stdout);

  /* memory alloc */
  hcode = (int*)calloc(nhex+1,sizeof(int));
  assert ( hcode );
  link  = adjahex;
  hsize = nhex;

  /* init */
  inival = INT_MAX;
  for (k=0; k<=nhex; k++)
    hcode[k] = -inival;

  /* build hash table */
  for ( k=1; k<=nhex; k++ ) {
    iadr = 9*k;
    for ( iii=0; iii<8; iii++ )
      ph[iii] = listhexa[iadr+iii];

    for ( i=0; i<6; i++ ) {
      i1 = H2T_hidir[i][0];
      i2 = H2T_hidir[i][1];
      i3 = H2T_hidir[i][2];
      i4 = H2T_hidir[i][3];

      H2T_MINVAL_AND_LOC(ph[i1],ph[i2],ph[i3],ph[i4],mins,imins);
      opps = ph[H2T_hopp[i][imins]];

      /* compute key */
      sum = ph[i1] + ph[i2] + ph[i3] + ph[i4];
      key = H2T_KA*mins + H2T_KB*opps + H2T_KC*sum;
      key = key % hsize + 1;

      /* insert */
      iadr = 6*(k-1) + i+1;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
  }

  /* set adjacency */
  for (l=6*nhex; l>0; l--) {

    if ( link[l] >= 0 )  continue;

    /* current element */
    k = ((l-1) /6) + 1;
    i = (l-1) % 6;
    i1 = H2T_hidir[i][0];
    i2 = H2T_hidir[i][1];
    i3 = H2T_hidir[i][2];
    i4 = H2T_hidir[i][3];
    iadr = 9*k;
    for(iii=0 ; iii<8 ; iii++)
      ph[iii] = listhexa[iadr+iii];
    sum  = ph[i1] + ph[i2] + ph[i3] + ph[i4];

    H2T_MINVAL_AND_LOC(ph[i1],ph[i2],ph[i3],ph[i4],mins,imins);
    opps = ph[H2T_hopp[i][imins]];

    /* accross link */
    ll = -link[l];
    pp = 0;
    link[l] = 0;
    while ( ll != inival ) {
      kk = ((ll-1) /6) + 1;
      ii = (ll-1) % 6;
      i1 = H2T_hidir[ii][0];
      i2 = H2T_hidir[ii][1];
      i3 = H2T_hidir[ii][2];
      i4 = H2T_hidir[ii][3];
      iadr = 9*kk;
      for(iii=0 ; iii<8 ; iii++)
        ph1[iii] = listhexa[iadr+iii];
      sum1 = ph1[i1] + ph1[i2] + ph1[i3] + ph1[i4];
      if ( sum1 == sum ) {

        H2T_MINVAL_AND_LOC(ph1[i1],ph1[i2],ph1[i3],ph1[i4],mins1,imins1);
        opps1 = ph1[H2T_hopp[ii][imins1]];

        if ( mins1 == mins ) {
          if ( opps1 == opps ) {
            /* adjacent found */
            if ( pp != 0 )  link[pp] = link[ll];
            link[l] = 6*kk + ii;
            link[ll]= 6*k + i;
            break;
          }
        }
      }
      pp = ll;
      ll = -link[ll];
    }
  }

  free(hcode); hcode=NULL;

  return 1;
}
