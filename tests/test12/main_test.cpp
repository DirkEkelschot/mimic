#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include "../../src/adapt_parops.h"
#include "../../src/adapt_output.h"
#include "../../src/adapt_boundary.h"
#include "../../src/adapt_distri_parstate.h"
#include <iomanip>

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))

// This is basically textbook recursive merge sort using std::merge_inplace
// but it considers the offsets of segments that are already sorted
void merge_indexed(int data[], const int offsets[], size_t index_begin, size_t index_end)
{
    if (index_end - index_begin > 1) {
        auto index_middle = index_begin + (index_end - index_begin) / 2;
        merge_indexed(data, offsets, index_begin, index_middle);
        merge_indexed(data, offsets, index_middle, index_end);
        std::inplace_merge(&data[offsets[index_begin]], &data[offsets[index_middle]], &data[offsets[index_end]]);
    }
}




int main(int argc, char** argv)
{
    MPI_Init(NULL, NULL);
    FILE            *inm;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j,k;
    
    int ier,opt;
    
    const char* fn_grid="../test_mesh/cylinder_hybrid/grid.h5";
    const char* fn_conn="../test_mesh/cylinder_hybrid/conn.h5";
    const char* fn_data="../test_mesh/cylinder_hybrid/data.h5";
    const char* fn_metric = "metric.inp";
    
    std::vector<double> metric_inputs = ReadMetricInputs(fn_metric);
    
    int ReadFromStats = 0;
    if(metric_inputs.size()==6)
    {
        ReadFromStats=metric_inputs[5];
    }
    
    US3D* us3d    = ReadUS3DData(fn_conn,fn_grid,fn_data,ReadFromStats,comm,info);
    int Nve = us3d->xcn->getNglob();
    
    int Nel_part = us3d->ien->getNrow();
    
    Array<double>* xcn_ref = ReadDataSetFromFile<double>(fn_grid,"xcn");
    Array<int>* ien_ref = ReadDataSetFromFile<int>(fn_conn,"ien");

    Array<double>* Ui = new Array<double>(Nel_part,1);
    int varia = 4;
    double rhoState,uState,vState,wState,TState,VtotState,aState,MState;
    for(int i=0;i<Nel_part;i++)
    {
        rhoState = us3d->interior->getVal(i,0);
        uState   = us3d->interior->getVal(i,1);
        vState   = us3d->interior->getVal(i,2);
        wState   = us3d->interior->getVal(i,3);
        TState   = us3d->interior->getVal(i,4);
        VtotState = sqrt(uState*uState+vState*vState+wState*wState);
        aState   = sqrt(1.4*287.05*TState);
        MState = VtotState/aState;
        Ui->setVal(i,0,MState);
    }
    
    delete us3d->interior;

    
    ParallelState* ien_pstate               = new ParallelState(us3d->ien->getNglob(),comm);
    ParallelState* ife_pstate               = new ParallelState(us3d->ifn->getNglob(),comm);
    ParallelState_Parmetis* parmetis_pstate = new ParallelState_Parmetis(us3d->ien,us3d->elTypes,us3d->ie_Nv,comm);
    ParallelState* xcn_pstate               = new ParallelState(us3d->xcn->getNglob(),comm);
    
    clock_t t;
    double tn = 0.0;
    t = clock();
      
    Partition* P = new Partition(us3d->ien, us3d->iee, us3d->ief, us3d->ie_Nv , us3d->ie_Nf,
                                 us3d->ifn, us3d->ife, us3d->if_ref, us3d->if_Nv,
                                 parmetis_pstate, ien_pstate, ife_pstate,
                                 us3d->xcn, xcn_pstate, Ui, us3d->ie_tetCnt, comm);

    
    //std::cout << "understand bitshift " << world_rank << "  " << (world_rank & (~(1 << 3)))<< std::endl;
    //std::vector<int> LocElem     = P->getLocElem();

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
     Domain* pDom = P->getPartitionDomain();
    std::map<int,std::vector<int> > tetras     = pDom->GTetras;
    std::map<int,std::vector<int> > Gtetras     = pDom->GTetras;

    int nTetras = tetras.size();
    int* Elplease = new int[world_size];
    int* red_Elplease = new int[world_size];
    std::vector<int> rankNoEls;
    for(int u=0;u<world_size;u++)
    {
        Elplease[u]     = 0;
        red_Elplease[u] = 0;
    }
    if(tetras.size()==0)
    {
        Elplease[world_rank] =  1;
    }
    else
    {
        Elplease[world_rank] =  0;
    }
    
    MPI_Allreduce(Elplease, red_Elplease, world_size, MPI_INT, MPI_SUM, comm);
    
    
    int NRankNoTets = 0;
    int rankalloc = 0;
    for(int u=0;u<world_size;u++)
    {
        NRankNoTets = NRankNoTets + red_Elplease[u];
    }
    
    
    std::vector<int> loc_part_verts = pDom->loc_part_verts;
    std::vector<Vert*> Verts        = P->getLocalVerts();
    std::map<int,int> lpartv2gv     = pDom->lpartv2gv;
    
    
    /*
    std::string filename = "checkPart_" + std::to_string(world_rank) + ".dat";
    std::ofstream myfile;
    myfile.open(filename);
    myfile << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    myfile <<"ZONE N = " << loc_part_verts.size() << ", E = " << tetras.size() << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;

    for(int i=0;i<loc_part_verts.size();i++)
    {
        int loc_vid  = loc_part_verts[i];
        int glob_vid = lpartv2gv[loc_vid];
        myfile << Verts[loc_vid]->x << " " << Verts[loc_vid]->y << " " << Verts[loc_vid]->z << std::endl;
    }
    int gv0,gv1,gv2,gv3,gv4,gv5,gv6,gv7;
    int lv0,lv1,lv2,lv3,lv4,lv5,lv6,lv7;
    std::map<int,std::vector<int> >::iterator itertet;
    for(itertet=tetras.begin();itertet!=tetras.end();itertet++)
    {
        int glob_id = itertet->first;
        myfile <<   itertet->second[0]+1 << "  " <<
                    itertet->second[1]+1 << "  " <<
                    itertet->second[2]+1 << "  " <<
                    itertet->second[3]+1 << "  " << std::endl;
    }


    myfile.close();
     */
     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //====================================================================================
  
    //std::vector<int> tris_ref   = pDom->faces_ref;

//    int nTet= 0;
//    std::vector<int> tetcnt = P->getTetCnt();
//    std::map<int,int> mapTets;
//    for(int i=0;i<LocElem.size();i++)
//    {
//        if(tetcnt[i]!=-1)
//        {
//            int gEid = LocElem[i];
//            mapTets[gEid] = tetcnt[i];
//            nTet++;
//        }
//     }
     
    //std::cout << "after partitioning rank = " << world_rank << " #tets = " << tetras.size() << " #prisms " << prisms.size() << std::endl;
    
    i_part_map* if_Nv_part_map                 = P->getIF_Nvpartmap();
    i_part_map* ifn_part_map                   = P->getIFNpartmap();
    i_part_map* ife_part_map                   = P->getIFEpartmap();
    i_part_map* ief_part_map                   = P->getIEFpartmap();
    i_part_map* ien_part_map                   = P->getIENpartmap();
    i_part_map* if_ref_part_map                = P->getIFREFpartmap();
    Array<int>* part_global                    = P->getGlobalPartition();
    std::set<int> ushell                       = pDom->ushell;
    // map local tets to local tets for parmmg;
    
    int gvid,gfid;
    
    std::set<int> uverts;
    std::vector<int> uverts_vec;
    std::set<int> ufaces;
    std::vector<int> ufaces_vec;
    
    std::map<int,std::vector<int> >::iterator ite;
    int r0,r1,el0,el1,pos,ra;
    std::map<int,int> gv2lpv;
    std::map<int,int> lv2gpv;
    int ref,nbcfaces;
    int lcv = 0;
    std::map<int,std::vector<int> > ref2bcface;
    std::vector<int> loc_tet_verts;
    std::map<int,std::vector<int> > face_color_map;
    std::set<int> ordered_rank;

    int Nvface;
    std::set<int> utets;
    std::vector<int> utetv;
    
    int telref = 0;
    int duppe1 = 0;
    int duppe2 = 0;
    int duppe3 = 0;
    
    std::set<int> ufaces_comp;
    std::vector<int> uInterFaces;
    std::vector<int> ufaceOnRank;
    std::vector<int> ufaceOnRankRef;
    std::vector<int> ufaceOnRankLoc;

    std::map<int,std::vector<int> > recv_SharedFaces;
    std::map<int,std::vector<int> > send_SharedFaces;
    
    std::map<int,std::vector<int> > recv_SharedFacesLoc;
    std::map<int,std::vector<int> > send_SharedFacesLoc;
    int shf = 0;
    int lf  = 0;
    
    std::map<int,int> sharedFaces;
    std::map<int,int> sharedFaces_notOwned;

    int locOwned  = 0;
    int NewlyRecv = 0;
    
//    int* red_ielement_nlocs   = new int[world_size];
    std::vector<int> red_ielement_nlocs(world_size);
    int* ielement_offsets     = new int[world_size];
    int* ielement_nlocs     = new int[world_size];

    std::map<int,int> gEl2LtetEl;
    int elTel = 0;
    std::map<int,int> lV2gV_tets;
    std::map<int,int> gV2lV_tets;
    std::set<int> gvid_set;
    int lvid = 0;
    
    std::map<int,int> lF2gF_tets;
    std::map<int,int> gF2lF_tets;
    std::set<int> gfid_set;
    int lfid = 0;
    
    int shfn = 0;
    for(ite=tetras.begin();ite!=tetras.end();ite++)
    {
        //key = GID, value global node nmber;
        int gEl             = ite->first;
        gEl2LtetEl[elTel]   = gEl;
        
        for(int j=0;j<4;j++)
        {
            gfid = ief_part_map->i_map[gEl][j];
            if(gfid_set.find(gfid)==gfid_set.end())
            {
                gfid_set.insert(gfid);
                gF2lF_tets[gvid] = lfid;
                lF2gF_tets[lvid] = gfid;
                lfid++;
            }
            
            for(int k=0;k<3;k++)
            {
                gvid = ifn_part_map->i_map[gfid][k];
                
                if(gvid_set.find(gvid)==gvid_set.end())
                {
                    gvid_set.insert(gvid);
                    gV2lV_tets[gvid] = lvid;
                    lV2gV_tets[lvid] = gvid;
                    lvid++;
                }
            }
            
            
            if(ushell.find(gfid)!=ushell.end())
            {
                ref = 13;
            }
            else
            {
                ref = if_ref_part_map->i_map[gfid][0];
            }

            if(ufaces.find(gfid)==ufaces.end())
            {
                ufaces.insert(gfid);
                ufaces_vec.push_back(gfid);
                Nvface = if_Nv_part_map->i_map[gfid][0];

                el0    = ife_part_map->i_map[gfid][0];
                el1    = ife_part_map->i_map[gfid][1];

                if(ref==2)
                {
                    r0 = part_global->getVal(el0,0);
                    r1 = part_global->getVal(el1,0);
                    
                    if(r0==world_rank && r1==world_rank)
                    {
                        ufaceOnRank.push_back(gfid);
                        ufaceOnRankRef.push_back(ref);
                        ufaceOnRankLoc.push_back(lf);
                        lf++;
                    }
                    if(r0==world_rank && r1!=world_rank)
                    {
                        sharedFaces[gfid] = r0;
                        shf++;
                    }
                    if(r0!=world_rank && r1==world_rank)
                    {
                        sharedFaces[gfid] = r1;
                        shf++;
                    }
                }
                
                if(ref!=2)
                {
                    ref2bcface[ref].push_back(gfid);
                    ufaceOnRank.push_back(gfid);
                    ufaceOnRankRef.push_back(ref);
                    nbcfaces++;
                    ufaceOnRankLoc.push_back(lf);
                    lf++;
                }
            }
        }
        
        elTel++;
    }
    
    int nSharedFaces   = sharedFaces.size();
    int nInteriorFaces = ufaceOnRank.size();
    int nLocalVerts    = gV2lV_tets.size();
    int nLocalFaces    = gF2lF_tets.size();

    DistributedParallelState* distSharedFaces = new DistributedParallelState(nSharedFaces,comm);
    DistributedParallelState* distLocalVerts  = new DistributedParallelState(nLocalVerts,comm);
    DistributedParallelState* distLocalFaces  = new DistributedParallelState(nLocalFaces,comm);

    int Nt_shFaces              = distSharedFaces->getNel();
    int* shFace_offsets         = distSharedFaces->getOffsets();
    int* shFace_nlocs           = distSharedFaces->getNlocs();
    int* shFaces_key            = new int[nSharedFaces];
    int* shFaces_val            = new int[nSharedFaces];
    int* TotalSharedFaces       = new int[Nt_shFaces];
    int* TotalSharedFacesRank   = new int[Nt_shFaces];
    
    DistributedParallelState* distInteriorFaces = new DistributedParallelState(nInteriorFaces,comm);
    int Nt_interiorFaces        = distInteriorFaces->getNel();
    int* intFace_offsets        = distInteriorFaces->getOffsets();
    int* intFace_nlocs          = distInteriorFaces->getNlocs();
    int* intFaces_key           = new int[nInteriorFaces];
    int* intFaces_val           = new int[nInteriorFaces];
    int* TotalInteriorFaces_key = new int[Nt_interiorFaces];
    int* TotalInteriorFaces_val = new int[Nt_interiorFaces];

    int iter         = 0;
    std::set<int> sharedVerts_set;
    std::vector<int> sharedVerts;
    std::vector<int> sharedVertsRank;
    
    std::map<int,int>::iterator itsf;
    int lvrtid = 0;
    int tel = shFace_offsets[world_rank];
    for(itsf=sharedFaces.begin();itsf!=sharedFaces.end();itsf++)
    {
        shFaces_key[iter] = itsf->first;
        shFaces_val[iter] = itsf->second;
        gfid = itsf->first;

        for(int q=0;q<3;q++)
        {
            gvid   = ifn_part_map->i_map[gfid][q];
            
            if(sharedVerts_set.find(gvid)==sharedVerts_set.end())
            {
                sharedVerts_set.insert(gvid);
                sharedVerts.push_back(gvid);
                sharedVertsRank.push_back(world_rank);
                lvrtid++;
            }
        }
        tel++;
        iter++;
    }

    int nSharedVerts = sharedVerts.size();

    DistributedParallelState* distSharedVerts = new DistributedParallelState(nSharedVerts,comm);
    
    int Nt_shVerts              = distSharedVerts->getNel();
    int* shVerts_nlocs          = distSharedVerts->getNlocs();
    int* shVerts_offsets        = distSharedVerts->getOffsets();
    int* TotalSharedVerts       = new int[Nt_shVerts];
    int* TotalSharedVertsRank   = new int[Nt_shVerts];

    // Communicate vert map to all ranks.
    MPI_Allgatherv(&sharedVerts[0],
                   nSharedVerts,
                   MPI_INT,
                   TotalSharedVerts,
                   shVerts_nlocs,
                   shVerts_offsets,
                   MPI_INT, comm);
    
    MPI_Allgatherv(&sharedVertsRank[0],
                   nSharedVerts,
                   MPI_INT,
                   TotalSharedVertsRank,
                   shVerts_nlocs,
                   shVerts_offsets,
                   MPI_INT, comm);
    
    // Communicate face map to all ranks.
    MPI_Allgatherv(shFaces_key,
                   nSharedFaces,
                   MPI_INT,
                   TotalSharedFaces,
                   shFace_nlocs,
                   shFace_offsets,
                   MPI_INT, comm);
    
    MPI_Allgatherv(shFaces_val,
                   nSharedFaces,
                   MPI_INT,
                   TotalSharedFacesRank,
                   shFace_nlocs,
                   shFace_offsets,
                   MPI_INT, comm);
    
    int* interiorfaces_key = new int[nInteriorFaces];
    int* interiorfaces_val = new int[nInteriorFaces];
    
    for(int i=0;i<nInteriorFaces;i++)
    {
        interiorfaces_key[i] = intFace_offsets[world_rank]+i;
        interiorfaces_val[i] = interiorfaces_val[i];
    }
    
    int tmp;
    std::map<int,int> f2r;
    std::set<int> f2r_s;
    int* owned_faces = new int[world_size];
    for(int u=0;u<world_size;u++)
    {
        owned_faces[u] = 0;
    }
    
    int* NewGlobVertCountPerRank = new int[world_size];
    int* NewGlobFaceCountPerRank = new int[world_size];

    for(int u=0;u<world_size;u++)
    {
        NewGlobVertCountPerRank[u] = 0;
        NewGlobFaceCountPerRank[u] = 0;
    }
        
    for(int i=0;i<Nt_shFaces;i++)
    {
        int key = TotalSharedFaces[i];
        int val = TotalSharedFacesRank[i];
        
        if(f2r_s.find(key)==f2r_s.end())
        {
            f2r_s.insert(key);
            f2r[key]=val;
            NewGlobFaceCountPerRank[val]=NewGlobFaceCountPerRank[val]+1;
        }
        else
        {
            tmp = f2r[key];
            if(val<tmp)
            {
                f2r[key]=val;
                owned_faces[val]=owned_faces[val]+1;

            }
            if(val>tmp)
            {
                f2r[key]=tmp;
                owned_faces[tmp]=owned_faces[tmp]+1;
            }
        }
    }
    
    std::map<int,int> v2r;
    std::set<int> v2r_s;
    int* owned_verts = new int[world_size];
    for(int u=0;u<world_size;u++)
    {
        owned_verts[u] = 0;
    }
    for(int i=0;i<Nt_shVerts;i++)
    {
        int key = TotalSharedVerts[i];
        int val = TotalSharedVertsRank[i];
        
        if(v2r_s.find(key)==v2r_s.end())
        {
            v2r_s.insert(key);
            v2r[key]=val;
            NewGlobVertCountPerRank[val]=NewGlobVertCountPerRank[val]+1;
        }
        else
        {
            tmp = v2r[key];
            if(val<tmp)
            {
                v2r[key]=val;
                owned_verts[val]=owned_verts[val]+1;

            }
            if(val>tmp)
            {
                v2r[key]=tmp;
                owned_verts[tmp]=owned_verts[tmp]+1;
            }
        }
    }

    int* NewGlobVertOffset = new int[world_size];
    int* NewGlobFaceOffset = new int[world_size];

    for(int u=1;u<world_size;u++)
    {
        NewGlobVertOffset[u] = NewGlobVertOffset[u-1]+NewGlobVertCountPerRank[u-1];
        NewGlobFaceOffset[u] = NewGlobFaceOffset[u-1]+NewGlobFaceCountPerRank[u-1];
    }
     
    std::map<int,std::vector<int> >::iterator itv;

    std::map<int,int> sharedVertsGlobal;
    
    for(int u=0;u<world_size;u++)
    {
        NewGlobVertCountPerRank[u] = 0;
        NewGlobFaceCountPerRank[u] = 0;
    }
    
    int iVshared = distLocalVerts->getNel()-nSharedVerts;
    
    std::map<int,int >::iterator itvv;
    std::map<int,int> sharedVmap;
    for(itvv=v2r.begin();itvv!=v2r.end();itvv++)
    {
        sharedVmap[itvv->first] = iVshared;
        iVshared++;
    }
    
    std::map<int,int> sharedFmap;
    int iFshared = distLocalFaces->getNel()-nSharedFaces;

    for(itvv=f2r.begin();itvv!=f2r.end();itvv++)
    {
        sharedFmap[itvv->first] = iFshared;
        iFshared++;
    }
    
    int nNonSharedVerts = nLocalVerts-owned_verts[world_rank];
    int nNonSharedFaces = nLocalFaces-owned_faces[world_rank];

    int* nNonSharedArray         = new int[world_size];
    int* nNonSharedArrayRed      = new int[world_size];
    int* nNonSharedVertsArrayOff = new int[world_size];

    int* nNonSharedFacesArray    = new int[world_size];
    int* nNonSharedFacesArrayRed = new int[world_size];
    int* nNonSharedFacesArrayOff = new int[world_size];

    for(i=0;i<world_size;i++)
    {
        nNonSharedArray[i]      = 0;
        nNonSharedFacesArray[i] = 0;
        if(i==world_rank)
        {
            nNonSharedArray[i] = nNonSharedVerts;
            nNonSharedFacesArray[i] = nNonSharedFaces;
        }
    }
    
    MPI_Allreduce(nNonSharedArray,
                  nNonSharedArrayRed,
                  world_size,
                  MPI_INT, MPI_SUM, comm);
    
    MPI_Allreduce(nNonSharedFacesArray,
                  nNonSharedFacesArrayRed,
                  world_size,
                  MPI_INT, MPI_SUM, comm);
    
    int nonSharedOff = 0;
    int nonFacesSharedOff = 0;
    for(i=0;i<world_size;i++)
    {
        nNonSharedVertsArrayOff[i] = nonSharedOff;
        nNonSharedFacesArrayOff[i] = nonFacesSharedOff;
        
        nonSharedOff = nonSharedOff + nNonSharedArrayRed[i];
        nonFacesSharedOff = nonFacesSharedOff + nNonSharedFacesArrayRed[i];
    }
    
    int nTet = 0;
    int allocRank = std::max_element(red_ielement_nlocs.begin(),red_ielement_nlocs.end())-red_ielement_nlocs.begin();
    
    if(world_rank == allocRank)
    {
        nTet = nTetras-NRankNoTets;
    }
    else
    {
        nTet=nTetras;
    }
    
    if(nTetras==0)
    {
        nTet = 1;
    }
    
    Array<int>* new_ien     = new Array<int>(nTet,4);
    Array<int>* new_ief     = new Array<int>(nTet,4);

    Array<int>* new_ien_or  = new Array<int>(nTet,4);
    Array<int>* new_ief_or  = new Array<int>(nTet,4);

    int elloc   = 0;
    int lbvid   = nNonSharedVertsArrayOff[world_rank];
    int lbfid   = nNonSharedFacesArrayOff[world_rank];
    int ltetvid = 0;
    
    std::set<int> gl_set;
    std::map<int,int> gl_map;
    int lbvids  = 0;
    int lvid2   = 0;
    
    std::set<int> glf_set;
    std::map<int,int> glf_map;
    int lbfids  = 0;
    int lfid2   = 0;
    
    int dd      = 0;
    int Nsend   = 0;
    int nv      = 4;
    if(world_rank == allocRank)
    {
        int u = 0;
        for(ite=tetras.begin();ite!=tetras.end();ite++)
        {
            int gEl = ite->first;
            
            int* ien   = new int[nv];
            int* ien_o = new int[nv];
            int* ief   = new int[nv];
            int* ief_o = new int[nv];
            
            for(int q=0;q<4;q++)
            {
                gvid = ite->second[q];
                
                if(v2r.find(gvid)!=v2r.end())
                {
                    lvid2 = sharedVmap[gvid];
                    ien[q] = lvid2;
                }
                else
                {
                    if(gl_set.find(gvid)==gl_set.end())
                    {
                        gl_set.insert(gvid);
                        gl_map[gvid] = lbvid;
                        lbvid = lbvid + 1;
                        ien[q] = lbvid;
                        lbvids++;
                    }
                    else
                    {
                        int lbbvid = gl_map[gvid];
                        ien[q] = lvid2;
                    }
                }
            }
            
            
            for(int q=0;q<4;q++)
            {
                gfid = ief_part_map->i_map[gEl][q];
                
                if(f2r.find(gfid)!=f2r.end())
                {
                    lfid2 = sharedFmap[gfid];
                    ief[q] = lfid2;
                }
                else
                {
                    if(glf_set.find(gfid)==glf_set.end())
                    {
                        glf_set.insert(gfid);
                        glf_map[gfid] = lbfid;
                        lbfid = lbfid + 1;
                        lbfids++;
                    }
                    else
                    {
                        int lbbfid = glf_map[gfid];
                        ief[q] = lbbfid;
                    }
                }
            }
            
            if(u<world_size && red_Elplease[u]==1)
            {
                ien_o[0]   = ite->second[0];
                ien_o[1]   = ite->second[1];
                ien_o[2]   = ite->second[2];
                ien_o[3]   = ite->second[3];
                
                ief_o[0]   = ief_part_map->i_map[gEl][0];
                ief_o[1]   = ief_part_map->i_map[gEl][1];
                ief_o[2]   = ief_part_map->i_map[gEl][2];
                ief_o[3]   = ief_part_map->i_map[gEl][3];
                
                MPI_Send(&ien[0],   nv, MPI_INT, u,   u*64975,  comm);
                MPI_Send(&ien_o[0], nv, MPI_INT, u,   u*84975,  comm);
                MPI_Send(&ief[0],   nv, MPI_INT, u,   u*94975,  comm);
                MPI_Send(&ief_o[0], nv, MPI_INT, u,  u*104975, comm);
                
                Nsend++;
            }
            else
            {
                new_ien->setVal(elloc,0,ien[0]);
                new_ien->setVal(elloc,1,ien[1]);
                new_ien->setVal(elloc,2,ien[2]);
                new_ien->setVal(elloc,3,ien[3]);
                
                new_ief->setVal(elloc,0,ief[0]);
                new_ief->setVal(elloc,1,ief[1]);
                new_ief->setVal(elloc,2,ief[2]);
                new_ief->setVal(elloc,3,ief[3]);
                
                elloc++;
            }
            u++;
        }
        
        std::cout << "elloc " << elloc << " " << nTetras << " " << NRankNoTets << std::endl;
    }
    else
    {
        std::vector<int> ien(4);
        std::vector<int> ien_o(4);
        
        std::vector<int> ief(4);
        std::vector<int> ief_o(4);

        if(red_Elplease[world_rank]==1)
        {
            MPI_Recv(&ien[0], 4, MPI_INT, allocRank, world_rank*64975, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&ien_o[0], 4, MPI_INT, allocRank, world_rank*84975, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&ief[0], 4, MPI_INT, allocRank, world_rank*94975, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&ief_o[0], 4, MPI_INT, allocRank, world_rank*104975, comm, MPI_STATUS_IGNORE);
            
            new_ien->setVal(0,0,ien[0]);
            new_ien->setVal(0,1,ien[1]);
            new_ien->setVal(0,2,ien[2]);
            new_ien->setVal(0,3,ien[3]);
            
            new_ien_or->setVal(0,0,ien_o[0]);
            new_ien_or->setVal(0,1,ien_o[1]);
            new_ien_or->setVal(0,2,ien_o[2]);
            new_ien_or->setVal(0,3,ien_o[3]);
            
            new_ief->setVal(0,0,ief[0]);
            new_ief->setVal(0,1,ief[1]);
            new_ief->setVal(0,2,ief[2]);
            new_ief->setVal(0,3,ief[3]);
            
            new_ief_or->setVal(0,0,ief_o[0]);
            new_ief_or->setVal(0,1,ief_o[1]);
            new_ief_or->setVal(0,2,ief_o[2]);
            new_ief_or->setVal(0,3,ief_o[3]);
        }
        else
        {
            for(ite=tetras.begin();ite!=tetras.end();ite++)
            {
                int gEl = ite->first;
                
                for(int q=0;q<4;q++)
                {
                    gvid = ite->second[q];
                    
                    if(v2r.find(gvid)!=v2r.end())
                    {
                        lvid2 = sharedVmap[gvid];
                        new_ien->setVal(elloc,q,lvid2);
                    }
                    else
                    {
                        if(gl_set.find(gvid)==gl_set.end())
                        {
                            gl_set.insert(gvid);
                            new_ien->setVal(elloc,q,lbvid);
                            gl_map[gvid] = lbvid;
                            lbvid = lbvid + 1;
                            lbvids++;
                        }
                        else
                        {
                            int lbbvid = gl_map[gvid];
                            new_ien->setVal(elloc,q,lbbvid);
                        }
                    }
                }
                
                
                for(int q=0;q<4;q++)
                {
                    gfid = ief_part_map->i_map[gEl][q];
                    
                    if(f2r.find(gfid)!=f2r.end())
                    {
                        lfid2 = sharedFmap[gfid];
                        new_ief->setVal(elloc,q,lfid2);
                    }
                    else
                    {
                        if(glf_set.find(gfid)==glf_set.end())
                        {
                            glf_set.insert(gfid);
                            new_ief->setVal(elloc,q,lbfid);
                            glf_map[gfid] = lbfid;
                            lbfid = lbfid + 1;
                            lbfids++;
                        }
                        else
                        {
                            int lbbfid = glf_map[gfid];
                            new_ief->setVal(elloc,q,lbbfid);
                        }
                    }
                }
                elloc++;
            }
        }
    }
    
    int* elmdist              = new int[world_size+1];

    for(i=0;i<world_size;i++)
    {
        ielement_nlocs[i]     = 0;
        red_ielement_nlocs[i] = 0;

        if(i==world_rank)
        {
            ielement_nlocs[i] = nTet;
        }
        else
        {
            ielement_nlocs[i] = 0;
        }
    }
    
    MPI_Allreduce(ielement_nlocs, &red_ielement_nlocs[0], world_size, MPI_INT, MPI_SUM, comm);
    
    int o_ie = 0;
    
    for(i=0;i<world_size;i++)
    {
        ielement_offsets[i] = o_ie;
        ielement_nlocs[i] = red_ielement_nlocs[i];
        elmdist[i]          = o_ie;
        o_ie                = o_ie+red_ielement_nlocs[i];
    }
    
    elmdist[world_size] = o_ie;

    int* eptr     = new int[nTet+1];
    int* eind     = new int[nTet*4];

    eptr[0]  = 0;
    for(int i=0;i<nTet;i++)
    {
        eptr[i+1] = eptr[i]+4;
        for(int j=eptr[i];j<eptr[i+1];j++)
        {
            eind[j]    = new_ien->data[j];
        }
    }

    idx_t numflag_[] = {0};
    idx_t *numflag = numflag_;
    idx_t ncommonnodes_[] = {3};
    idx_t *ncommonnodes = ncommonnodes_;
    idx_t *xadj_par      = NULL;
    idx_t *adjncy_par    = NULL;
    int np           = world_size;
    idx_t ncon_[]    = {1};
    int edgecut      = 0;
    idx_t *ncon      = ncon_;
    real_t *tpwgts   = new real_t[np*ncon[0]];
    idx_t wgtflag_[] = {2};
    idx_t *wgtflag   = wgtflag_;
    real_t ubvec_[]  = {1.1};
    real_t *ubvec    = ubvec_;
    idx_t options_[] = {0, 0, 0};
    idx_t *options   = options_;
    int* part_arr = new int[nTet];
    idx_t nparts_[] = {np};
    idx_t *nparts = nparts_;
    for(int i=0; i<np*ncon[0]; i++)
    {
        tpwgts[i] = 1.0/np;
    }

    int *elmwgt = new int[nTet];
    for(int i=0;i<nTet;i++)
    {
        elmwgt[i] = 1;
    }

    ParMETIS_V3_Mesh2Dual(elmdist,
                          eptr,
                          eind,
                          numflag,ncommonnodes,
                          &xadj_par,&adjncy_par,&comm);

    ParMETIS_V3_PartKway(elmdist,
                         xadj_par,
                         adjncy_par,
                         elmwgt, NULL, wgtflag, numflag,
                         ncon, nparts,
                         tpwgts, ubvec, options,
                         &edgecut, part_arr, &comm);

    Array<int>* part_global_new = new Array<int>(o_ie,1);
    Array<int>* part_new = new Array<int>(nTet,1);

    part_new->data = part_arr;
    
    MPI_Allgatherv(&part_new->data[0],
                   red_ielement_nlocs[world_rank], MPI_INT,
                   &part_global_new->data[0],
                   ielement_nlocs,
                   ielement_offsets,
                   MPI_INT,comm);
    
    if(world_rank == 0)
    {
        int* nElpRank = new int[world_size];
        for(int u=0;u<world_size;u++)
        {
            nElpRank[u] = 0;
        }
        for(int u=0;u<o_ie;u++)
        {
            int pid = part_global_new->getVal(u,0);
            nElpRank[pid] = nElpRank[pid]+1;
        }
        
        std::cout << "Printing out the new distribution ::" << std::endl;
        for(int u=0;u<world_size;u++)
        {
            std::cout << nElpRank[u] << " " << o_ie <<  std::endl;
        }
    }
    
    
    /**/
    /**/
    
    
//    if(world_rank == 3)
//    {
//        std::map<int,int >::iterator itsf;
//        for(itsf=sharedFacesGlobal.begin();itsf!=sharedFacesGlobal.end();itsf++)
//        {
//            std::cout << itsf->first << " " << itsf->second << std::endl;
//        }
//    }
    
    //std::map<int,std::vector<int> >::iterator itsf;
    
//    for(itsf=recv_SharedFaces.begin();itsf!=recv_SharedFaces.end();itsf++)
//    {
//        std::cout << "recv sizing " << world_rank << " " << itsf->first << " " << itsf->second.size() << std::endl;
//    }
    //std::cout << "recv sizing " << world_rank << " " << recv_SharedFaces.size() << " " << recv_SharedFaces.size() << std::endl;
    
//    ScheduleObj* sObj = DoScheduling(send_SharedFaces,comm);
//
//    std::map<int, std::set<int> >::iterator its;
    
//    if(world_rank == 0)
//    {
//        std::map<int,std::vector<int> >::iterator itm;
//        for(itm=send_SharedFaces.begin();itm!=send_SharedFaces.end();itm++)
//        {
//            std::cout << "Rank  " << world_rank << " sends the following elements: ";
////            for(int q=0;q<itm->second.size();q++)
////            {
////                std::cout <<itm->second[q] << " ";
////            }
//            std::cout << " to rank " << itm->first << std::endl;
//        }
//        //std::cout << "Send Schedule " << sObj->SendFromRank2Rank.size() <<std::endl;
//
//        for(its=sObj->SendFromRank2Rank.begin();its!=sObj->SendFromRank2Rank.end();its++)
//        {
//            std::cout << world_rank << " -> " << its->first << " gets " << its->second.size() << " elements ";
//            std::set<int>::iterator it;
////            for(it=its->second.begin();it!=its->second.end();it++)
////            {
////                std::cout << *it << " ";
////            }
////            std::cout << std::endl;
//        }
////
////
//        //std::cout << "Recv Schedule " << sObj->SendFromRank2Rank.size() <<std::endl;
//
//        for(its=sObj->RecvRankFromRank.begin();its!=sObj->RecvRankFromRank.end();its++)
//        {
//            std::cout << "Rank  " << world_rank << " recieved the following elements: ";
////            for(int q=0;q<itm->second.size();q++)
////            {
////                std::cout <<itm->second[q] << " ";
////            }
//            std::cout << " from rank " << itm->first << std::endl;
//        }
//    }
//
//
//    std::cout << "======================================================="<<std::endl;
//    for(itsf=send_SharedFaces.begin();itsf!=send_SharedFaces.end();itsf++)
//    {
//        std::cout << "send sizing " << world_rank << " " << itsf->first << " " << itsf->second.size() << std::endl;
//    }
    //std::cout << "send sizing " << world_rank << " " << recv_SharedFaces.size() << " " << recv_SharedFaces.size() << std::endl;
    
    
    
//    int* iface_nlocs       = new int[world_size];
//    int* red_iface_nlocs   = new int[world_size];
//    int* iface_offsets     = new int[world_size];
//
//    int* uface_nlocs       = new int[world_size];
//    int* red_uface_nlocs   = new int[world_size];
//    int* uface_offsets     = new int[world_size];
//
//    for(i=0;i<world_size;i++)
//    {
//        iface_nlocs[i]     = 0;
//        red_iface_nlocs[i] = 0;
//
//        uface_nlocs[i]     = 0;
//        red_uface_nlocs[i] = 0;
//
//        if(i==world_rank)
//        {
//            iface_nlocs[i] = uInterFaces.size();
//            uface_nlocs[i] = ufaceOnRank.size();
//        }
//        else
//        {
//            iface_nlocs[i] = 0;
//            uface_nlocs[i] = 0;
//
//        }
//    }
//
//    MPI_Allreduce(iface_nlocs, red_iface_nlocs, world_size, MPI_INT, MPI_SUM, comm);
//    MPI_Allreduce(uface_nlocs, red_uface_nlocs, world_size, MPI_INT, MPI_SUM, comm);
//
//    int o_if = 0;
//    int o_uf = 0;
//
//    for(i=0;i<world_size;i++)
//    {
//        iface_offsets[i] = o_if;
//        o_if = o_if+red_iface_nlocs[i];
//
//        uface_offsets[i] = o_uf;
//        o_uf = o_uf+red_uface_nlocs[i];
//
//    }
//
//    int n_glob_if = o_if;
//
//    std::vector<int> intFace(n_glob_if);
//
//    MPI_Allgatherv(&uSharedFaces[0],
//                   uSharedFaces.size(),
//                   MPI_INT,
//                   &intFace[0],
//                   red_iface_nlocs,
//                   iface_offsets,
//                   MPI_INT, comm);
//
//    std::map<int,int> interFaceMap;
//    int val = 0;
//    for(int u=0;u<intFace.size();u++)
//    {
//        int key = intFace[u];
//        if(interFaceMap.find(key)==interFaceMap.end())
//        {
//            interFaceMap[key] = val;
//            val++;
//        }
//    }
    /* */
    // Communicate the face_color_map;
    //MeshTransfer* TrObj = GetMeshTransfer(tetras,ufaces_vec,uverts_vec,comm,info);
    
//    for(int q=0;q<TrObj->UniqueFacesVec.size();q++)
//    {
//        //std::cout << q << " " << TrObj->UniqueFacesVec[q]<<std::endl;
//    }
//    std::map<int,int> glob2loc_PartInterF = TrObj->iFaces->glob2loc_PartInterEntity;
//    std::set<int> uniqueFaces = TrObj->iFaces->UniqueFaces;
//    tetras     = pDom->GTetras;
//    int uf = 0;
//    ufaces.clear();
//    std::set<int> ufaces2;
//    int duppe3 = 0;
//    for(ite=tetras.begin();ite!=tetras.end();ite++)
//    {
//        int gEl = ite->first;
//
//        for(int j=0;j<4;j++)
//        {
//            gfid = ief_part_map->i_map[gEl][j];
//
//            if(glob2loc_PartInterF.find(gfid)!=glob2loc_PartInterF.end())
//            {
//                duppe3++;
//                ufaces2.insert(gfid);
//            }
//
//        }
//    }
//
//
//    if(world_rank == 0)
//    {
////        std::map<int,int>::iterator its;
////        int dupl = 0;
////        for(its=glob2loc_PartInterF.begin();its!=glob2loc_PartInterF.end();its++)
////        {
////            std::cout << its->first << " " << its->second << std::endl;
////        }
////
////        std::set<int>::iterator itss;
////        int ig = 0;
////        for(itss=ufaces_comp.begin();itss!=ufaces_comp.end();itss++)
////        {
////            std::cout << world_rank << " " << ig << " " << *itss << std::endl;
////            ig++;
////        }
//
//    }
//
//
//
//
//    std::cout << "uf = " << world_rank <<  " " << ufaces.size() << " " << glob2loc_PartInterF.size() << " " << ufaces_VV.size() << " " << ufaces2.size() << " " << duppe3 << " " << ufaces2.size() << " dupl = " << " " << ufaces_comp.size() << " " << uniqueFaces.size() << std::endl;
//
//
    
    
//    std::map<int,int> glob2loc_PartInterV = TrObj->iVerts->glob2loc_PartInterEntity;
//
//    std::cout << "rankert = " << world_rank << " " << duppe << " " << glob2loc_PartInterF.size() << std::endl;
    
    //Array<int>* xcn_tet = new Array<int>(uverts_vec.size(),3);
//    for(ite=tetras.begin();ite!=uverts_vec.end();ite++)
//    {
//        int vid = uverts_vec[i];
//
//        xcn_tet->setVal(i,0,Verts[vid]->x);
//        xcn_tet->setVal(i,1,Verts[vid]->y);
//        xcn_tet->setVal(i,2,Verts[vid]->z);
//    }
    
    
//    Array<int>* ien_tet = new Array<int>(tetras.size(),4);
//    Array<int>* iee_tet = new Array<int>(tetras.size(),6);
//    Array<int>* ief_tet = new Array<int>(tetras.size(),4);
//    ufaces.clear();
//    uverts.clear();
//    int lnfid,lnvid;
//    i = 0;
//    for(ite=tetras.begin();ite!=tetras.end();ite++)
//    {
//        int gEl = ite->first;
//
//        for(int j=0;j<4;j++)
//        {
//            gfid = ief_part_map->i_map[gEl][j];
//            lnfid = glob2loc_PartInterF[gfid];
//            ief_tet->setVal(i,j,lnfid);
//
//            gvid = ien_part_map->i_map[gEl][j];
//            lnvid = glob2loc_PartInterV[gfid];
//            ien_tet->setVal(i,j,lnvid);
//        }
//
//        i++;
//    }
    
    
    
//
    
//+
    
    
    //    Array<int>* ien_tet = GatherTetrahedraOnRoot(tetras,comm,info);

    
    
    /**/
//    int loc_iftri = 0;
//    std::map<int,std::vector<int> >::iterator itm;
//    std::vector<int> interfaces;
//    for(itm  = face_color_map.begin();
//        itm != face_color_map.end();
//        itm++)
//    {
//
//        //std::cout << "R = " << world_rank << " " << itm->second.size() << std::endl;
//
//
//        loc_iftri = loc_iftri + itm->second.size();
//
//        for(int u=0;u<itm->second.size();u++)
//        {
//            interfaces.push_back(itm->second[u]);
//        }
//
//    }
//    //std::cout << "rank " << world_rank << " " << loc_iftri << std::endl;
//
//    std::vector<int> duple = FindDuplicates(interfaces);
    //std::cout << "RR " << world_rank << " " << duple.size() << std::endl;
    
    
    
    /*
    
    std::map<int,std::vector<int> >::iterator itb;
    
    int nbcface     = 0;
    int loc_bftri   = 0;
    int loc_bfqua   = 0;
    int bfaceid     = 0;
    int NvpF        = 0;
    
    // Counting boundary triangles.
    int gv0,gv1,gv2;
    std::vector<std::vector<int> > faces_part;
    std::vector<int> faces_ref;
    for(itb=ref2bcface.begin();itb!=ref2bcface.end();itb++)
    {
        ref     = itb->first;

        for(int fid=0;fid<itb->second.size();fid++)
        {
            bfaceid = itb->second[fid];
            NvpF    = if_Nv_part_map->i_map[bfaceid][0];
            gv0 = gv2lpv[ifn_part_map->i_map[bfaceid][0]];
            gv1 = gv2lpv[ifn_part_map->i_map[bfaceid][1]];
            gv2 = gv2lpv[ifn_part_map->i_map[bfaceid][2]];
            
            std::vector<int> face(3);
            face[0] = gv0;
            face[1] = gv1;
            face[2] = gv2;
            faces_part.push_back(face);
            faces_ref.push_back(ref);

            loc_bftri++;
        }
    }
    
    
    
    int* bftris     = new int[world_size];
    int* red_bftris = new int[world_size];
    int* off_bftris = new int[world_size];
    
    
    int* iftris     = new int[world_size];
    int* red_iftris = new int[world_size];
    int* off_iftris = new int[world_size];

    
    int offset_inter;
    int loc_iftri = 0;
    
    // Counting internal partition triangles.

    std::map<int,std::vector<int> >::iterator itm;
    for(itm  = face_color_map.begin();
        itm != face_color_map.end();
        itm++)
    {
        loc_iftri = loc_iftri + itm->second.size();
    }

    for(int i=0;i<world_size;i++)
    {
        bftris[i]=0;
        iftris[i]=0;

        
        if(i==world_rank)
        {
            bftris[i] = loc_bftri;
            iftris[i] = loc_iftri;
        }
    }
    
    
    MPI_Allreduce(iftris,  red_iftris,   world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(bftris,  red_bftris,   world_size, MPI_INT, MPI_SUM, comm);
    
    int offset_btri = 0;
    int offset_itri = 0;
    
    for(int i=0;i<world_size;i++)
    {
        off_bftris[i] = offset_btri;
        off_iftris[i] = offset_itri;
        
        offset_btri   = offset_btri    +  red_bftris[i];
        offset_itri   = offset_itri   +   red_iftris[i];
    }

    int offset_tri_glob = offset_btri;
    int ncomm           = ordered_rank.size();
    int* color_face     = new int[ncomm];
    std::vector<int* > ifc_tria_glob(ncomm);
    std::vector<int* > ifc_tria_loc(ncomm);

    int* ntifc = new int[ncomm];
    std::set<int>::iterator it;
    int u = 0;
    for(itm  = face_color_map.begin();
        itm != face_color_map.end();
        itm++)
    {
        ra              = itm->first;
        pos             = std::distance(ordered_rank.begin(),ordered_rank.find(ra));
        color_face[pos] = ra;
        ntifc[pos]      = itm->second.size();

        int* ifc_tria_glob_pos = new int[ntifc[pos]];
        int* ifc_tria_loc_pos = new int[ntifc[pos]];

//        std::cout << "rank =  " << rank << " ncomm = "
//         << ncomm << " pos = " << pos << " ra = " << ra
//         << " nF = " << itm->second.size() << std::endl;
        
        for(int i=0;i<itm->second.size();i++)
        {
            std::vector<int> face(3);
            face[0] = gv2lpv[ifn_part_map->i_map[itm->second[i]][0]];
            face[1] = gv2lpv[ifn_part_map->i_map[itm->second[i]][1]];
            face[2] = gv2lpv[ifn_part_map->i_map[itm->second[i]][2]];
            faces_part.push_back(face);
            faces_ref.push_back(2);
            ifc_tria_glob_pos[i] = offset_tri_glob+off_iftris[world_rank]+u;
            ifc_tria_loc_pos[i]  = red_bftris[world_rank]+u;
            u++;
        }
        ifc_tria_glob[pos] = ifc_tria_glob_pos;
        ifc_tria_loc[pos] = ifc_tria_loc_pos;
    }
    
    
    
    
    
    PMMG_pParMesh   parmesh;
    PMMG_Init_parMesh(PMMG_ARG_start,
                      PMMG_ARG_ppParMesh,&parmesh,
                      PMMG_ARG_pMesh,PMMG_ARG_pMet,
                      PMMG_ARG_dim,3,PMMG_ARG_MPIComm,MPI_COMM_WORLD,
                      PMMG_ARG_end);
    */
    
    
    
    

//==================================================================================
//==================================================================================
//==================================================================================
//    int nVertices       = lv2gpv.size();
//    int nTetrahedra     = tetras.size();
//    int nTriangles      = tris.size();
//    int nEdges          = 0;
//    int nPrisms         = 0;
//    int nQuadrilaterals = bfqua;
//
//    std::map<int,int> gv2lpartv = pDom->gv2lpartv;
//    //std::cout << world_rank << " " << bftri << " " << bfqua << std::endl;
//
//    if ( PMMG_Set_meshSize(parmesh,nVertices,nTetrahedra,nPrisms,nTriangles,
//                              nQuadrilaterals,nEdges) != 1 ) {
//      MPI_Finalize();
//      exit(EXIT_FAILURE);
//    }
//    const int nSolsAtVertices = 1; // 3 solutions per vertex
//    int solType[1] = {MMG5_Tensor};
//
//    if ( PMMG_Set_solsAtVerticesSize(parmesh,nSolsAtVertices,nVertices,solType) != 1 ) {
//       MPI_Finalize();
//       exit(EXIT_FAILURE);
//    }
//
//
//    double* tensor = new double[6];
//    for ( k=0; k<nVertices; ++k )
//    {
//          int gvid = lv2gpv[k];
//          int lvid_p = gv2lpartv[gvid];
//
//          double vx = verts[lvid_p]->x;
//          double vy = verts[lvid_p]->y;
//          double vz = verts[lvid_p]->z;
//
//          if ( PMMG_Set_vertex(parmesh,vx,vy,vz, 1.0, k+1) != 1 )
//          {
//            MPI_Finalize();
//            exit(EXIT_FAILURE);
//          }
//
//
//
//          tensor[0] = hess_vmap[gvid]->getVal(0,0);
//          tensor[1] = hess_vmap[gvid]->getVal(0,1);
//          tensor[2] = hess_vmap[gvid]->getVal(0,2);
//          tensor[3] = hess_vmap[gvid]->getVal(1,1);
//          tensor[4] = hess_vmap[gvid]->getVal(1,2);
//          tensor[5] = hess_vmap[gvid]->getVal(2,2);
//
//
//          if ( PMMG_Set_ithSol_inSolsAtVertices(parmesh,1,tensor,k+1) != 1 )
//          {
//             MPI_Finalize();
//             exit(EXIT_FAILURE);
//          }
//
//
//
//    }
//
//
//
//
//    std::map<int,std::vector<int> >::iterator ittet;
//    k = 0;
//    for ( ittet=tetras.begin(); ittet!=tetras.end(); ittet++ )
//    {
//    	if ( PMMG_Set_tetrahedron(parmesh,ittet->second[0]+1,ittet->second[1]+1,
//        		ittet->second[2]+1,ittet->second[3]+1,1.0,k+1) != 1 ) {
//            MPI_Finalize();
//            exit(EXIT_FAILURE);
//          }
//        k++;
//    }
//
//
//
//
//
//
//
//    std::map<int,std::vector<int> >::iterator itm;
//    k = 0;
//    for(int q=0;q<tris.size();q++)
//    {
//        int vt0 = tris[q][0];
//        int vt1 = tris[q][1];
//        int vt2 = tris[q][2];
//        int ref = tris_ref[q];
//
//        if ( PMMG_Set_triangle(parmesh, vt0+1,vt1+1,vt2+1, ref,k+1) != 1 )
//        {
//            MPI_Finalize();
//            exit(EXIT_FAILURE);
//        }
//        k++;
//    }
//
//    int API_mode = PMMG_APIDISTRIB_faces;
//
//    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_APImode, 0 ) )
//    {
//       MPI_Finalize();
//       exit(EXIT_FAILURE);
//    };
//
//    if( !world_rank )
//    printf("\n--- API mode: Setting face communicators\n");
//
//    /* Set the number of interfaces */
//    //ier = PMMG_Set_numberOfFaceCommunicators(parmesh, ncomm);
//    int n_face_comm;
//    //ier = PMMG_Set_numberOfFaceCommunicators(parmesh, n_face_comm);
//    int* color_face     = pDom->color_face;
//    std::vector<int* > ifc_tria_loc  = pDom->ifc_tria_loc;
//    std::vector<int* > ifc_tria_glob = pDom->ifc_tria_glob;
//    int* ntifc          = pDom->ntifc;
//    int ncomm           = pDom->ncomm;
//    ier = PMMG_Set_numberOfFaceCommunicators(parmesh, ncomm);
//
//    //std::cout << n_face_comm << std::endl;
//    /* Loop on each interface (proc pair) seen by the current rank) */
//    for(int icomm=0; icomm<ncomm; icomm++ )
//    {
//
//        std::cout << "rank = " << world_rank << " " << color_face[icomm] << " " << ntifc[icomm] << " " << nTriangles << " bftrie = " << bftri << std::endl;
//    	/* Set nb. of entities on interface and rank of the outward proc */
//    	ier = PMMG_Set_ithFaceCommunicatorSize(parmesh, icomm,
//										   color_face[icomm],
//										   ntifc[icomm]);
//
//
////        // Set local and global index for each entity on the interface
//    	ier = PMMG_Set_ithFaceCommunicator_faces(parmesh, icomm,
//                                                 ifc_tria_loc[icomm],
//                                                 ifc_tria_glob[icomm], 1 );
//    }
//
//    int niter = 1;
//
//    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_niter, niter ) )
//    {
//       MPI_Finalize();
//       exit(EXIT_FAILURE);
//    };
//
//
//    std::cout << "starting the parallel remeshing..." << std::endl;
//    int ierlib = PMMG_parmmglib_distributed( parmesh );
//    std::cout << "finishing the parallel remeshing..." << std::endl;
//
//
    
    MPI_Finalize();
    
}

