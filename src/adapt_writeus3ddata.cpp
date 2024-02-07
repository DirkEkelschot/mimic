#include "adapt_writeus3ddata.h"



void WritePrismsUS3DFormat(MPI_Comm comm, 
                           RepartitionObject* prism_repart,
                           PrismTetraTrace* pttrace,
                           std::map<int,int> tracetagV2globalV,
                           std::vector<int> &ifn_P,
                           std::map<int,std::vector<int> > ranges_id)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    std::map<int,std::vector<int> >::iterator itmiv;
    std::map<int,std::vector<int> > InternalFaces_P         = prism_repart->getInteriorFaceMap();
    std::map<int,std::vector<int> > SharedFaces_P           = prism_repart->getSharedFaceMap();
    std::map<int, std::vector<double> > LocalVertsMap_P     = prism_repart->getLocalVertsMap();
    // std::map<int,int> tag2globvID_P                      = prism_repart->getUpdatedTag2GlobalVMap();
    std::map<int,int> SharedVertsNotOwned                   = prism_repart->getSharedVertsNotOwned();
    std::map<int,int> NonSharedVertsOwned                   = prism_repart->getNonSharedVertsOwned();
    std::map<int,int> SharedVertsOwned                      = prism_repart->getSharedVertsOwned();
    std::map<int,std::vector<int> > ElementTag2VertexTag_P  = prism_repart->getElement2VertexMap();

    std::map<int,int>  lh_P                                 = prism_repart->GetLeftHandFaceElementMap();
    std::map<int,int>  rh_P                                 = prism_repart->GetRightHandFaceElementMap();
    std::map<int,int> gE2tagE_P                             = prism_repart->getGlobalElement2ElementTag();
    
    std::map<int,int> NonSharedVertsOwned_P                 = prism_repart->getNonSharedVertsOwned();
    std::map<int,int> SharedVertsOwned_P                    = prism_repart->getSharedVertsOwned();


    std::map<int,int>::iterator itmii;
    int vert = 0;
    std::map<int,int> tag2globvID_P;
    // std::map<int,int>::iterator itmii;

    int nLocIntVrts         = NonSharedVertsOwned_P.size();
    int nLocShVrts          = SharedVertsOwned_P.size();
    DistributedParallelState* distnLocIntVrts = new DistributedParallelState(nLocIntVrts,comm);
    DistributedParallelState* distnLocShVrts  = new DistributedParallelState(nLocShVrts,comm);

    for(itmii=NonSharedVertsOwned_P.begin();itmii!=NonSharedVertsOwned_P.end();itmii++)
    {
        int tag = itmii->first;
        tag2globvID_P[tag]              = vert+distnLocIntVrts->getOffsets()[world_rank]+1;        
        vert++;
    }


    int vertsh = 0;
    
    for(itmii=SharedVertsOwned_P.begin();itmii!=SharedVertsOwned_P.end();itmii++)
    {
        int tag                         = itmii->first;
        tag2globvID_P[tag]              = SharedVertsOwned_P[tag];
        vertsh++;
    }

    std::cout << "SharedFaces_P " << SharedFaces_P.size() << " " << world_rank << std::endl; 
    std::cout << "lh_P " << lh_P.size() << " rank = " << world_rank << std::endl;
    std::cout << "lr_P " << rh_P.size() << " rank = " << world_rank << std::endl;
    std::vector<double> VcF(3,0);

    std::map<int,int> unique_trace_verts2refmap     = pttrace->GetUniqueTraceVerts2RefMap();
    int nLocFace_P = InternalFaces_P.size()+SharedFaces_P.size();
    ifn_P.resize(nLocFace_P*8);
    std::fill(ifn_P.begin(), ifn_P.end(), 0);

    int fptot = 0;
    for( itmiv=InternalFaces_P.begin();itmiv!=InternalFaces_P.end();itmiv++)
    {
        int gfid    = itmiv->first;
        int npf     = itmiv->second.size();
        std::vector<int> fce(npf);
        ifn_P[fptot*8+0] = npf;

        std::vector<std::vector<double> > Vfaces;
        VcF[0] = 0.0;
        VcF[1] = 0.0;
        VcF[2] = 0.0;
        

        for(int g=0;g<npf;g++)
        {
            int tagV = itmiv->second[g];
            
            std::vector<double> Vf(3,0.0);
            Vf[0] = LocalVertsMap_P[tagV][0];
            Vf[1] = LocalVertsMap_P[tagV][1];
            Vf[2] = LocalVertsMap_P[tagV][2];
            Vfaces.push_back(Vf);

            VcF[0] = VcF[0] + Vf[0];
            VcF[1] = VcF[1] + Vf[1];
            VcF[2] = VcF[2] + Vf[2];
            
            if(tag2globvID_P.find(tagV)!=tag2globvID_P.end() &&
                unique_trace_verts2refmap.find(tagV)==unique_trace_verts2refmap.end())
            {
                int globid  = tag2globvID_P[tagV];
                fce[g]      = globid;
            
            }
            else if(SharedVertsNotOwned.find(tagV)!=SharedVertsNotOwned.end() &&
                    unique_trace_verts2refmap.find(tagV)==unique_trace_verts2refmap.end())
            {
                int globid  = SharedVertsNotOwned[tagV];
                fce[g]      = globid;   
            }
            else if(unique_trace_verts2refmap.find(tagV)!=unique_trace_verts2refmap.end())
            {
                int ref     = unique_trace_verts2refmap[tagV];
                int globid  = tracetagV2globalV[ref];
                fce[g]      = globid;
            }
        }

        VcF[0] = VcF[0]/npf;
        VcF[1] = VcF[1]/npf;
        VcF[2] = VcF[2]/npf;

        if(rh_P[gfid] == 0 || lh_P[gfid] == 0)
        {
            std::cout << world_rank << " rhp[gfid] and lhp[gfid] are zero INTERNAL " << std::endl;
        }

        int leftEl  = lh_P[gfid];
        int leftTag = gE2tagE_P[leftEl];
        
        std::vector<double> Vijk(3,0);
        Vijk[0] = 0.0;
        Vijk[1] = 0.0;
        Vijk[2] = 0.0;
        // compute element center;
        int nvp = ElementTag2VertexTag_P[leftTag].size();

        for(int q=0;q<nvp;q++)
        {
            int tagV  = ElementTag2VertexTag_P[leftTag][q];

            Vijk[0] = Vijk[0] + LocalVertsMap_P[tagV][0];
            Vijk[1] = Vijk[1] + LocalVertsMap_P[tagV][1];
            Vijk[2] = Vijk[2] + LocalVertsMap_P[tagV][2];
        }

        Vijk[0] = Vijk[0]/nvp;
        Vijk[1] = Vijk[1]/nvp;
        Vijk[2] = Vijk[2]/nvp;

        double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);
        if(orient0 < 0.0)
        {
            if(npf == 3)
            {
                ifn_P[fptot*8+1] = fce[0];
                ifn_P[fptot*8+2] = fce[2];
                ifn_P[fptot*8+3] = fce[1];
                ifn_P[fptot*8+4] = 0;
            }
            if(npf == 4)
            {
                ifn_P[fptot*8+1] = fce[0];
                ifn_P[fptot*8+2] = fce[3];
                ifn_P[fptot*8+3] = fce[2];
                ifn_P[fptot*8+4] = fce[1];
            }
        }
        else
        {
            if(npf == 3)
            {
                ifn_P[fptot*8+1] = fce[0];
                ifn_P[fptot*8+2] = fce[1];
                ifn_P[fptot*8+3] = fce[2];
                ifn_P[fptot*8+4] = 0;
            }
            if(npf == 4)
            {
                ifn_P[fptot*8+1] = fce[0];
                ifn_P[fptot*8+2] = fce[1];
                ifn_P[fptot*8+3] = fce[2];
                ifn_P[fptot*8+4] = fce[3];
            }
        }

        // if(lh_P[gfid]>279070 || rh_P[gfid]>279070)
        // {
        //     std::cout << lh_P[gfid] << " " << rh_P[gfid] << " InternalFaces_P pp " << std::endl;
        // }
        
        for(int h=0;h<fce.size();h++)
        {
            if(fce[h]>162713)
            {
                std::cout << "Error line 167 Interior " << fce[h] << " " << 162713 << std::endl;
            }
        }

        ifn_P[fptot*8+5] = rh_P[gfid];
        ifn_P[fptot*8+6] = lh_P[gfid];
        ifn_P[fptot*8+7] = 2;

        fptot++;
    }


    std::cout << "InteriorFaces_P " << InternalFaces_P.size() << std::endl;






    int ishare = 0;

    std::cout << "SharedFaces_P " << SharedFaces_P.size() << std::endl;

    for( itmiv=SharedFaces_P.begin();itmiv!=SharedFaces_P.end();itmiv++)
    {
        int gfid    = itmiv->first;
        int npf     = itmiv->second.size();
        std::vector<int> fce(npf);
        ifn_P[fptot*8+0] = npf;

        std::vector<std::vector<double> > Vfaces;
        VcF[0] = 0.0;
        VcF[1] = 0.0;
        VcF[2] = 0.0;


        for(int g=0;g<npf;g++)
        {
            int tagV = itmiv->second[g];
            
            std::vector<double> Vf(3,0.0);
            Vf[0] = LocalVertsMap_P[tagV][0];
            Vf[1] = LocalVertsMap_P[tagV][1];
            Vf[2] = LocalVertsMap_P[tagV][2];
            Vfaces.push_back(Vf);
            
            VcF[0] = VcF[0] + Vf[0];
            VcF[1] = VcF[1] + Vf[1];
            VcF[2] = VcF[2] + Vf[2];
            
            if(tag2globvID_P.find(tagV)!=tag2globvID_P.end() &&
                unique_trace_verts2refmap.find(tagV)==unique_trace_verts2refmap.end())
            {
                int globid  = tag2globvID_P[tagV];
                fce[g]      = globid;
            
            }
            else if(SharedVertsNotOwned.find(tagV)!=SharedVertsNotOwned.end() &&
                    unique_trace_verts2refmap.find(tagV)==unique_trace_verts2refmap.end())
            {
                int globid  = SharedVertsNotOwned[tagV];
                fce[g]      = globid;   
            }
            else if(unique_trace_verts2refmap.find(tagV)!=unique_trace_verts2refmap.end())
            {
                int ref     = unique_trace_verts2refmap[tagV];
                int globid  = tracetagV2globalV[ref];
                fce[g]      = globid;
            }
        }

        VcF[0] = VcF[0]/npf;
        VcF[1] = VcF[1]/npf;
        VcF[2] = VcF[2]/npf;

        if(rh_P[gfid] == 0 || lh_P[gfid] == 0)
        {
            std::cout << world_rank << " rhp[gfid] and lhp[gfid] are zero SHARED " << std::endl;
            ishare++;
        }

        int leftEl  = lh_P[gfid];
        int leftTag = gE2tagE_P[leftEl];
        std::vector<double> Vijk(3,0);
        Vijk[0] = 0.0;
        Vijk[1] = 0.0;
        Vijk[2] = 0.0;
        // compute element center;
        int nvp = ElementTag2VertexTag_P[leftTag].size();

        for(int q=0;q<nvp;q++)
        {
            int tagV  = ElementTag2VertexTag_P[leftTag][q];

            Vijk[0] = Vijk[0] + LocalVertsMap_P[tagV][0];
            Vijk[1] = Vijk[1] + LocalVertsMap_P[tagV][1];
            Vijk[2] = Vijk[2] + LocalVertsMap_P[tagV][2];
        }

        Vijk[0] = Vijk[0]/nvp;
        Vijk[1] = Vijk[1]/nvp;
        Vijk[2] = Vijk[2]/nvp;

        double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);

        if(orient0 < 0.0)
        {
            if(npf == 3)
            {
                ifn_P[fptot*8+1] = fce[0];
                ifn_P[fptot*8+2] = fce[2];
                ifn_P[fptot*8+3] = fce[1];
                ifn_P[fptot*8+4] = 0;
            }
            if(npf == 4)
            {
                ifn_P[fptot*8+1] = fce[0];
                ifn_P[fptot*8+2] = fce[3];
                ifn_P[fptot*8+3] = fce[2];
                ifn_P[fptot*8+4] = fce[1];
            }
        }
        else
        {
            if(npf == 3)
            {
                ifn_P[fptot*8+1] = fce[0];
                ifn_P[fptot*8+2] = fce[1];
                ifn_P[fptot*8+3] = fce[2];
                ifn_P[fptot*8+4] = 0;
            }
            if(npf == 4)
            {
                ifn_P[fptot*8+1] = fce[0];
                ifn_P[fptot*8+2] = fce[1];
                ifn_P[fptot*8+3] = fce[2];
                ifn_P[fptot*8+4] = fce[3];
            }
        }

        for(int h=0;h<fce.size();h++)
        {
            if(fce[h]>162713)
            {
                std::cout << "Error line 303 Shared " << fce[h] << " " << 162713 << std::endl;
            }
        }


        
        ifn_P[fptot*8+5] = rh_P[gfid];
        ifn_P[fptot*8+6] = lh_P[gfid];
        ifn_P[fptot*8+7] = 2;

        // if(lh_P[gfid]>279070 || rh_P[gfid]>279070)
        // {
        //     std::cout << lh_P[gfid] << " " << rh_P[gfid] << " SharedFaces_P " << std::endl;
        // }


        if(ifn_P[fptot*8+0]==4 && ifn_P[fptot*8+1]==368974 && ifn_P[fptot*8+2]==369090 && ifn_P[fptot*8+3]==369091 &&
           ifn_P[fptot*8+4]==368975 && ifn_P[fptot*8+5]==279071 && ifn_P[fptot*8+6]==254992 && ifn_P[fptot*8+7]==2)
        {
            std::cout << "FOUND ERROR FACE " << std::endl;
        }

        fptot++;
    }

    std::cout << "ishare " << ishare << " " << SharedFaces_P.size() << " " << world_rank << std::endl;
}


void WriteBoundaryDataUS3DFormat(MPI_Comm comm, 
                                RepartitionObject* prism_repart,
                                PrismTetraTrace* pttrace,
                                std::vector<int> ifn_T,
                                std::vector<int> ifn_P,
                                std::map<int,std::vector<int> > bcref2bcface_T,
                                std::map<int,int> unique_trace_verts2refmap,
                                std::map<int,int> tracetagV2globalV,
                                std::map<int,std::vector<int> > BoundaryFaces_T,
                                std::map<int,int> glob2locVid_T,
                                std::vector<std::vector<int> > new_tetrahedra_T,
                                std::vector<double> new_vertices_T,
                                std::map<int,int> LocationSharedVert_T,
                                std::map<int,int> oldglob2newglob_T,
                                std::map<int,int> lh_T_bc,
                                std::map<int,int> new_globE2locE_T,
                                std::vector<std::vector<int> > &bcArrays,
                                std::map<int,int> &bcsizing,
                                std::vector<std::vector<int> > zdefs,
                                std::map<int,int> zone2bcref,
                                std::map<int,char*> zone2name,
                                std::vector<std::vector<char> > znames,
                                std::vector<int> &bciTot_offsets,
                                int &nTotBCFaces,
                                std::set<int> &bcsToT,
                                std::vector<int> &nlbc,
                                std::vector<int> &bcid,
                                std::vector<int> &bci_offsets,
                                std::vector<int> &zdefs_new)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    std::map<int,std::vector<int> > BoundaryFaces_P         = prism_repart->getBoundaryFaceMap();
    std::map<int, std::vector<double> > LocalVertsMap_P     = prism_repart->getLocalVertsMap();
    // std::map<int,int> tag2globvID_P                      = prism_repart->getUpdatedTag2GlobalVMap();
    std::map<int,int> SharedVertsNotOwned                   = prism_repart->getSharedVertsNotOwned();
    std::map<int,int> gE2tagE_P                             = prism_repart->getGlobalElement2ElementTag();
    std::map<int,int>  lh_P                                 = prism_repart->GetLeftHandFaceElementMap();
    std::map<int,int>  rh_P                                 = prism_repart->GetRightHandFaceElementMap();
    std::map<int,std::vector<int> > ElementTag2VertexTag_P  = prism_repart->getElement2VertexMap();

    std::map<int,int> tag2globvID_P;
    std::map<int,int> NonSharedVertsOwned_P                 = prism_repart->getNonSharedVertsOwned();
    std::map<int,int> SharedVertsOwned_P                    = prism_repart->getSharedVertsOwned();

    int nLocIntVrts         = NonSharedVertsOwned_P.size();
    int nLocShVrts          = SharedVertsOwned_P.size();
    DistributedParallelState* distnLocIntVrts = new DistributedParallelState(nLocIntVrts,comm);
    DistributedParallelState* distnLocShVrts  = new DistributedParallelState(nLocShVrts,comm);
    int vert = 0;
    std::map<int,int>::iterator itmii;
    for(itmii=NonSharedVertsOwned_P.begin();itmii!=NonSharedVertsOwned_P.end();itmii++)
    {
        int tag = itmii->first;
        tag2globvID_P[tag] = vert+distnLocIntVrts->getOffsets()[world_rank]+1;
        
        vert++;
    }


    int vertsh = 0;
    
    for(itmii=SharedVertsOwned_P.begin();itmii!=SharedVertsOwned_P.end();itmii++)
    {
        int tag                      = itmii->first;
        
        tag2globvID_P[tag] = SharedVertsOwned_P[tag];
    }

    // DistributedParallelState* distFaceTetraTot = new DistributedParallelState(ifn_T.size(),comm);
    // DistributedParallelState* distFacePrismTot = new DistributedParallelState(ifn_P.size(),comm);

    // int nTotInteriorFaces_tetra          = distFaceTetraTot->getNel();    
    // int nTotInteriorFaces_prism          = distFacePrismTot->getNel();
   
    std::map<int,std::vector<int> >::iterator bit;
    std::set<int> sorted_BCid;
    std::map<int,int> sorted_BCid_map;
    std::map<int,int> sorted_NBCid_map;
    int ii = 0;
    int Nbf;

    std::map<int,std::vector<int> > bcref2bcface_P = prism_repart->getZone2boundaryFaceID();
    int nloc_bcs_p = bcref2bcface_P.size();
    std::set<int> bcids_tot;
    
    std::vector<int> Lbcs;
    
    for(bit=bcref2bcface_T.begin();bit!=bcref2bcface_T.end();bit++)
    {
        if(bcids_tot.find(bit->first)==bcids_tot.end())
        {
            bcids_tot.insert(bit->first);
            Lbcs.push_back(bit->first);
        }
    }
    
    for(bit=bcref2bcface_P.begin();bit!=bcref2bcface_P.end();bit++)
    {
        std::cout << "bit->first " << bit->first << std::endl; 
        if(bcids_tot.find(bit->first)==bcids_tot.end())
        {
            bcids_tot.insert(bit->first);
            Lbcs.push_back(bit->first);
        }
    }
    
    int nloc_bcs   = Lbcs.size();

    DistributedParallelState* distLocBCs = new DistributedParallelState(nloc_bcs,comm);
     
    int Nt_BCs                 = distLocBCs->getNel();
    int* BCs_offsets           = distLocBCs->getOffsets();
    int* BCs_nlocs             = distLocBCs->getNlocs();
    int* BCs_arr               = new int[Nt_BCs];
    
    MPI_Allgatherv(&Lbcs[0],
                    nloc_bcs,
                    MPI_INT,
                    BCs_arr,
                    BCs_nlocs,
                    BCs_offsets,
                    MPI_INT, comm);

    
    std::map<int,int> bcentry;
    for(int i=0;i<Nt_BCs;i++)
    {
        if(bcsToT.find(BCs_arr[i])==bcsToT.end())
        {
            bcsToT.insert(BCs_arr[i]);
        }
    }
    
    // int* bcid = new int[bcsToT.size()];
    // int* nlbc = new int[bcsToT.size()];
    bcid.resize(bcsToT.size());
    std::fill(bcid.begin(), bcid.end(), 0);
    nlbc.resize(bcsToT.size());
    std::fill(nlbc.begin(), nlbc.end(), 0);

    std::set<int>::iterator entry;
    int cnt = 0;

    int lac;
    int q = 0;
    for(entry=bcsToT.begin();entry!=bcsToT.end();entry++)
    {
        int ee   = *entry;
        int Nbft = 0;
        int Nbfp = 0;
        
        if(bcref2bcface_T.find(ee)!=bcref2bcface_T.end())
        {
            Nbft = bcref2bcface_T[ee].size();
        }
        if(bcref2bcface_P.find(ee)!=bcref2bcface_P.end())
        {
            Nbfp = bcref2bcface_P[ee].size();
        }

        bcid[q]  = ee;
        nlbc[q]  = Nbft+Nbfp;
        q++;
    }



    
    nTotBCFaces = 0;
    int nTotBCFaces_offset = 0;

    int failbc = 0;
    int globalVid;
    std::vector<double> VcF(3,0);
    std::vector<double> Vijk(3,0);
    std::vector<int> fq(3);
    int tel = 0;
    
    std::cout << "zone2bcref " << zone2bcref.size() << " bcsToT " << bcsToT.size() << std::endl; 
    for(int i=0;i<bcsToT.size();i++)
    {
        int bc_id = bcid[i];
        std::cout << "bc_id " << bc_id << std::endl;
        int bc_reference=-1;
        
        if(zone2bcref.find(bc_id)!=zone2bcref.end())
        {
            bc_reference = zone2bcref[bc_id];
        }
    
        
        DistributedParallelState* distBCi = new DistributedParallelState(nlbc[i],comm);

        int NelLoc_bci = nlbc[i];
        int NelTot_bci = distBCi->getNel();
        std::vector<int> ifn_bc(NelLoc_bci*8,0);
        int offsetbci = distBCi->getOffsets()[world_rank];
        int fbc  = 0;
        int Nbft = 0;
        int Nbfp = 0;
        
        if(bcref2bcface_T.find(bc_id)!=bcref2bcface_T.end())
        {
            Nbft = bcref2bcface_T[bc_id].size();
        }
        if(bcref2bcface_P.find(bc_id)!=bcref2bcface_P.end())
        {
            Nbfp = bcref2bcface_P[bc_id].size();
        }


        if(Nbfp!=0)
        {
            int sk = 0;
            
            for(int q=0;q<Nbfp;q++)
            {
                int bcface = bcref2bcface_P[bc_id][q];
                int flag = -1;
                
                int nppf = BoundaryFaces_P[bcface].size();
                
                ifn_bc[fbc*8+0] = nppf;
                std::vector<int> face_tmp(nppf);
                std::vector<int> fce(nppf,0);
                std::vector<std::vector<double> > Vfaces;
                
                VcF[0] = 0.0;
                VcF[1] = 0.0;
                VcF[2] = 0.0;
                
                for(int g=0;g<nppf;g++)
                {
                    int tagV = BoundaryFaces_P[bcface][g];

                    std::vector<double> Vf(3,0.0);
                    Vf[0] = LocalVertsMap_P[tagV][0];
                    Vf[1] = LocalVertsMap_P[tagV][1];
                    Vf[2] = LocalVertsMap_P[tagV][2];
                    Vfaces.push_back(Vf);

                    VcF[0] = VcF[0] + Vf[0];
                    VcF[1] = VcF[1] + Vf[1];
                    VcF[2] = VcF[2] + Vf[2];
                    
                    if(tag2globvID_P.find(tagV)!=tag2globvID_P.end() &&
                        unique_trace_verts2refmap.find(tagV)==unique_trace_verts2refmap.end())
                    {
                        int globid = tag2globvID_P[tagV];
                        fce[g]     = globid;
                    }
                    else if(SharedVertsNotOwned.find(tagV)!=SharedVertsNotOwned.end() &&
                            unique_trace_verts2refmap.find(tagV)==unique_trace_verts2refmap.end())
                    {
                        int globid = SharedVertsNotOwned[tagV];
                        fce[g]     = globid;
                    }
                    else if(unique_trace_verts2refmap.find(tagV)!=unique_trace_verts2refmap.end())
                    {
                        int ref     = unique_trace_verts2refmap[tagV];
                        int globid  = tracetagV2globalV[ref];
                        fce[g]      = globid;
                    }
                    // std::cout << " (" << tagV << " " << fce[g] << ") ";
                }
                
                // std::cout << std::endl;
                VcF[0] = VcF[0]/nppf;
                VcF[1] = VcF[1]/nppf;
                VcF[2] = VcF[2]/nppf;
                
                int leftEl  = lh_P[bcface];
                int leftTag = gE2tagE_P[leftEl];
                
                Vijk[0] = 0.0;
                Vijk[1] = 0.0;
                Vijk[2] = 0.0;
                // compute element center;
                int nvp = ElementTag2VertexTag_P[leftTag].size();

                for(int q=0;q<nvp;q++)
                {
                    int tagV  = ElementTag2VertexTag_P[leftTag][q];

                    Vijk[0] = Vijk[0] + LocalVertsMap_P[tagV][0];
                    Vijk[1] = Vijk[1] + LocalVertsMap_P[tagV][1];
                    Vijk[2] = Vijk[2] + LocalVertsMap_P[tagV][2];
                }

                Vijk[0] = Vijk[0]/nvp;
                Vijk[1] = Vijk[1]/nvp;
                Vijk[2] = Vijk[2]/nvp;

                double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);
                
                if(orient0 < 0.0)
                {
                    if(nppf == 3)
                    {
                        ifn_bc[fbc*8+1] = fce[0];
                        ifn_bc[fbc*8+2] = fce[2];
                        ifn_bc[fbc*8+3] = fce[1];
                        ifn_bc[fbc*8+4] = 0;
                    }
                    if(nppf == 4)
                    {
                        ifn_bc[fbc*8+1] = fce[0];
                        ifn_bc[fbc*8+2] = fce[3];
                        ifn_bc[fbc*8+3] = fce[2];
                        ifn_bc[fbc*8+4] = fce[1];
                    }
                    
                }
                else
                {
                    if(nppf == 3)
                    {
                        ifn_bc[fbc*8+1] = fce[0];
                        ifn_bc[fbc*8+2] = fce[1];
                        ifn_bc[fbc*8+3] = fce[2];
                        ifn_bc[fbc*8+4] = 0;
                    }
                    if(nppf == 4)
                    {
                        ifn_bc[fbc*8+1] = fce[0];
                        ifn_bc[fbc*8+2] = fce[1];
                        ifn_bc[fbc*8+3] = fce[2];
                        ifn_bc[fbc*8+4] = fce[3];
                    }
                }

                for(int h=0;h<fce.size();h++)
                {
                    if(fce[h]>162713)
                    {
                        std::cout << "Error line 643 Boundary Prism " << fce[h] << " " << 162713 << std::endl;
                    }
                }
                
                ifn_bc[fbc*8+5] = 0;
                ifn_bc[fbc*8+6] = lh_P[bcface];
                ifn_bc[fbc*8+7] = bc_reference;

                fbc++;
            }   
        }
        
        if(Nbft!=0 )
        {

            for(int q=0;q<Nbft;q++)
            {
                int bcface      = bcref2bcface_T[bc_id][q];
                ifn_bc[fbc*8+0] = 3;
                std::vector<std::vector<double> > Vfaces;
                
                VcF[0] = 0.0;
                VcF[1] = 0.0;
                VcF[2] = 0.0;

                for(int w=0;w<3;w++)
                {
                    int vertexID = BoundaryFaces_T[bcface][w];
                    int lvert = glob2locVid_T[vertexID];
                    
                    std::vector<double> Vf(3);
                    Vf[0]  = new_vertices_T[(lvert-1)*3+0];
                    Vf[1]  = new_vertices_T[(lvert-1)*3+1];
                    Vf[2]  = new_vertices_T[(lvert-1)*3+2];
                    
                    VcF[0] = VcF[0] + Vf[0];
                    VcF[1] = VcF[1] + Vf[1];
                    VcF[2] = VcF[2] + Vf[2];
                    
                    Vfaces.push_back(Vf);
                    
                    if(LocationSharedVert_T.find(vertexID)!=LocationSharedVert_T.end())
                    {
                        fq[w] = LocationSharedVert_T[vertexID];
                        ifn_bc[fbc*8+(w+1)] = LocationSharedVert_T[vertexID];
                        
                    }
                    else
                    {
                        fq[w] = oldglob2newglob_T[vertexID];
                        ifn_bc[fbc*8+(w+1)] = oldglob2newglob_T[vertexID];
                        
                    }
                }
                
                //std::cout << std::endl;
                int elLh    = lh_T_bc[bcface];
                int locElID = new_globE2locE_T[elLh];
                Vijk[0] = 0.0;
                Vijk[1] = 0.0;
                Vijk[2] = 0.0;
                // compute element center;
                //ienOUT[curElID] = Elvrts
                for(int u=0;u<new_tetrahedra_T[locElID].size();u++)
                {   
                    // Vijk[0] = Vijk[0] + new_vertices_T[new_tetrahedra_T[locElID][u]-1][0];
                    // Vijk[1] = Vijk[1] + new_vertices_T[new_tetrahedra_T[locElID][u]-1][1];
                    // Vijk[2] = Vijk[2] + new_vertices_T[new_tetrahedra_T[locElID][u]-1][2];

                    Vijk[0] = Vijk[0] + new_vertices_T[(new_tetrahedra_T[locElID][u]-1)*3+0];
                    Vijk[1] = Vijk[1] + new_vertices_T[(new_tetrahedra_T[locElID][u]-1)*3+1];
                    Vijk[2] = Vijk[2] + new_vertices_T[(new_tetrahedra_T[locElID][u]-1)*3+2];
                }

                Vijk[0] = Vijk[0]/new_tetrahedra_T[locElID].size();
                Vijk[1] = Vijk[1]/new_tetrahedra_T[locElID].size();
                Vijk[2] = Vijk[2]/new_tetrahedra_T[locElID].size();
                
                double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);

                if(orient0 < 0.0)
                {
                    std::cout << "Weve got negative faces " << orient0 << " bc_reference " << bc_reference << " " << elLh << " " << tel << std::endl;
                    std::cout << "VcF " << VcF[0] << " " << VcF[1] << " " << VcF[2] << std::endl;
                    for(int i=0;i<Vfaces.size();i++)
                    {
                        std::cout << Vfaces[i][0] << " " << Vfaces[i][1] << " " << Vfaces[i][2] << std::endl;
                    }
                    std::cout << "======RANK=====" << world_rank << std::endl;
                    tel = tel + 1;
                }
                
                Vfaces.clear();
                //std::cout << "bc_reference "<< bc_reference << std::endl;
                ifn_bc[fbc*8+4] = 0;
                ifn_bc[fbc*8+5] = 0;
                ifn_bc[fbc*8+6] = lh_T_bc[bcface];
                ifn_bc[fbc*8+7] = bc_reference;
                
                for(int h=0;h<fq.size();h++)
                {
                    if(fq[h]>162713)
                    {
                        std::cout << "Error line 739 Boundary Tet " << fq[h] << " " << 162713 << std::endl;
                    }
                }


                fbc++;
                
            }
            
        }

        int ftott = (int)ifn_T.size()/8;
        int ftotp = (int)ifn_P.size()/8;

        int nbt = Nbft+Nbfp;
        
        bcsizing[bc_id] = NelTot_bci;
        bci_offsets.push_back(offsetbci);
        bciTot_offsets.push_back(nTotBCFaces_offset);
        bcArrays.push_back(ifn_bc);
        
        nTotBCFaces_offset = nTotBCFaces_offset + NelTot_bci;
        nTotBCFaces        = nTotBCFaces + NelTot_bci;   
        
    }
    int ftott = (int)ifn_T.size()/8;
    int ftotp = (int)ifn_P.size()/8;
    int nNonSharedVertsOwned_P                  = prism_repart->getNonSharedVertsOwned().size();
    int nSharedVertsOwned_P                     = prism_repart->getSharedVertsOwned().size();

    int nvertices = (int)new_vertices_T.size()/3;

    std::cout << "nvertices " << nvertices << std::endl;

    
    DistributedParallelState* distTetraVerts    = new DistributedParallelState(nvertices,comm);
    DistributedParallelState* distPrismVerts    = new DistributedParallelState(nNonSharedVertsOwned_P+nSharedVertsOwned_P,comm);
    
    DistributedParallelState* distTetra         = new DistributedParallelState(new_tetrahedra_T.size(),comm);
    DistributedParallelState* distPrism         = new DistributedParallelState(prism_repart->getLocElem().size(),comm);

    DistributedParallelState* distftott         = new DistributedParallelState(ftott,comm);
    DistributedParallelState* distftotp         = new DistributedParallelState(ftotp,comm);

    int ToTElements_prism                       = distPrism->getNel();
    int ToTElements                             = distTetra->getNel();
    int nTotElements                            = ToTElements_prism+ToTElements;

    int nTotTetraVerts                          = distTetraVerts->getNel();
    int nTotPrismVerts                          = distPrismVerts->getNel();
    int nTotVertsPrismTetra                     = nTotPrismVerts+nTotTetraVerts;

    int nTotInteriorFaces_prism                 = distftotp->getNel();
    int nTotInteriorFaces_tetra                 = distftott->getNel();

    int nTotIntFaces                            = nTotInteriorFaces_prism + nTotInteriorFaces_tetra;

    std::cout << " nTotInteriorFaces_prism " << nTotInteriorFaces_prism << std::endl;
    std::cout << " nTotInteriorFaces_tetra " << nTotInteriorFaces_tetra << std::endl;
    int nbo = bcArrays.size();
    //std::cout << "-- Constructing the zdefs array..."<<std::endl;
    // Array<int>* adapt_zdefs = new Array<int>(3+nbo,7);

    // Collect node data (10) . Starting index-ending index Nodes
    // std::vector<std::vector<int> > zdefs_new(3+nbo);
    std::cout << "3+nbo " << 3+nbo << " " << nbo << std::endl;
    std::vector<std::vector<int> > zdefs_new_copy;
    std::vector<int> zdefs_row0(3+nbo,0);
    zdefs_row0[0] = 10;
    zdefs_row0[1] = -1;
    zdefs_row0[2] = 1;
    zdefs_row0[3] = 1;
    zdefs_row0[4] = nTotVertsPrismTetra;
    zdefs_row0[5] = zdefs[0][5];
    zdefs_row0[6] = zdefs[0][6];
    for(int i = 0;i<(3+nbo);i++)
    {
        zdefs_new.push_back(zdefs_row0[i]);
    }
    zdefs_new_copy.push_back(zdefs_row0);
    // Collect element data (12) . Starting index-ending index Element
    std::vector<int> zdefs_row1(3+nbo,0);
    zdefs_row1[0] = 12;
    zdefs_row1[1] = -1;
    zdefs_row1[2] = 2;
    zdefs_row1[3] = 1;
    zdefs_row1[4] = nTotElements;
    zdefs_row1[5] = zdefs[1][5];
    zdefs_row1[6] = 2;
    for(int i = 0;i<(3+nbo);i++)
    {
        zdefs_new.push_back(zdefs_row1[i]);
    }
    zdefs_new_copy.push_back(zdefs_row1);
    // Collect internal face data (13) . Starting index-ending index internal face.
    std::vector<int> zdefs_row2(3+nbo,0);
    zdefs_row2[0] = 13;
    zdefs_row2[1] = -1;
    zdefs_row2[2] = 3;
    zdefs_row2[3] = 1;
    zdefs_row2[4] = nTotIntFaces;
    zdefs_row2[5] = zdefs[2][5];
    zdefs_row2[6] = 2;
    for(int i = 0;i<(3+nbo);i++)
    {
        zdefs_new.push_back(zdefs_row2[i]);
    }
    zdefs_new_copy.push_back(zdefs_row2);

    int qq  = 1;
    int nb  = 0;
    int face_start = nTotIntFaces+1;
    int face_end;
    int bnd_zone;
    std::map<int,int>::iterator itr;
    std::map<int,char*> newzone_2_name;
    std::map<int,int> newzone2oldzone;
    for(itr=bcsizing.begin();itr!=bcsizing.end();itr++)
    {
        int bnd_zone  		  = itr->first;
        int bnd_ref   		  = bnd_zone;
        newzone2oldzone[nb+3] = bnd_zone;
        
        if(zone2bcref.find(bnd_zone)!=zone2bcref.end())
        {
            bnd_ref = zone2bcref[bnd_ref];
            newzone_2_name[nb+3]=zone2name[bnd_zone];
        }
        
        int bnd_size = itr->second;
        face_end = face_start+bnd_size-1;
        std::vector<int> zdefs_row(3+nbo,0);
        zdefs_row[0] = 13;
        zdefs_row[1] = -1;
        zdefs_row[2] = 3+qq;
        zdefs_row[3] = face_start;
        zdefs_row[4] = face_end;
        zdefs_row[5] = bnd_ref;
        zdefs_row[6] = 2;
        for(int i = 0;i<(3+nbo);i++)
        {
            zdefs_new.push_back(zdefs_row[i]);
        }
        face_start = face_end+1;
        zdefs_new_copy.push_back(zdefs_row);
        
        nb++;
        qq++;
    }            
    
    if(world_rank == 0)
    {
        PlotBoundaryData_Lite(znames,zdefs_new_copy);
    }
}