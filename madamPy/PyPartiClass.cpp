#include <Python.h>
#include <mpi.h>
#include <mpi4py/mpi4py.h>
#define MPICH_SKIP_MPICXX 1
#define OMPI_SKIP_MPICXX  1
#include <numpy/arrayobject.h>
#include <CGAL/Simple_cartesian.h>
#include <cstdio>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <unordered_map>

//~ #include "structmember.h"

#include "../src/adapt_io.h"
#include "../src/adapt_partition.h"
#include "../src/adapt_recongrad.h"
//using transformation::Parti;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;



typedef struct {
    PyObject_HEAD
    MPI_Comm commu;
    Partition * ptrObj;
    Array<double>* gB;
} PyPartition;










template<typename T>
typename std::map<int,T > getMapFromPyDict(PyObject* pyDict)
{
    PyObject* pKeys       = PyDict_Keys(pyDict);
    PyObject* pValues     = PyDict_Values(pyDict);
    int len               = (int)PyList_Size(pValues);
    PyObject* firstItem   = PyList_GET_ITEM(pValues, 0);

    //PyDict_Next(index, &pos, &key, &value);
    
    std::map<int,T > mapCpp;
    
    if (PyLong_Check(firstItem))
    {
        for (Py_ssize_t i = 0; i < PyDict_Size(pyDict); ++i)
        {
            int key = (int) PyLong_AsLong(PyList_GetItem(pKeys,   i));
            int val = (int) PyLong_AsLong(PyList_GetItem(pValues,   i));

            mapCpp.insert( {key, val} );
        }
    }
    if (PyFloat_Check(firstItem))
    {
        for (Py_ssize_t i = 0; i < PyDict_Size(pyDict); ++i)
        {
            int key = (int) PyLong_AsLong(PyList_GetItem(pKeys,   i));
            double val = (double) PyFloat_AsDouble(PyList_GetItem(pValues,   i));

            mapCpp.insert( {key, val} );
        }
    }
    
    
    return mapCpp;
    
}










template<typename T>
typename std::map<int,Array<T>* > getMapOfArraysFromPyDict(PyObject* pyDict)
{

    std::map<int,Array<T>* > output;


    PyObject* pKeys       = PyDict_Keys(pyDict);
    PyObject* pValues     = PyDict_Values(pyDict);
    int n                 = (int)PyList_Size(pValues);
    PyObject* firstItem   = PyList_GET_ITEM(pValues, 0);

    int ncol;
    T row_val;
    if (PyLong_Check(firstItem))
    {
        for (int i = 0; i < n; i++)
        {
            int key = (int) PyLong_AsLong(PyList_GetItem(pKeys,   i));
            PyObject* value = PyList_GetItem(pValues, i);
            Array<T>* row = new Array<T>(1,1);
            int val_entry = (int) PyLong_AsLong(value);
            row->setVal(0,0,val_entry);
            output[key] = row;
        }
    }
//                                                                                                                            
    if (PyFloat_Check(firstItem))
    {
        for (int i = 0; i < n; i++)
        {
            int key = (int) PyLong_AsLong(PyList_GetItem(pKeys,   i));
            PyObject* value = PyList_GetItem(pValues, i);
            Array<T>* row = new Array<T>(1,1);
            double val_entry = (double) PyFloat_AsDouble(value);
            row->setVal(0,0,val_entry);
            output[key] = row;
            //std::cout << "key val = " << key << " " << val_entry << std::endl;                                              
        }
    }

    return output;
}















static int PyPartition_init(PyPartition* self, PyObject *args, PyObject *kwds)
// initialize PyParti Object
{
    int argc = PyTuple_GET_SIZE(args);
    
    std::cout << "Initializing the parititions " << std::endl;
  
    PyObject *py_comm		=	NULL;
    PyObject *py_comm_i		=	NULL;

    MPI_Comm *comm_p		=	NULL;
    MPI_Info *comm_p_info	=	NULL;
    
    const char* gname;
    const char* cname;
    const char* dname;
    const char* vname;
    int ReadFromStats;
    int inti;
    PyObject *py_dict;
    if (! PyArg_ParseTuple(args, "sssisOO", &cname, &gname, &dname, &ReadFromStats, &vname, &py_comm, &py_comm_i))
        return -1;
    
    comm_p 		 = PyMPIComm_Get(py_comm);
    comm_p_info  = PyMPIInfo_Get(py_comm_i);
    int world_size;
    MPI_Comm_size(*comm_p, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(*comm_p, &world_rank);
    self->commu = *comm_p;
    
    US3D* us3d   = ReadUS3DData(cname,gname,dname,ReadFromStats,*comm_p,*comm_p_info);
    
    int Nve      = us3d->xcn->getNglob();
    int Nel_part = us3d->ien->getNrow();
    
    ParallelState* ien_pstate               = new ParallelState(us3d->ien->getNglob(),*comm_p);
    ParallelState* ife_pstate               = new ParallelState(us3d->ifn->getNglob(),*comm_p);
    ParallelState_Parmetis* parmetis_pstate = new ParallelState_Parmetis(us3d->ien,us3d->elTypes,us3d->ie_Nv,*comm_p);
    ParallelState* xcn_pstate               = new ParallelState(us3d->xcn->getNglob(),*comm_p);

    Array<double>* Uivar = new Array<double>(Nel_part,1);
    double rhoState,uState,vState,wState,TState,VtotState,aState,MState;
    Array<double>* gB = new Array<double>(us3d->ghost->getNrow(),1);
//
//
    if (strcmp(vname, "rho") == 0)
    {
        for(int i=0;i<Nel_part;i++)
        {
          rhoState = us3d->interior->getVal(i,0);
          Uivar->setVal(i,0,rhoState);
        }
        for(int i=0;i<us3d->ghost->getNrow();i++)
        {
          gB->setVal(i,0,us3d->ghost->getVal(i,0));
        }
    }
    if (strcmp(vname, "u") == 0)
    {
        for(int i=0;i<Nel_part;i++)
        {
          uState = us3d->interior->getVal(i,1);
          Uivar->setVal(i,0,uState);
        }
        for(int i=0;i<us3d->ghost->getNrow();i++)
        {
          gB->setVal(i,0,us3d->ghost->getVal(i,1));
        }
    }
    if (strcmp(vname, "v") == 0)
    {
      for(int i=0;i<Nel_part;i++)
      {
          uState = us3d->interior->getVal(i,2);
          Uivar->setVal(i,0,uState);
      }
      for(int i=0;i<us3d->ghost->getNrow();i++)
      {
          gB->setVal(i,0,us3d->ghost->getVal(i,2));
      }
    }
    if (strcmp(vname, "w") == 0)
    {
      for(int i=0;i<Nel_part;i++)
      {
          uState = us3d->interior->getVal(i,3);
          Uivar->setVal(i,0,uState);
      }
      for(int i=0;i<us3d->ghost->getNrow();i++)
      {
          gB->setVal(i,0,us3d->ghost->getVal(i,3));
      }
    }
    if (strcmp(vname, "T") == 0)
    {
      for(int i=0;i<Nel_part;i++)
      {
          uState = us3d->interior->getVal(i,4);
          Uivar->setVal(i,0,uState);
      }
      for(int i=0;i<us3d->ghost->getNrow();i++)
      {
          gB->setVal(i,0,us3d->ghost->getVal(i,4));
      }
    }
    if (strcmp(vname, "Vt") == 0)
    {
      for(int i=0;i<Nel_part;i++)
      {
          uState   = us3d->interior->getVal(i,1);
          vState   = us3d->interior->getVal(i,2);
          wState   = us3d->interior->getVal(i,3);
          VtotState = sqrt(uState*uState+vState*vState+wState*wState);
          Uivar->setVal(i,0,VtotState);
      }
      for(int i=0;i<us3d->ghost->getNrow();i++)
      {
          uState   = us3d->ghost->getVal(i,1);
          vState   = us3d->ghost->getVal(i,2);
          wState   = us3d->ghost->getVal(i,3);
          VtotState = sqrt(uState*uState+vState*vState+wState*wState);
          gB->setVal(i,0,VtotState);
      }
    }
    if (strcmp(vname, "Mach") == 0)
    {
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
          Uivar->setVal(i,0,MState);
        }

        for(int i=0;i<us3d->ghost->getNrow();i++)
        {
            rhoState = us3d->ghost->getVal(i,0);
            uState   = us3d->ghost->getVal(i,1);
            vState   = us3d->ghost->getVal(i,2);
            wState   = us3d->ghost->getVal(i,3);
            TState   = us3d->ghost->getVal(i,4);
            VtotState = sqrt(uState*uState+vState*vState+wState*wState);
            aState   = sqrt(1.4*287.05*TState);
            MState = VtotState/aState;
            gB->setVal(i,0,MState);
        }
    }
//
    
    
    PyObject *py_ghost = PyList_New(us3d->ghost->getNrow());

    for (int q = 0;q<us3d->ghost->getNrow();q++)
    {
        PyObject *entryGhost = PyFloat_FromDouble(gB->getVal(q,0));

        PyList_SetItem(py_ghost,q,entryGhost);
    }
    
    PyObject *part_py = PyList_New(13);

    PyList_SET_ITEM(part_py, 12, py_ghost);
    self->gB = gB;
    delete us3d->interior;

    self->ptrObj = new Partition(us3d->ien, us3d->iee, us3d->ief, us3d->ie_Nv , us3d->ie_Nf,
                                 us3d->ifn, us3d->ife, us3d->if_ref, us3d->if_Nv,
                                 parmetis_pstate, ien_pstate, ife_pstate,
                                 us3d->xcn, xcn_pstate, Uivar, *comm_p);

    
    
    
    
    return 0;
}



static void PyPartition_dealloc(PyPartition * self)
// destruct the object
{
    delete self->ptrObj;
    Py_TYPE(self)->tp_free(self);
}










static PyObject * PyPartition_getVertices(PyPartition* self, PyObject* args)
{

    std::vector<Vert*> verts = (self->ptrObj)->getLocalVerts();
    int nvert = verts.size();
    PyObject *pList_verts = PyList_New((nvert));
    
    for(int i=0;i<nvert;i++)
    {
        PyObject *coords = PyTuple_New(3);
        PyTuple_SET_ITEM(coords, 0, PyFloat_FromDouble(verts[i]->x));
        PyTuple_SET_ITEM(coords, 1, PyFloat_FromDouble(verts[i]->y));
        PyTuple_SET_ITEM(coords, 2, PyFloat_FromDouble(verts[i]->z));
        PyList_SetItem(pList_verts,i,coords);

    }
    
    return pList_verts;
}









static PyObject * PyPartition_getUState(PyPartition* self, PyObject* args)
{
    std::vector<int> LocElem = (self->ptrObj)->getLocElem();
    std::vector<double> Ustate = (self->ptrObj)->getLocElemVaria();

    PyObject *pDict_UState = PyDict_New();
    for (int q = 0;q<LocElem.size();q++)
    {
        PyObject *key = PyLong_FromSsize_t(LocElem[q]);
	PyObject *val = PyFloat_FromDouble(Ustate[q]);

	PyDict_SetItem(pDict_UState,key,val);
        
    }
    
    return pDict_UState;
}










static PyObject * PyPartition_AddAdjUState(PyPartition* self, PyObject* args)
{
    PyObject* Umap;
    if (! PyArg_ParseTuple(args, "O", &Umap))
        return Py_False;

    std::map<int,Array<double>* > Uvaria_map = getMapOfArraysFromPyDict<double>(Umap);
    
    // convert dictionary and other inputs to c++
    self->ptrObj->AddStateVecForAdjacentElements(Uvaria_map,1,self->commu);

    PyObject *pDict_UState = PyDict_New();
    std::map<int,Array<double>* >::iterator it;
    for (it=Uvaria_map.begin();it!=Uvaria_map.end();it++)
    {
        PyObject *key = PyLong_FromSsize_t(it->first);
        PyObject *val = PyFloat_FromDouble(it->second->getVal(0,0));

        PyDict_SetItem(pDict_UState,key,val);

    }
   
    return pDict_UState;
}









static PyObject * PyPartition_computeGradU(PyPartition* self, PyObject* args)
{

    PyObject* Umap;
    if (! PyArg_ParseTuple(args, "O", &Umap))
        return Py_False;
    
    //Convert dictionary into a C++ map of arrays.

    std::map<int,Array<double>* > U_map = getMapOfArraysFromPyDict<double>(Umap);

    
    //Collect required datastructures to compute gradient.
    
    std::vector<int> LocElem              = (self->ptrObj)->getLocElem();
    std::vector<double> Ustate            = (self->ptrObj)->getLocElemVaria();
    std::vector<Vert*> verts              = (self->ptrObj)->getLocalVerts();
    std::map<int,int> gV2lV               = (self->ptrObj)->getGlobalVert2LocalVert();
    std::map<int,int> LocElem2Nf          = (self->ptrObj)->getLocElem2Nf();
    std::map<int,int> LocElem2Nv          = (self->ptrObj)->getLocElem2Nv();
    std::map<int,std::vector<int> > gE2lV = (self->ptrObj)->getGlobElem2LocVerts();
    i_part_map*  if_Nv_map                = (self->ptrObj)->getIF_Nvpartmap();
    i_part_map*  ifn_map                  = (self->ptrObj)->getIFNpartmap();
    i_part_map*  ief_map                  = (self->ptrObj)->getIEFpartmap();
    i_part_map*  iee_map                  = (self->ptrObj)->getIEEpartmap();
    int nGlob                             = (self->ptrObj)->getNglob_Elem();
    
    
    //Compute the gradient.
    
    
    std::map<int,Array<double>* > Ugrad = Py_ComputedUdx_LSQ_US3D(verts,
                                            gE2lV,gV2lV, LocElem,
                                            ifn_map->i_map,ief_map->i_map,
                                            iee_map->i_map,if_Nv_map->i_map,
                                            LocElem2Nf,LocElem2Nv,nGlob,
                                            U_map, self->gB, self->commu);
    
    
    // Convert the C++ map to a dictionary.
    
    PyObject *pDict_UState = PyDict_New();
    for (int q = 0;q<LocElem.size();q++)
    {
        PyObject *key = PyLong_FromSsize_t(LocElem[q]);
        PyObject *val = PyFloat_FromDouble(Ustate[q]);
        PyDict_SetItem(pDict_UState,key,val);
    }
    
    return pDict_UState;
}








static PyObject* PyPartition_gatherMeshOnRoot(PyPartition* self, PyObject* args)
{
    
    PyObject* Umap;
    if (! PyArg_ParseTuple(args, "O", &Umap))
        return Py_False;
    
    
    int size;
    MPI_Comm_size(self->commu, &size);
    int rank;
    MPI_Comm_rank(self->commu, &rank);
    
    std::vector<int> LocElem              = (self->ptrObj)->getLocElem();
    int nLoc                              = LocElem.size();
    std::vector<double> Ustate            = (self->ptrObj)->getLocElemVaria();
    std::vector<Vert*> verts              = (self->ptrObj)->getLocalVerts();
    std::map<int,int> gV2lV               = (self->ptrObj)->getGlobalVert2LocalVert();
    std::map<int,int> LocElem2Nf          = (self->ptrObj)->getLocElem2Nf();
    std::map<int,int> LocElem2Nv          = (self->ptrObj)->getLocElem2Nv();
    std::map<int,std::vector<int> > gE2lV = (self->ptrObj)->getGlobElem2LocVerts();
    i_part_map*  if_Nv_map                = (self->ptrObj)->getIF_Nvpartmap();
    i_part_map*  ifn_map                  = (self->ptrObj)->getIFNpartmap();
    i_part_map*  ief_map                  = (self->ptrObj)->getIEFpartmap();
    i_part_map*  iee_map                  = (self->ptrObj)->getIEEpartmap();
    int nElemGlob                         = (self->ptrObj)->getNglob_Elem();
    int nVertGlob                         = (self->ptrObj)->getNglob_Vert();

    
    int* proc_vnlocs  = new int[size];
    int* proc_voffset = new int[size];
    int* vnlocs       = new int[size];
    int* voffsets     = new int[size];

    int* proc_nlocs  = new int[size];
    int* proc_offset = new int[size];
    int* nlocs       = new int[size];
    int* offsets     = new int[size];
    std::cout << "rank := " << rank << " nv := " << verts.size() << " " << gV2lV.size() << std::endl;
    for(int i=0;i<size;i++)
    {
        nlocs[i]    = 0;
        offsets[i]  = 0;
        
        vnlocs[i]   = 0;
        voffsets[i] = 0;
        
        if(i==rank)
        {
            proc_nlocs[i]  = LocElem.size();
            proc_vnlocs[i] = verts.size();
        }
        else
        {
            proc_nlocs[i]  = 0;
            proc_vnlocs[i] = 0;

        }
    }
    MPI_Allreduce(proc_nlocs,  nlocs,     size, MPI_INT, MPI_SUM, self->commu);
    //MPI_Allreduce(proc_offset, offsets,   size, MPI_INT, MPI_SUM, self->commu);
    
    MPI_Allreduce(proc_vnlocs,  vnlocs,   size, MPI_INT, MPI_SUM, self->commu);
    //MPI_Allreduce(proc_voffset, voffsets, size, MPI_INT, MPI_SUM, self->commu);
    
    int offset = 0;
    int voffset = 0;
    for(int i=0;i<size;i++)
    {
        offsets[i] = offset;
        offset     = offset+nlocs[i];
        
        voffsets[i] = voffset;
        voffset     = voffset+vnlocs[i];
    }
    
    Array<int>*  lE2gE_g;
    Array<int>*  lV2gV_g;

    if(rank == 0)
    {
        lE2gE_g  = new Array<int>(nElemGlob,1);
        lV2gV_g  = new Array<int>(voffset,1);
    }
    else
    {
        lE2gE_g  = new Array<int>(1,1);
    }
    
    MPI_Gatherv(&LocElem[0],
                LocElem.size(),
                MPI_INT,
                &lE2gE_g->data[0],
                nlocs,
                offsets,
                MPI_INT, 0, self->commu);
    
    std::map<int,int>::iterator itmm;
    std::vector<int> lv2gv_vec(gV2lV.size());
    //gV2LV should be just a vector.
    for(itmm=gV2lV.begin();itmm!=gV2lV.end();itmm++)
    {
        lv2gv_vec[itmm->first]=itmm->second;
    }
    
    MPI_Gatherv(&lv2gv_vec[0],
                lv2gv_vec.size(),
                MPI_INT,
                &lV2gV_g->data[0],
                vnlocs,
                voffsets,
                MPI_INT, 0, self->commu);
    
    
    PyObject *pDict_UState = PyDict_New();

    
    return pDict_UState;
}

















static PyModuleDef PartitionModule = {
    PyModuleDef_HEAD_INIT,
    "PartiClass",
    "Example module that wrapped a C++ object",
    -1,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL
};











static PyMethodDef PyPartition_methods[] = {
    { "getVertices", (PyCFunction)PyPartition_getVertices,    METH_VARARGS,       "Get Ustate based on partitioning" },
  { "getUState", (PyCFunction)PyPartition_getUState,    METH_VARARGS,       "Get Ustate based on partitioning" },
    { "addAdjUState", (PyCFunction)PyPartition_AddAdjUState,    METH_VARARGS,       "Add State Solutions For Adjacent Elements" },
    { "computeGradU", (PyCFunction)PyPartition_computeGradU,    METH_VARARGS,       "Get Ustate based on partitioning" },
    {"gatherMeshOnRoot",   (PyCFunction)PyPartition_gatherMeshOnRoot,    METH_VARARGS,       "Get Ustate based on partitioning" },
    {NULL}  /* Sentinel */
};








static PyTypeObject PyPartitionType = { PyVarObject_HEAD_INIT(NULL, 0)
                                    "PartiClass.Partition"   /* tp_name */
                                };







PyMODINIT_FUNC PyInit_PartiClass(void)
// create the module
{
    PyObject* m;
    
    /* Initialize mpi4py's C-API */
    if (import_mpi4py() < 0) goto bad;

    PyPartitionType.tp_new = PyType_GenericNew;
    PyPartitionType.tp_basicsize=sizeof(PyPartition);
    PyPartitionType.tp_dealloc=(destructor) PyPartition_dealloc;
    PyPartitionType.tp_flags=Py_TPFLAGS_DEFAULT;
    PyPartitionType.tp_doc="Partition objects";
    PyPartitionType.tp_methods=PyPartition_methods;
    //~ PyPartiType.tp_members=Noddy_members;
    PyPartitionType.tp_init=(initproc)PyPartition_init;

    if (PyType_Ready(&PyPartitionType) < 0)
        return NULL;

    m = PyModule_Create(&PartitionModule);
    if (m == NULL)
        return NULL;

    Py_INCREF(&PyPartitionType);
    PyModule_AddObject(m, "Partition", (PyObject *)&PyPartitionType); // Add Parti object to the module
    return m;
    
    bad:
      return NULL;
}
