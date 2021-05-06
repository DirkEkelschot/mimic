#define MPICH_SKIP_MPICXX 1
#define OMPI_SKIP_MPICXX  1
#include <Python.h>
#include <mpi4py/mpi4py.h>
#include <numpy/arrayobject.h>
#include "energy.h"
#include "../src/adapt_io.h"
#include "../src/adapt_partition.h"
#include "../src/adapt_recongrad.h"
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <unordered_map>
/* -------------------------------------------------------------------------- */

double **ptrvector(long n)  {
    double **v;
    v=(double **)malloc((size_t) (n*sizeof(double)));
    if (!v)   {
        printf("In **ptrvector. Allocation of memory for double array failed.");
        exit(0);  }
    return v;
}

int **ptrintvector(long n)  {
    int **v;
    v=(int **)malloc((size_t) (n*sizeof(int)));
    if (!v)   {
        printf("In **ptrintvector. Allocation of memory for int array failed.");
        exit(0);  }
    return v;
}

void free_Cdoublearrayptrs(double **v)  {
    free((char*) v);
}

void free_Cintarrayptrs(int **v)  {
    free((char*) v);
}

double **pydouble2Darray_to_Carrayptrs(PyArrayObject *arrayin)  {
    double **c, *a;
    int i,n,m;
    
    n=arrayin->dimensions[0];
    m=arrayin->dimensions[1];
    c=ptrvector(n);
    a=(double *) arrayin->data;  /* pointer to arrayin data as double */
    for ( i=0; i<n; i++)  {
        c[i]=a+i*m;
    }
    return c;
}

double **pyint2Darray_to_Carrayptrs(PyArrayObject *arrayin)  {
    double **c, *a;
    int i,n,m;
    
    n=arrayin->dimensions[0];
    m=arrayin->dimensions[1];
    c=ptrvector(n);
    a=(double *) arrayin->data;  /* pointer to arrayin data as int */
    for ( i=0; i<n; i++)  {
        c[i]=a+i*m;
    }
    return c;
}



template<typename T>
typename std::vector<T> GetVectorFromPyList(PyObject* pyList, int n)
{
    std::vector<T> output(n);
    PyObject* item   = PyList_GET_ITEM(pyList, 0);

    if (PyLong_Check(item))
    {
        for (int i = 0; i < n; i++)
        {
            item = PyList_GetItem(pyList, i);

            output[i] = PyLong_AsLong(item);

        }
    }
//
    if (PyFloat_Check(item))
    {
        for (int i = 0; i < n; i++)
        {
            item = PyList_GetItem(pyList, i);

            output[i] = PyFloat_AsDouble(item);
        }
    }
    
    return output;
}



template<typename T>
Array<T>* GetArrayFromPyList(PyObject* pyList)
{
    
    int n = (int)PyList_Size(pyList);

    Array<T>* output= new Array<T>(n,1);
    PyObject* item   = PyList_GET_ITEM(pyList, 0);

    if (PyLong_Check(item))
    {
        for (int i = 0; i < n; i++)
        {
            item = PyList_GetItem(pyList, i);

            output->setVal(i,0,PyLong_AsLong(item));

        }
    }
//
    if (PyFloat_Check(item))
    {
        for (int i = 0; i < n; i++)
        {
            item = PyList_GetItem(pyList, i);

            output->setVal(i,0,PyFloat_AsDouble(item));
        }
    }
    
    return output;
}








template<typename T>
typename std::vector<std::vector<T> > GetVectorOfVectorsFromPyList(PyObject* pyList, int n)
{
    std::vector<std::vector<T> > output(n);
    PyObject* itemTuple   = PyList_GET_ITEM(pyList, 0);
    PyObject* item        = PyTuple_GET_ITEM(itemTuple, 0);
    int ncol;
    if (PyLong_Check(item))
    {
        for (int i = 0; i < n; i++)
        {
            item = PyList_GetItem(pyList, i);
            ncol = PyTuple_Size(item);
            std::vector<T> row(ncol);
            for(int j = 0; j < ncol ;  j++)
            {
                row[j] = (int) PyLong_AsLong(PyTuple_GetItem(item, j));
            }


            output[i] = row;

        }
    }
//
    if (PyFloat_Check(item))
    {
        for (int i = 0; i < n; i++)
        {
            item = PyList_GetItem(pyList, i);
            ncol = PyTuple_Size(item);
            std::vector<T> row(ncol);
            for(int j = 0; j < ncol ;  j++)
            {
                row[j] = (double) PyFloat_AsDouble(PyTuple_GetItem(item, j));
            }


            output[i] = row;

        }
    }

    return output;
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







std::map<int,std::vector<int> > getJaggedMapFromPyDict(PyObject* pyDict)
{
    PyObject* pKeys       = PyDict_Keys(pyDict);
    PyObject* pValues     = PyDict_Values(pyDict);
    int len               = (int)PyList_Size(pValues);
    PyObject* firstItem   = PyList_GET_ITEM(pValues, 0);

    //PyDict_Next(index, &pos, &key, &value);
    
    std::map<int,std::vector<int> > mapCpp;
    for (Py_ssize_t i = 0; i < PyDict_Size(pyDict); ++i)
    {
        int key = (int) PyLong_AsLong(PyList_GetItem(pKeys,   i));
        PyObject* value = PyList_GetItem(pValues, i);
        // PyString_AsString returns a char*
        std::vector<int> valvec(PyTuple_Size(value));
        for(int j = 0;j<PyTuple_Size(value);j++)
        {
            int val_entry = (int) PyLong_AsLong(PyTuple_GetItem(value,j));
            valvec[j] = val_entry;
        }
        mapCpp.insert( {key, valvec} );
    }
    
    return mapCpp;
}

//std::map<int,int > getMapFromPyDict(PyObject* pyDict)
//{
//    PyObject* pKeys       = PyDict_Keys(pyDict);
//    PyObject* pValues     = PyDict_Values(pyDict);
//    int len               = (int)PyList_Size(pValues);
//    PyObject* firstItem   = PyList_GET_ITEM(pValues, 0);
//
//    //PyDict_Next(index, &pos, &key, &value);
//
//    std::map<int,int > mapCpp;
//    for (Py_ssize_t i = 0; i < PyDict_Size(pyDict); ++i)
//    {
//        int key = (int) PyLong_AsLong(PyList_GetItem(pKeys,   i));
//        PyObject* value = PyList_GetItem(pValues, i);
//        // PyString_AsString returns a char*
//        std::vector<int> valvec(PyTuple_Size(value));
//        for(int j = 0;j<PyTuple_Size(value);j++)
//        {
//            int val_entry = (int) PyLong_AsLong(PyTuple_GetItem(value,j));
//            valvec[j] = val_entry;
//        }
//        mapCpp.insert( {key, valvec} );
//    }
//
//    return mapCpp;
//}



/*
This function reads the arguments that is called in python
and stored initially as PyObject. The PyObjects are then used
to derive a corresponding C-Object from it. In the case of the
Communicator, we use the PyMPIComm_Get interface function that
generates C-datatype MPI_Comm from a PyObject.
*/

static PyObject *
madam_CalculateEnergy(PyObject *self, PyObject *args)
{
  int argc = PyTuple_GET_SIZE(args);
  
  PyObject *py_comm = NULL;
  PyArrayObject *B;
  PyArrayObject *I;
  int nB, nI;
  double En = 0.0;
  // Reads the python MPIComm, numpy array, integer and string.
  // The O refers to a PyObject or PyArrayObject as an argument.
  // i refers to an integer.
  // s refers to a string.
  
  if (!PyArg_ParseTuple(args, "iOiO:madam", &nB, &B, &nI, &I))
    return NULL;

  // Grab the data from the PyArrayObject.
  
  double *B_c                  = (double *)    B->data;
  double *I_c                  = (double *)    I->data;
  int i;
  
  En = ComputeEnergy(nB, B_c, nI, I_c);
  
  return PyFloat_FromDouble(En);
}









static PyObject *
madam_ReadDomainData(PyObject *self, PyObject *args)
{
    int argc = PyTuple_GET_SIZE(args);
  
    PyObject *py_comm     = NULL;
    PyObject *py_comm_i   = NULL;

    MPI_Comm *comm_p       = NULL;
    MPI_Info *comm_p_info  = NULL;
    
    // Reads the python MPIComm, numpy array, integer and string.
    // The O refers to a PyObject or PyArrayObject as an argument.
    // i refers to an integer.
    // s refers to a string.
  
    const char* gname;
    const char* cname;
    const char* dname;
    
    int ReadFromStats;
    if (!PyArg_ParseTuple(args, "sssiOO:madam", &cname, &gname, &dname, &ReadFromStats,&py_comm,&py_comm_i))
    return NULL;

    comm_p       = PyMPIComm_Get(py_comm);
    comm_p_info  = PyMPIInfo_Get(py_comm_i);
    
    US3D* us3d   = ReadUS3DData(cname,gname,dname,ReadFromStats,*comm_p,*comm_p_info);
  
    PyObject *mesh_io_py = PyList_New((9));
    
    PyObject *xcn_py  = PyList_New((us3d->xcn->getNrow()*us3d->xcn->getNcol()));
    if (!xcn_py)
        return NULL;
    for(int i=0;i<us3d->xcn->getNrow();i++)
    {
        for(int j=0;j<us3d->xcn->getNcol();j++)
        {
            PyObject *xcn = PyLong_FromSsize_t(us3d->xcn->getVal(i,j));
            PyList_SET_ITEM(xcn_py, i*us3d->xcn->getNcol()+j, xcn);
            if (!xcn)
            {
                Py_DECREF(xcn_py);
                return NULL;
            }
        }
    }
    
    PyObject *ien_py  = PyList_New((us3d->ien->getNrow()*us3d->ien->getNcol()));
    PyObject *ief_py  = PyList_New((us3d->ief->getNrow()*us3d->ief->getNcol()));
    PyObject *iee_py  = PyList_New((us3d->ief->getNrow()*us3d->iee->getNcol()));
    PyObject *iet_py  = PyList_New((us3d->iet->getNrow()*1));

    if (!ien_py)
        return NULL;
    for(int i=0;i<us3d->ien->getNrow();i++)
    {
        PyObject *iet = PyLong_FromSsize_t(us3d->iet->getVal(i,0));
        PyList_SET_ITEM(iet_py, i, iet);
      
        for(int j=0;j<us3d->ien->getNcol();j++)
        {
            PyObject *ien = PyLong_FromSsize_t(us3d->ien->getVal(i,j));
            PyList_SET_ITEM(ien_py, i*us3d->ien->getNcol()+j, ien);
          
            if (!ien)
            {
                Py_DECREF(ien_py);
                return NULL;
            }
            if(j<6)
            {
                PyObject *ief = PyLong_FromSsize_t(us3d->ief->getVal(i,j));
                PyObject *iee = PyLong_FromSsize_t(us3d->iee->getVal(i,j));

                PyList_SET_ITEM(ief_py, i*us3d->ief->getNcol()+j, ief);
                PyList_SET_ITEM(iee_py, i*us3d->iee->getNcol()+j, ief);
                
                if (!ief)
                {
                    Py_DECREF(ief_py);
                    return NULL;
                }
                if (!iee)
                {
                    Py_DECREF(iee_py);
                    return NULL;
                }
            }
        }
    }
    
    PyObject *ife_py     = PyList_New((us3d->ief->getNrow()*us3d->iee->getNcol()));
    PyObject *ifn_py     = PyList_New((us3d->ief->getNrow()*us3d->iee->getNcol()));
    PyObject *if_ref_py  = PyList_New((us3d->if_ref->getNrow()*us3d->iee->getNcol()));
    PyObject *if_Nv_py   = PyList_New((us3d->if_Nv->getNrow()*us3d->iee->getNcol()));
    
    if (!ien_py)
        return NULL;
    for(int i=0;i<us3d->ifn->getNrow();i++)
    {
        PyObject *if_ref = PyLong_FromSsize_t(us3d->if_ref->getVal(i,0));
        PyList_SET_ITEM(if_ref_py, i, if_ref);
            
        if (!if_ref)
        {
            Py_DECREF(if_ref_py);
            return NULL;
        }

        PyObject *if_Nv = PyLong_FromSsize_t(us3d->if_Nv->getVal(i,0));
        PyList_SET_ITEM(if_Nv_py, i, if_Nv);
           
        if (!if_Nv)
        {
            Py_DECREF(if_Nv_py);
            return NULL;
        }

        for(int j=0;j<us3d->ifn->getNcol();j++)
        {
            PyObject *ifn = PyLong_FromSsize_t(us3d->ifn->getVal(i,j));
            PyList_SET_ITEM(ien_py, i*us3d->ifn->getNcol()+j, ifn);
                
            if (!ifn)
            {
                Py_DECREF(ifn_py);
                return NULL;
            }
            if(j<2)
            {
                PyObject *ife = PyLong_FromSsize_t(us3d->ife->getVal(i,j));

                PyList_SET_ITEM(ife_py, i*us3d->ife->getNcol()+j, ife);
                if (!ife)
                {
                    Py_DECREF(ife_py);
                    return NULL;
                }
            }
        }
    }
    
    PyList_SET_ITEM(mesh_io_py, 0, xcn_py);
    PyList_SET_ITEM(mesh_io_py, 1, iet_py);
    PyList_SET_ITEM(mesh_io_py, 2, ien_py);
    PyList_SET_ITEM(mesh_io_py, 3, ief_py);
    PyList_SET_ITEM(mesh_io_py, 4, iee_py);
    
    PyList_SET_ITEM(mesh_io_py, 5, ifn_py);
    PyList_SET_ITEM(mesh_io_py, 6, ife_py);
    PyList_SET_ITEM(mesh_io_py, 7, if_ref_py);
    PyList_SET_ITEM(mesh_io_py, 8, if_Nv_py);

    
    delete us3d->xcn;
    
    delete us3d->ien;
    delete us3d->ief;
    delete us3d->iee;
    
    delete us3d->ifn;
    delete us3d->ife;
    delete us3d->if_ref;
    delete us3d->if_Nv;
    //
    return mesh_io_py;
}




static PyObject *
madam_Partition(PyObject *self, PyObject *args)
{
    int argc = PyTuple_GET_SIZE(args);
  
    PyObject *py_comm     = NULL;
    PyObject *py_comm_i   = NULL;

    MPI_Comm *comm_p       = NULL;
    MPI_Info *comm_p_info  = NULL;
    
    // Reads the python MPIComm, numpy array, integer and string.
    // The O refers to a PyObject or PyArrayObject as an argument.
    // i refers to an integer.
    // s refers to a string.
  
    const char* gname;
    const char* cname;
    const char* dname;
    const char* vname;
    int ReadFromStats;
    if (!PyArg_ParseTuple(args, "sssisOO:madam", &cname, &gname, &dname, &ReadFromStats,&vname,&py_comm,&py_comm_i))
    return NULL;

    comm_p       = PyMPIComm_Get(py_comm);
    comm_p_info  = PyMPIInfo_Get(py_comm_i);
    int world_size;
    MPI_Comm_size(*comm_p, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(*comm_p, &world_rank);
    
    std::cout << "Done converting communicator" << std::endl;

    
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

    delete gB;
    delete us3d->interior;

    Partition* P = new Partition(us3d->ien, us3d->iee, us3d->ief, us3d->ie_Nv , us3d->ie_Nf,
                                 us3d->ifn, us3d->ife, us3d->if_ref, us3d->if_Nv,
                                 parmetis_pstate, ien_pstate, ife_pstate,
                                 us3d->xcn, xcn_pstate, Uivar, *comm_p);
      
    

    
    std::map<int,int >::iterator itmap;
    std::map<int,std::vector<int> >::iterator itm;

    //========================================================================
    std::vector<Vert*> verts_part         = P->getLocalVerts();
    PyObject *py_partVerts = PyList_New(verts_part.size());
    for (int q = 0;q<verts_part.size();q++)
    {
        PyObject *row = PyTuple_New(3);
        PyTuple_SET_ITEM(row, 0, PyFloat_FromDouble(verts_part[q]->x));
        PyTuple_SET_ITEM(row, 1, PyFloat_FromDouble(verts_part[q]->y));
        PyTuple_SET_ITEM(row, 2, PyFloat_FromDouble(verts_part[q]->z));

        PyList_SetItem(py_partVerts,q,row);
    }
    PyList_SET_ITEM(part_py, 0, py_partVerts);
    verts_part.clear();
    //========================================================================
    std::vector<int> Loc_Elem_id          = P->getLocElem();
    PyObject *py_locElem_id = PyList_New(Loc_Elem_id.size());
    PyObject *pDict_UState = PyDict_New();
    std::map<int,Array<double>*> Uvaria_map;
    double UvariaV = 0.0;
    
    for (int q = 0;q<Loc_Elem_id.size();q++)
    {
        PyObject *entry = PyLong_FromSsize_t(Loc_Elem_id[q]);
        PyList_SetItem(py_locElem_id,q,entry);
        int gid = Loc_Elem_id[q];
        UvariaV   = Uivar->getVal(q,0);
        Array<double>* Uarr = new Array<double>(1,1);
        Uarr->setVal(0,0,UvariaV);
        Uvaria_map[gid] = Uarr;
    }
    PyList_SET_ITEM(part_py, 1, py_locElem_id);
    Loc_Elem_id.clear();
    //========================================================================
    std::map<int,Array<double>* >::iterator itmapArr;

    for (itmapArr=Uvaria_map.begin();itmapArr!=Uvaria_map.end();itmapArr++)
    {
        PyObject *key = PyLong_FromSsize_t(itmapArr->first);
        PyObject *val = PyFloat_FromDouble(itmapArr->second->getVal(0,0));
        PyDict_SetItem(pDict_UState,key,val);
    }
    PyList_SET_ITEM(part_py, 2, pDict_UState);
    Uvaria_map.clear();
    //========================================================================
    std::map<int,int> gV2lV               = P->getGlobalVert2LocalVert();
    PyObject *pDict_gV2lV = PyDict_New();
    for (itmap=gV2lV.begin();itmap!=gV2lV.end();itmap++)
    {
        PyObject *key = PyLong_FromSsize_t(itmap->first);
        PyObject *val = PyLong_FromSsize_t(itmap->second);
        PyDict_SetItem(pDict_gV2lV,key,val);
    }
    PyList_SET_ITEM(part_py, 3, pDict_gV2lV);
    gV2lV.clear();
//    //========================================================================
    std::map<int,int> LocElem2Nf          = P->getLocElem2Nf();
    PyObject *pDict_LocElem2Nf = PyDict_New();
    for (itmap=LocElem2Nf.begin();itmap!=LocElem2Nf.end();itmap++)
    {
        PyObject *key = PyLong_FromSsize_t(itmap->first);
        PyObject *val = PyLong_FromSsize_t(itmap->second);
        PyDict_SetItem(pDict_LocElem2Nf,key,val);
    }
    PyList_SET_ITEM(part_py, 4, pDict_LocElem2Nf);
    LocElem2Nf.clear();
//    //========================================================================
    std::map<int,int> LocElem2Nv          = P->getLocElem2Nv();
    PyObject *pDict_LocElem2Nv = PyDict_New();
    for (itmap=LocElem2Nv.begin();itmap!=LocElem2Nv.end();itmap++)
    {
        PyObject *key = PyLong_FromSsize_t(itmap->first);
        PyObject *val = PyLong_FromSsize_t(itmap->second);
        PyDict_SetItem(pDict_LocElem2Nv,key,val);
    }
    PyList_SET_ITEM(part_py, 5, pDict_LocElem2Nv);
    LocElem2Nv.clear();
//    //========================================================================
    std::map<int,std::vector<int> > gE2lV = P->getGlobElem2LocVerts();
    PyObject *pDict_gE2lV = PyDict_New();
    for (itm=gE2lV.begin();itm!=gE2lV.end();itm++)
    {
        PyObject *key = PyLong_FromSsize_t(itm->first);
        Py_ssize_t row_size = itm->second.size();
        PyObject *val = PyTuple_New(row_size);
        
        for(int j=0;j<row_size;j++)
        {
            PyTuple_SET_ITEM(val, j, PyLong_FromSsize_t(itm->second[j]));
        }
        PyDict_SetItem(pDict_gE2lV,key,val);
    }
    PyList_SET_ITEM(part_py, 6, pDict_gE2lV);
    gE2lV.clear();
    //========================================================================
    i_part_map*  iee_map                  = P->getIEEpartmap();
    PyObject *pDict_iee = PyDict_New();
    for (itm=iee_map->i_map.begin();itm!=iee_map->i_map.end();itm++)
    {
        PyObject *key = PyLong_FromSsize_t(itm->first);

        Py_ssize_t row_size = itm->second.size();
        PyObject *val = PyTuple_New(row_size);
        for(int j=0;j<row_size;j++)
        {
            PyTuple_SET_ITEM(val, j, PyLong_FromSsize_t(itm->second[j]));
        }
        PyDict_SetItem(pDict_iee,key,val);
    }
    PyList_SET_ITEM(part_py, 7, pDict_iee);
    delete iee_map;
    //========================================================================
    i_part_map*  ief_map                  = P->getIEFpartmap();
    PyObject *pDict_ief = PyDict_New();
    for (itm=ief_map->i_map.begin();itm!=ief_map->i_map.end();itm++)
    {
        PyObject *key = PyLong_FromSsize_t(itm->first);

        Py_ssize_t row_size = itm->second.size();
        PyObject *val = PyTuple_New(row_size);
        for(int j=0;j<row_size;j++)
        {
            PyTuple_SET_ITEM(val, j, PyLong_FromSsize_t(itm->second[j]));
        }
        PyDict_SetItem(pDict_ief,key,val);
    }
    PyList_SET_ITEM(part_py, 8, pDict_ief);
    delete ief_map;
    //========================================================================
    i_part_map*  ifn_map                  = P->getIFNpartmap();
    PyObject *pDict_ifn = PyDict_New();
    for (itm=ifn_map->i_map.begin();itm!=ifn_map->i_map.end();itm++)
    {
        PyObject *key = PyLong_FromSsize_t(itm->first);

        Py_ssize_t row_size = itm->second.size();
        PyObject *val = PyTuple_New(row_size);
        for(int j=0;j<row_size;j++)
        {
            PyTuple_SET_ITEM(val, j, PyLong_FromSsize_t(itm->second[j]));
        }
        PyDict_SetItem(pDict_ifn,key,val);
    }
    PyList_SET_ITEM(part_py, 9, pDict_ifn);
    delete ifn_map;
    //========================================================================
    i_part_map*  if_Nv_map                = P->getIF_Nvpartmap();
    PyObject *pDict_if_Nv = PyDict_New();
    for (itm=if_Nv_map->i_map.begin();itm!=if_Nv_map->i_map.end();itm++)
    {
        PyObject *key = PyLong_FromSsize_t(itm->first);

        Py_ssize_t row_size = itm->second.size();
        PyObject *val = PyTuple_New(row_size);
        for(int j=0;j<row_size;j++)
        {
            PyTuple_SET_ITEM(val, j, PyLong_FromSsize_t(itm->second[j]));
        }
        PyDict_SetItem(pDict_if_Nv,key,val);
    }
    PyList_SET_ITEM(part_py, 10, pDict_if_Nv);
    delete if_Nv_map;
    //========================================================================
    PyList_SET_ITEM(part_py, 11, PyLong_FromSsize_t(P->getNglob_Elem()));


//
//
//
    delete P;

    delete us3d->xcn;

    delete us3d->ien;
    delete us3d->ief;
    delete us3d->iee;

    delete us3d->ifn;
    delete us3d->ife;
    delete us3d->if_ref;
    delete us3d->if_Nv;
    //



    return part_py;
    

}





static PyObject *
madam_ReconstructdUdXi(PyObject *self, PyObject *args)
{
    int argc = PyTuple_GET_SIZE(args);
  
    PyObject *py_comm     = NULL;
    PyObject *py_comm_i   = NULL;

    MPI_Comm *comm_p       = NULL;
    MPI_Info *comm_p_info  = NULL;
    

    int nB, nI;
    // Reads the python MPIComm, numpy array, integer and string.
    // The O refers to a PyObject or PyArrayObject as an argument.
    // i refers to an integer.
    // s refers to a string.
  
    const char* gname;
    const char* cname;
    const char* dname;
    const char* vname;
    double En = 0.0;
    int ReadFromStats;
    PyObject* py_part;
    PyObject* py_Ustate;
    PyObject* py_ghost;
    if (!PyArg_ParseTuple(args, "OOOO:madam", &py_part, &py_Ustate, &py_comm, &py_comm_i))
      return NULL;
    
    comm_p       = PyMPIComm_Get(py_comm);
    comm_p_info  = PyMPIInfo_Get(py_comm_i);
    // Unpakcing the py_part object.
    PyObject *iter = PyObject_GetIter(py_part);
    if (!iter) {
      // error not iterator
    }
    PyObject *key;
    PyObject *value;
    Py_ssize_t pos = 0;
    int entry = 0;
    
    std::map<int,Array<double>* > Ustate_map = getMapOfArraysFromPyDict<double>(py_Ustate);
    //P->AddStateVecForAdjacentElements(Ustate_map,1,*comm_p);
    //Unpacking the partitioning object in Python.
    PyObject* pyDict_LocVerts               = PyList_GET_ITEM(py_part, 0);
    int nV                                  = (int)PyList_Size(pyDict_LocVerts);
    std::vector<std::vector<double> > verts = GetVectorOfVectorsFromPyList<double>(pyDict_LocVerts,nV);
    PyObject* pyDict_LocElem                = PyList_GET_ITEM(py_part, 1);
    int nL                                  = (int)PyList_Size(pyDict_LocElem);
    std::vector<int>   locElem_py           = GetVectorFromPyList<int>(pyDict_LocElem,nL);
    std::map<int,Array<double>* >::iterator itmp;
    // ================================================================================
    PyObject* pyDict_gV2lV = PyList_GET_ITEM(py_part, 3);
    std::map<int,int > gV2lV_py = getMapFromPyDict<int>(pyDict_gV2lV);
    PyObject* pyDict_LocElem2Nf = PyList_GET_ITEM(py_part, 4);
    std::map<int,int > LocElem2Nf_py = getMapFromPyDict<int>(pyDict_LocElem2Nf);
    PyObject* pyDict_LocElem2Nv = PyList_GET_ITEM(py_part, 5);
    std::map<int,int > LocElem2Nv_py = getMapFromPyDict<int>(pyDict_LocElem2Nv);
    // ================================================================================
    PyObject* pyDict_gE2lV = PyList_GET_ITEM(py_part, 6);
    std::map<int,std::vector<int> > gE2lV_py = getJaggedMapFromPyDict(pyDict_gE2lV);
    PyObject* pyDict_iee = PyList_GET_ITEM(py_part, 7);
    std::map<int,std::vector<int> > iee_py = getJaggedMapFromPyDict(pyDict_iee);
    PyObject* pyDict_ief = PyList_GET_ITEM(py_part, 8);
    std::map<int,std::vector<int> > ief_py =  getJaggedMapFromPyDict(pyDict_ief);
    PyObject* pyDict_ifn = PyList_GET_ITEM(py_part, 9);
    std::map<int,std::vector<int> > ifn_py =  getJaggedMapFromPyDict(pyDict_ifn);
    PyObject* pyDict_if_Nv = PyList_GET_ITEM(py_part, 10);
    std::map<int,std::vector<int> > if_Nv_py =  getJaggedMapFromPyDict(pyDict_if_Nv);
    // ================================================================================
    PyObject* nEl_glob = PyList_GET_ITEM(py_part, 11);
    int nGlob = (int) PyLong_AsLong(nEl_glob);

    PyObject* ghost_py = PyList_GET_ITEM(py_part, 12);

    Array<double>* gB = GetArrayFromPyList<double>(ghost_py);
    
    std::map<int,Array<double>* > Ugrad = Py_ComputedUdx_LSQ_US3D(verts,
                                            gE2lV_py,gV2lV_py, locElem_py,
                                            ifn_py,ief_py, iee_py,if_Nv_py,
                                            LocElem2Nf_py,LocElem2Nv_py,nGlob,
                                            Ustate_map, gB, *comm_p);
    
    PyObject *pDict_Ugrad = PyDict_New();
    std::map<int,Array<double>* >::iterator itm;
    for (itm=Ugrad.begin();itm!=Ugrad.end();itm++)
    {
        PyObject *key = PyLong_FromSsize_t(itm->first);
        Py_ssize_t row_size = itm->second->getNrow();
        
        
        PyObject *val = PyTuple_New(row_size);
        for(int j=0;j<row_size;j++)
        {
            PyTuple_SET_ITEM(val, j, PyFloat_FromDouble(itm->second->getVal(j,0)));
        }
        PyDict_SetItem(pDict_Ugrad,key,val);
    }

    return pDict_Ugrad;
}


/*
This function reads the arguments that is called in python
and stored initially as PyObject. The PyObjects are then used
to derive a corresponding C-Object from it. In the case of the
Communicator, we use the PyMPIComm_Get interface function that
generates C-datatype MPI_Comm from a PyObject.
*/

static PyObject *
madam_ReconstructGradient_LSQ(PyObject *self, PyObject *args)
{
  int argc = PyTuple_GET_SIZE(args);
  
  PyObject *py_comm     = NULL;
  PyObject *py_comm_i   = NULL;

  MPI_Comm *comm_p       = NULL;
  MPI_Info *comm_p_info  = NULL;
    

  int nB, nI;
  double En = 0.0;
  // Reads the python MPIComm, numpy array, integer and string.
  // The O refers to a PyObject or PyArrayObject as an argument.
  // i refers to an integer.
  // s refers to a string.
  
  const char* gname;
  const char* cname;
  const char* dname;
  const char* vname;
    
  int ReadFromStats;
  if (!PyArg_ParseTuple(args, "sssisOO:madam", &cname, &gname, &dname, &ReadFromStats,&vname,&py_comm,&py_comm_i))
    return NULL;

  comm_p       = PyMPIComm_Get(py_comm);
  comm_p_info  = PyMPIInfo_Get(py_comm_i);
    
  US3D* us3d   = ReadUS3DData(cname,gname,dname,ReadFromStats,*comm_p,*comm_p_info);
    
  int Nve      = us3d->xcn->getNglob();
  int Nel_part = us3d->ien->getNrow();
  
  
  // ParallelState is an object that allows the user to get Offsets and Nlocs for that input array.
    
  ParallelState* ien_pstate               = new ParallelState(us3d->ien->getNglob(),*comm_p);
  ParallelState* ife_pstate               = new ParallelState(us3d->ifn->getNglob(),*comm_p);
  ParallelState_Parmetis* parmetis_pstate = new ParallelState_Parmetis(us3d->ien,us3d->elTypes,us3d->ie_Nv,*comm_p);
  ParallelState* xcn_pstate               = new ParallelState(us3d->xcn->getNglob(),*comm_p);
    
    
    
    
  Array<double>* Uivar = new Array<double>(Nel_part,1);
  double rhoState,uState,vState,wState,TState,VtotState,aState,MState;
  Array<double>* gB = new Array<double>(us3d->ghost->getNrow(),1);
  
    
    
    
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
    
  delete us3d->interior;

  clock_t t,t1;
  double tmax = 0.0;
  double tn = 0.0;
  t = clock();
    
  // ien -> element2node    map coming from parallel reading.
  // iee -> element2element map coming from parallel reading.
  // ief -> element2face    map coming from parallel reading.
  // ifn -> face2node       map coming from parallel reading.
  // ife -> face2element    map coming from parallel reading.
  //std::cout << "Starting to partition..."<<std::endl;
    
    
  Partition* P = new Partition(us3d->ien, us3d->iee, us3d->ief, us3d->ie_Nv , us3d->ie_Nf,us3d->ifn, us3d->ife, us3d->if_ref, us3d->if_Nv, parmetis_pstate, ien_pstate, ife_pstate, us3d->xcn, xcn_pstate, Uivar, *comm_p);
    

  std::vector<int> LocElem = P->getLocElem();
  std::vector<double> Uvaria  = P->getLocElemVaria();
  std::map<int,Array<double>*> Uvaria_map;
  double UvariaV = 0.0;

  for(int i=0;i<LocElem.size();i++)
  {
      int gid = LocElem[i];
      UvariaV   = Uvaria[i];
      Array<double>* Uarr = new Array<double>(1,1);
      Uarr->setVal(0,0,UvariaV);
      Uvaria_map[gid] = Uarr;
  }
    
  P->AddStateVecForAdjacentElements(Uvaria_map,1,*comm_p);

//  delete us3d->ghost;
  std::map<int,Array<double>* > dUdXi = ComputedUdx_LSQ_US3D(P,Uvaria_map,gB,*comm_p);

  PyObject *pDict = PyDict_New();
    
  std::map<int,Array<double>* >::iterator itm;
  for (itm=dUdXi.begin();itm!=dUdXi.end();itm++)
  {
      PyObject *key = PyLong_FromSsize_t(itm->first);
      PyObject *val = PyFloat_FromDouble(itm->second->getVal(0,0));

      PyDict_SetItem(pDict,key,val);
  }
    
    return pDict;
    
}




static struct PyMethodDef methods[] = {
  {"CalculateEnergy", (PyCFunction)madam_CalculateEnergy, METH_VARARGS, NULL},
  {"ReconstructGradient_LSQ", (PyCFunction)madam_ReconstructGradient_LSQ, METH_VARARGS, NULL},
  {"ReadDomainData", (PyCFunction)madam_ReadDomainData, METH_VARARGS, NULL},
  {"Partition", (PyCFunction)madam_Partition, METH_VARARGS, NULL},
  {"GradU", (PyCFunction)madam_ReconstructdUdXi, METH_VARARGS, NULL},

    
  {NULL, NULL, 0, NULL} /* sentinel */
};


/* --- Python 3 --- */
// This is the default structure of a python module.
// For now we only add methods which is the PartitionData function.
static struct PyModuleDef module = {
  PyModuleDef_HEAD_INIT,
  "madam",       // m_name
  NULL,         // m_doc
  -1,           // m_size
  methods,      // m_methods
  NULL,         // m_reload
  NULL,         // m_traverse
  NULL,         // m_clear
  NULL          // m_free
};

// This constructs the PyModuleDef that is defined above.
PyMODINIT_FUNC
PyInit_madam(void)
{
  PyObject *m = NULL;

  /* Initialize mpi4py's C-API */
  if (import_mpi4py() < 0) goto bad;

  /* Module initialization  */
  m = PyModule_Create(&module);
  if (m == NULL) goto bad;

  return m;

 bad:
  return NULL;
}
