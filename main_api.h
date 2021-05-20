//This needs to be done when linking C programs to Fortran. The reason is name mangling
extern "C" {
int madam_api(const char* fn_grid,const char* fn_conn,const char* fn_data, MPI_Comm comm, int world_size, int world_rank);
}
