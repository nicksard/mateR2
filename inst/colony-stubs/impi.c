/* Serial stand-in for Intel MPI's Fortran bindings. simu2.exe imports exactly
   seven symbols; a single-rank world satisfies all of them. Fortran passes
   every argument by reference. */
#include <windows.h>
__declspec(dllexport) void MPI_INIT(int *ierr){ *ierr = 0; }
__declspec(dllexport) void MPI_FINALIZE(int *ierr){ *ierr = 0; }
__declspec(dllexport) void MPI_BARRIER(int *comm, int *ierr){ (void)comm; *ierr = 0; }
__declspec(dllexport) void MPI_COMM_RANK(int *comm, int *rank, int *ierr){ (void)comm; *rank = 0; *ierr = 0; }
__declspec(dllexport) void MPI_COMM_SIZE(int *comm, int *size, int *ierr){ (void)comm; *size = 1; *ierr = 0; }
__declspec(dllexport) void MPI_INITIALIZED(int *flag, int *ierr){ *flag = 1; *ierr = 0; }
__declspec(dllexport) void MPI_ABORT(int *comm, int *code, int *ierr){ (void)comm; (void)ierr; ExitProcess(*code); }
BOOL WINAPI DllMain(HINSTANCE h, DWORD r, LPVOID v){ (void)h;(void)r;(void)v; return TRUE; }
