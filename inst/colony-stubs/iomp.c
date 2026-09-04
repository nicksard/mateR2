/* Serial stand-in for the four libiomp5md.dll ordinals simu2.exe imports.
   Identified from the call sites: 110 and 120 run at entry and exit with a
   source-location pointer (__kmpc_begin / __kmpc_end); 747 and 726 sit on an
   allocate/free pair guarded by a runtime flag. */
#include <windows.h>
#include <stdlib.h>
__declspec(dllexport) void  kmpc_begin(void *loc, int flags){ (void)loc; (void)flags; }
__declspec(dllexport) void  kmpc_end(void *loc){ (void)loc; }
__declspec(dllexport) void *kmpc_malloc(size_t n){ return malloc(n); }
__declspec(dllexport) void  kmpc_free(void *p){ free(p); }
BOOL WINAPI DllMain(HINSTANCE h, DWORD r, LPVOID v){ (void)h;(void)r;(void)v; return TRUE; }
