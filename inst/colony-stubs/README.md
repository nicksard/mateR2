# Running Simu2.exe without Intel's runtime

`simu2.exe` imports seven Fortran MPI entry points from `impi.dll` and four
ordinals (110, 120, 726, 747) from `libiomp5md.dll`. On a machine with COLONY
installed, those DLLs sit in the Colony program folder and should simply be
copied next to `Input3.Par`, as the user guide says.

These stubs exist for the case where they are not available, such as testing the
output of `write_colony_sim()` on Linux under wine. They implement a single-rank
MPI world and, for the OpenMP ordinals, the four functions the call sites reveal:
`__kmpc_begin` and `__kmpc_end` at program entry and exit, and a malloc/free
pair. Nothing in `simu2.exe` calls `__kmpc_fork_call`, so serial stand-ins are
sufficient and results are unaffected.

    x86_64-w64-mingw32-gcc -O2 -shared -o impi.dll        impi.c -Wl,--kill-at
    x86_64-w64-mingw32-gcc -O2 -shared -o libiomp5md.dll  iomp.c iomp.def

Put both DLLs and `simu2.exe` in the directory holding `Input3.Par` and run it.
This does not replace the real Intel runtime for a production COLONY analysis.
