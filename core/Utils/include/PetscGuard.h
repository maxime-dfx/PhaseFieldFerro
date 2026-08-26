#pragma once

#ifdef USE_PETSC
#include <petscsys.h>
#endif

class PetscGuard {
public:
    PetscGuard([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
#ifdef USE_PETSC
        PetscInitialize(&argc, &argv, nullptr, nullptr);
#endif
    }
    
    ~PetscGuard() {
#ifdef USE_PETSC
        PetscFinalize();
#endif
    }
};