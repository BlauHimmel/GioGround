#pragma once

#include <time.h>

// Used by function

// The function just use this parameter. This variable might be modified, depending on the keyword "const"
#define In

// The function will return value using parameter, but the buffer should be allocated in advance
#define Out

// The function first use this parameter, then modify it, and finally return the modified value using it
#define InOut

// This parameter is optional. It would be ignored the it was null.
#define Optional

// Used by constructor of class and struct

// The memory is managed by the class, which meas the ownership has been moved
#define Managed

// The class just hold the reference of it, but do not manage the memory
#define Ref

// The constructor just use this parameter but do not hold it.
#define Used

#define TIMER_START(Var) clock_t __Start##Var, __Finish##Var; double Var; __Start##Var = clock();
#define TIMER_END(Var) __Finish##Var = clock();  Var = double(__Finish##Var - __Start##Var) / (CLOCKS_PER_SEC);