#include "spartan.h"
#include <stdio.h>
void spartan_nullprint(const char *s){
  return;
}
void spartan_printout(const char *s){
  printf("%s\n",s);
}

#ifndef DEBUG 
spartan_printfunc spartan_printfunction = spartan_nullprint;
#else 
spartan_printfunc spartan_printfunction = spartan_printout; 
#endif 

void spartan_setprintfunction(spartan_printfunc printfunction){
  spartan_printfunction = printfunction;
}

void spartan_print(const char *s){
  spartan_printfunction(s);
}

