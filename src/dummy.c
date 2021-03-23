/* dummy file with no function but creates symbols to pass R CMD check */

extern void R_registerRoutines();
extern void R_useDynamicSymbols();

void dummy_() {
  R_registerRoutines();
  R_useDynamicSymbols();
}
