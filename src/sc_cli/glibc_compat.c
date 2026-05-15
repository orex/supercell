/*
 * Optional glibc symbol-version compatibility shim.
 *
 * Enabled via the LEGACY_GLIBC_EXP CMake option (OFF by default). When ON, the
 * supercell executable is linked with -Wl,--wrap=exp, redirecting every call
 * site of exp(double) to __wrap_exp below, which in turn calls exp's older
 * GLIBC_2.2.5-versioned symbol. Without this shim, building on glibc >= 2.29
 * pins the resulting binary to glibc 2.29 (the default version of exp moved
 * there in Feb 2019). With the shim, that lower bound drops to 2.17.
 *
 * Only meaningful on Linux/x86_64 (the platform whose glibc bumped exp@@). On
 * any other target the file compiles to an empty TU.
 *
 * No detectable numerical impact for supercell's use of exp -- the only call
 * site is the Ewald-summation Gaussian exp(-K^2 / (4*eta)), where the two
 * versions agree to within a ULP and the algorithm has its own tolerance many
 * orders of magnitude looser.
 */

#if defined(__GNUC__) && defined(__linux__) && defined(__x86_64__)

__asm__(".symver exp_glibc_2_2_5, exp@GLIBC_2.2.5");
extern double exp_glibc_2_2_5(double);

double __wrap_exp(double x);
double __wrap_exp(double x) { return exp_glibc_2_2_5(x); }

#endif
