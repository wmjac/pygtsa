#ifndef __RANDOM_H__
#define __RANDOM_H__

#include <gsl/gsl_rng.h>

extern gsl_rng * rng;
#define random_uniform() (gsl_rng_uniform (rng))

void random_init ();

#endif /* __RANDOM_H__ */
