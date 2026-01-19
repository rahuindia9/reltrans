#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// TODO: is there a header file for these function definitions?

// either way, forward declare this as an extern, and let it be the linker's
// problem
extern void tdreltransdcp_(float *, int *, float *, int *, float *);

#define LOG_ERR(...)                                                           \
  fprintf(stderr, "ERR  : " __VA_ARGS__);                                      \
  putc('\n', stderr)

typedef struct reltrans_dcp_parameters {
  float h,    // lamp post height
      a,      // spin
      inc,    // inclination
      rin,    // inner radius
      rout,   // outer radius
      zcos,   // cosmological redshift
      gamma,  // photon index
      logxi,  // logξ ionisation parameter
      afe,    // iron abundance
      lognep, // electron abundance (?)
      kte,    // electron temperature in observer frame
      nh,     // hydrogen column density
      boost,  // boosting factor (ad-hoc normalisation)
      mass,   // black hole mass in solar units
      flo_hz, // lowest frequency in band
      fhi_hz, // highest frequency in band
      re_im,  // 1 -> Re, 2 -> Im, 3 -> modulus, 4 -> time lag, 5 -> folded
              // modulus, 6 -> folded time lag
      del_a, del_ab, g, telescope_response;
} RT_DCP_Params;

RT_DCP_Params default_parameters() {
  RT_DCP_Params params = {
      .h = 6.0,
      .a = 0.998,
      .inc = 30.0,
      .rin = 1.0,
      .rout = 1e3,
      .zcos = 0.0,
      .gamma = 2.0,
      .logxi = 3.0,
      .afe = 1.0,
      .lognep = 15,
      .kte = 60.0,
      .nh = 0.0,
      .boost = 1.0,
      .mass = 4.6e7,
      .flo_hz = 0.0,
      .fhi_hz = 0.0,
      .re_im = 1.0,
      .del_a = 0.0,
      .del_ab = 0.0,
      .g = 0.0,
      .telescope_response = 1,
  };
  return params;
}

int main() {
  const char *output_file = "output.txt";
  const float e_min = 0.1;
  const float e_max = 1000.0;
  int e_num = 1001;

  float *energy = malloc(sizeof(float) * e_num);
  if (energy == NULL) {
    LOG_ERR("Failed to allocate energy array");
    return 1;
  }

  float *output = malloc(sizeof(float) * e_num);
  if (output == NULL) {
    LOG_ERR("Failed to allocate output array");
    return 1;
  }

  RT_DCP_Params params = default_parameters();

  // logarithmic energy grid
  for (int i = 0; i < e_num; ++i) {
    energy[i] = e_min * powf(e_max / e_min, ((float)i) / e_num);
  }

  // zero the output buffer
  memset(output, 0, e_num);

  clock_t time = clock();
  int ifl = 1;
  e_num -= 1;
  tdreltransdcp_(energy, &e_num, (float *)&params, &ifl, output);
  time = clock() - time;

  printf("Total Call: %.6f seconds\n", ((double)time) / CLOCKS_PER_SEC);

  printf("Serialising output to file...\n");
  FILE *fp = fopen(output_file, "w");
  if (fp == NULL) {
    LOG_ERR("Failed to open/create output file %s", output_file);
    return 1;
  }

  for (int i = 0; i < e_num - 1; ++i) {
    fprintf(fp, "%f %f\n", energy[i], output[i]);
  }

  fclose(fp);
  free(output);
  free(energy);

  return 0;
}
