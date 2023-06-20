
#include <stdio.h>
#include <string>
#include "CellSyn.h" //Pair

using namespace std;

//prints all cells voltages
void print_occ(Pair *cell_sizes);
void print_used_mem();
//prints one column/row of cells voltages in a weird format
void print_freq(double **cx_base_v_SOMA, double **cxa_base_v_SOMA,Pair *cell_sizes, double const t);
void close_files(FILE **field_file, int num_field_layers);
void open_files(string output_location,FILE **field_file, int num_field_layers);

//read input from either input file or commandline
//TODO in long term we should have just single parameter class instead of passing zillion args.
void load_input_params(
  int argc,char *argv[],

  double& tmax,
  double&  t3D,
  double&  ttime,
  int&  num_mp_threads,
  int&  print_c_sten,
  int&  fre_print_cs,
  int&  LFP_local_field_effect,
  double&  LFP_lfp_scale,
  int&  LFP_num_field_layers,
  int&  homeo_boost,
  double&  homeo_amp_boost,
  double&  homeo_con_boost,
  double&  homeo_fre_boost,
  double&  homeo_target_f,
  int&  homeo_fre_window,
  int&  homeo_num_regions,
  double& stim_cx_start,
  double& stim_cx_end,
  double& stim_cx_strength,
  int& stim_cx_start_neuron,
  int& stim_cx_end_neuron,
  double& stim_in_start,
  double& stim_in_end,
  double& stim_in_strength,
  int& stim_in_start_neuron,
  int& stim_in_end_neuron,
  double& stim_tc_start,
  double& stim_tc_end,
  double& stim_tc_strength,
  int& stim_tc_start_neuron,
  int& stim_tc_end_neuron,
  double& stim_re_start,
  double& stim_re_end,
  double& stim_re_strength,
  int& stim_re_start_neuron,
  int& stim_re_end_neuron);



