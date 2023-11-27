#include "io.h"
#include <sstream>
#include "params.h"
#include <string.h>

double print_receive(int m , int n, enum Cell_Type type); //still in main

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
  int& stim_re_end_neuron,
  double& ach_cx_awake,
  double& ach_th_awake,
  double& ha_awake,
  double& dc_tc

 ){
  
  add_double_param(  tmax  );
  add_double_param(  t3D   );
  add_double_param(  ttime );
  add_int_param(  num_mp_threads );
  add_int_param(  print_c_sten   );
  add_int_param(  fre_print_cs   );
  add_int_param(  LFP_local_field_effect );
  add_double_param(  LFP_lfp_scale );
  add_int_param(  LFP_num_field_layers );
  add_int_param(  homeo_boost  );
  add_double_param(  homeo_amp_boost  );
  add_double_param(  homeo_con_boost  );
  add_double_param(  homeo_fre_boost  );
  add_double_param(  homeo_target_f   );
  add_int_param(  homeo_fre_window    );
  add_int_param(  homeo_num_regions   );

  add_double_param(  stim_cx_start  );
  add_double_param(  stim_cx_end  );
  add_double_param(  stim_cx_strength  );
  add_int_param(  stim_cx_start_neuron   );
  add_int_param(  stim_cx_end_neuron   );
 
  add_double_param(  stim_in_start  );
  add_double_param(  stim_in_end  );
  add_double_param(  stim_in_strength  );
  add_int_param(  stim_in_start_neuron   );
  add_int_param(  stim_in_end_neuron   );

  add_double_param(  stim_tc_start  );
  add_double_param(  stim_tc_end  );
  add_double_param(  stim_tc_strength  );
  add_int_param(  stim_tc_start_neuron   );
  add_int_param(  stim_tc_end_neuron   );

  add_double_param(  stim_re_start  );
  add_double_param(  stim_re_end  );
  add_double_param(  stim_re_strength  );
  add_int_param(  stim_re_start_neuron   );
  add_int_param(  stim_re_end_neuron   );
  add_double_param(  ach_cx_awake  );
  add_double_param(  ach_th_awake  );
  add_double_param(  ha_awake  );
  add_double_param(  dc_tc  );

  assert(load_parameters(argv[1]));
  assert(cmdline_parameters(argc,argv));
  print_parameters();

} //load_input_params

//giant list of all the output files some are probably not used
FILE *flocal, *f2, *f3, *f4, *f6, *f7, *f8, *f9, *f10, *f11, *f12, *f13, *f14, *f15, *f16, *f17, *f18, *f19,*f20,*f21,*f22,*f23,*f24,*f25,*f26,*f27,*f28; 

void print_occ(Pair *cell_sizes){

  int i =0;
  int j =0;

  for(i = 0; i < cell_sizes[E_RE].x; ++i){
    for(j = 0; j < cell_sizes[E_RE].y; ++j){
      fprintf(f6,"%lf ", print_receive(i,j,E_RE));
    }
    fprintf(f6,"\n");
  }

  for(i = 0; i < cell_sizes[E_TC].x; ++i){
    for(j = 0; j < cell_sizes[E_TC].y; ++j){
      fprintf(f8,"%lf ", print_receive(i,j,E_TC));
    }
    fprintf(f8,"\n");
  }

  for(i = 0; i < cell_sizes[E_TCa].x; ++i){
    for(j = 0; j < cell_sizes[E_TC].y; ++j){
      fprintf(f16,"%lf ", print_receive(i,j,E_TCa));
    }
    fprintf(f16,"\n");
  }

  for(i = 0; i < cell_sizes[E_CX].x; ++i){
    for(j = 0; j < cell_sizes[E_CX].y; ++j){
      fprintf(f10,"%lf ", print_receive(i,j,E_CX));
    }
    fprintf(f10,"\n");
  }
  
  for(i = 0; i < cell_sizes[E_CXa].x; ++i){
    for(j = 0; j < cell_sizes[E_CXa].y; ++j){
      fprintf(f14,"%lf ", print_receive(i,j,E_CXa));
    }
    fprintf(f14,"\n");
  }

  for(i = 0; i < cell_sizes[E_IN].x; ++i){
    for(j = 0; j < cell_sizes[E_IN].y; ++j){
      fprintf(f12,"%lf ", print_receive(i,j,E_IN));
    }
    fprintf(f12,"\n");
  }

  // TB IN layers
  for(i = 0; i < cell_sizes[E_INa].x; ++i){
    for(j = 0; j < cell_sizes[E_INa].y; ++j){
      fprintf(f18,"%lf ", print_receive(i,j,E_INa));
    }
    fprintf(f18,"\n");
  }
  
  for(i = 0; i < cell_sizes[E_CX6].x; ++i){
    for(j = 0; j < cell_sizes[E_CX6].y; ++j){
      fprintf(f24,"%lf ", print_receive(i,j,E_CX6));
    }
    fprintf(f24,"\n");
  }

  for(i = 0; i < cell_sizes[E_IN6].x; ++i){
    for(j = 0; j < cell_sizes[E_IN6].y; ++j){
      fprintf(f26,"%lf ", print_receive(i,j,E_IN6));
    }
    fprintf(f26,"\n");
  }
}

void print_freq(double **cx_base_v_SOMA, double **cxa_base_v_SOMA, Pair *cell_sizes, double const t){

  int i = 0;
  //int j = 0;

  //for(i=0; i <cell_sizes[3].x; i++){
  //for(j=0; j<cell_sizes[3].y; j++){
  //cx_base_v_SOMA[i][j] = print_receive(i,j,E_CX);
  //}
  //}
  
  //for(i=0; i <cell_sizes[4].x; i++){
  //for(j=0; j<cell_sizes[4].y; j++){
  //cxa_base_v_SOMA[i][j] = print_receive(i,j,E_CXa);
  //}
  //}
  
  // av=0;
  //for(j = 0; j < cell_sizes[3].y; ++j){
  //for(i = 0; i < cell_sizes[3].x; ++i){
  //av=av + cx_base_v_SOMA[i][j];
  //}
  //}
      
  //ava=0;
  //for(j = 0; j < cell_sizes[4].y; ++j){
  //for(i = 0; i < cell_sizes[4].x; ++i){
  //  ava=ava + cxa_base_v_SOMA[i][j];
  //}
  //}

  //synAMPA=0;
  //synNMDA=0;
  //synGABA=0;
  //synAMPAtc=0;

  //for(j = 0; j < cell_sizes[3].y; ++j)
  //for(i = 0; i < cell_sizes[3].x; ++i){
      //order of print_receives is important
      //double a_tc_cx_I =print_receive(i,j,E_CX);
      //double ga_in_cx_I =print_receive(i,j,E_CX);
      //double nmda_cx_cx_I =print_receive(i,j,E_CX);
      //double nmda_cxa_cx_I =print_receive(i,j,E_CX);
      //double a_cxa_cx_I =print_receive(i,j,E_CX);
      //double a_cx_cx_I = print_receive(i,j,E_CX);
      //synAMPA=synAMPA+a_cx_cx_I+a_cxa_cx_I;
      //synNMDA=synNMDA+nmda_cx_cx_I+nmda_cxa_cx_I;
      //synGABA=synGABA+ga_in_cx_I;
      //synAMPAtc=synAMPAtc+a_tc_cx_I;
      //}
  
  fprintf(f7,"%lf ", t);
  for(i = 0; i < cell_sizes[E_RE].x; ++i){
    fprintf(f7,"%lf ",print_receive(i,cell_sizes[E_RE].y/2,E_RE));
  }
  fprintf(f7,"\n");
  
  fprintf(f9,"%lf ", t);
  for(i = 0; i < cell_sizes[E_TC].x; ++i){
    fprintf(f9,"%lf ", print_receive(i,cell_sizes[E_TC].y/2,E_TC));
  }
  fprintf(f9,"\n");

  fprintf(f17,"%lf ", t);
  for(i = 0; i < cell_sizes[E_TCa].x; ++i){
    fprintf(f17,"%lf ", print_receive(i,cell_sizes[E_TCa].y/2,E_TCa));
  }
  fprintf(f17,"\n");
  
  //print time
  fprintf(f11,"%lf ", t);
  
  //fprintf(f11,"%lf %lf %lf %lf %lf %lf ",t, av/(cell_sizes[3].x*cell_sizes[3].y), synAMPA/(cell_sizes[3].x*cell_sizes[3].y), synNMDA/(cell_sizes[3].x*cell_sizes[3].y), synGABA/(cell_sizes[3].x*cell_sizes[3].y), synAMPAtc/(cell_sizes[3].x*cell_sizes[3].y)); 
  //for(i = 0; i < cell_sizes[3].x; ++i){
  //fprintf(f11,"%lf ", cx_base_v_SOMA[i][cell_sizes[3].y/2]);
  //}
  //fprintf(f11,"\n");

  // double *temp_cx  = new double [cell_sizes[E_CX].x];
  double total = 0.0;
  double temp_cx=0.0;

  for(i = 0; i < cell_sizes[E_CX].x; ++i){
    temp_cx = print_receive(i,cell_sizes[E_CX].y/2,E_CX);
    
    //print individual cell data
    fprintf(f11,"%lf ", temp_cx);
    total = total + temp_cx;
  }

  //print average 
  fprintf(f11,"%lf ", total/cell_sizes[E_CX].x);

  fprintf(f11,"\n");
  //fflush(f11);


  //print time
  fprintf(f15,"%lf ", t);
  
  //fprintf(f15,"%lf %lf ",t,ava/(cell_sizes[4].x*cell_sizes[4].y));
  //for(i = 0; i < cell_sizes[4].x; ++i){
  //fprintf(f15,"%lf ", cxa_base_v_SOMA[i][cell_sizes[4].y/2]);
  //}
  //fprintf(f15,"\n");

  
  double temp_cxa  = 0.0;
  double total_cxa = 0.0;

  
  for(i = 0; i < cell_sizes[E_CXa].x; ++i){
    temp_cxa = print_receive(i,cell_sizes[E_CXa].y/2,E_CXa);
    total_cxa = total_cxa + temp_cxa;

    //print individual cell data
    fprintf(f15,"%lf ", temp_cxa);
  }
  //print average 
  // fprintf(f15,"%lf 0.0 0.0 0.0 0.0 ", total_cxa/cell_sizes[E_CXa].x);
  fprintf(f15,"%lf ", total_cxa/cell_sizes[E_CXa].x);
  fprintf(f15,"\n");

  fprintf(f13,"%lf ", t);
  for(i = 0; i < cell_sizes[E_IN].x; ++i){
    fprintf(f13,"%lf ", print_receive(i,cell_sizes[E_IN].y/2,E_IN));
  }
  fprintf(f13,"\n");

  fprintf(f19,"%lf ", t); // TB IN layers
  for(i = 0; i < cell_sizes[E_INa].x; ++i){
    fprintf(f19,"%lf ", print_receive(i,cell_sizes[E_INa].y/2,E_INa));  // TB IN layers
  }
  fprintf(f19,"\n"); // TB IN layers

  fprintf(f25,"%lf ", t); // TB IN layers
  for(i = 0; i < cell_sizes[E_CX6].x; ++i){
    fprintf(f25,"%lf ", print_receive(i,cell_sizes[E_CX6].y/2,E_CX6));  // TB IN layers
  }
  fprintf(f25,"\n"); // TB IN layers

  
  fprintf(f27,"%lf ", t); // TB IN layers
  for(i = 0; i < cell_sizes[E_IN6].x; ++i){
    fprintf(f27,"%lf ", print_receive(i,cell_sizes[E_IN6].y/2,E_IN6));  // TB IN layers
  }
  fprintf(f27,"\n"); // TB IN layers
}

void open_files(string output_location,FILE **field_file, int num_field_layers){
  printf("open files\n");
  int i = 0;
  for(i=0; i<num_field_layers; i++){
    stringstream ss;
    ss << i;
    field_file[i] = fopen((output_location+"field_file_" + ss.str()  ).c_str(), "w");
  }
  
  f2 = fopen((output_location+"dat").c_str(), "w");
  f6 = fopen((output_location+"graf_re").c_str(), "w");
  if (!(f7=fopen((output_location+"time_re").c_str(), "w"))){
      printf("probably out put folder doesn't exist\n");
      exit(1); 
  }
  //f7 = fopen((output_location+"time_re").c_str(), "w");
  f8 = fopen((output_location+"graf_tc").c_str(), "w");
  f9 = fopen((output_location+"time_tc").c_str(), "w");
  f10 = fopen((output_location+"graf_cx").c_str(), "w");
  f11 = fopen((output_location+"time_cx").c_str(), "w");
  f12 = fopen((output_location+"graf_in").c_str(), "w");
  f13 = fopen((output_location+"time_in").c_str(), "w");
  //TANYA MODIF layers 
  f14 = fopen((output_location+"graf_cxa").c_str(), "w");
  f15 = fopen((output_location+"time_cxa").c_str(), "w");
  f16 = fopen((output_location+"graf_tca").c_str(), "w");
  f17 = fopen((output_location+"time_tca").c_str(), "w");
  f18 = fopen((output_location+"graf_ina").c_str(), "w"); // TB IN layers
  f19 = fopen((output_location+"time_ina").c_str(), "w"); // TB IN layers
  //END TANYA MODIF layers 
  f20 = fopen((output_location+"time_G_AMPA0_CX_CX").c_str(), "w");
  f21 = fopen((output_location+"time_G_AMPA0_CXa_CX").c_str(), "w");
  f22 = fopen((output_location+"time_G_AMPA0_CX_CXa").c_str(), "w"); // TB IN layers
  f23 = fopen((output_location+"time_G_AMPA0_CXa_CXa").c_str(), "w"); // TB IN layers

  f24 = fopen((output_location+"graf_cx6").c_str(), "w"); // TB IN layers
  f25 = fopen((output_location+"time_cx6").c_str(), "w"); // TB IN layers
  f26 = fopen((output_location+"graf_in6").c_str(), "w"); // TB IN layers
  f27 = fopen((output_location+"time_in6").c_str(), "w"); // TB IN layers

  f28 = fopen((output_location+"cx_cx_g_ampa0").c_str(), "w");
  printf("files open for write\n");
}


//TODO convert this to something sensible
void close_files(FILE **field_file, int num_field_layers){

  int i = 0;
  for(i=0; i<num_field_layers; i++){
    fclose(field_file[i]);
  }

  if(f2!=NULL){
    fclose(f2);
  }
  if(f6!=NULL){
    fclose(f6);
  }
  if(f7!=NULL){
    fclose(f7);
  }
  if(f8!=NULL){
    fclose(f8);
  }
  if(f9!=NULL){
    fclose(f9);
  }
  if(f10!=NULL){
    fclose(f10);
  }
  if(f11!=NULL){
    fclose(f11);
  }
  if(f12!=NULL){
    fclose(f12);
  }
  if(f13!=NULL){
    fclose(f13);
  }
  if(f14!=NULL){
    fclose(f14);
  }
  if(f15!=NULL){
    fclose(f15);
  }
  if(f16!=NULL){
    fclose(f16);
  }
  if(f17!=NULL){
    fclose(f17);
  }
  if(f18!=NULL){
    fclose(f18);
  }  // TB IN layers
  if(f19!=NULL){
    fclose(f19);
  }
  if(f20!=NULL){
    fclose(f20);
  }
  if(f21!=NULL){
    fclose(f21);
  }
  if(f22!=NULL){
    fclose(f22);
  }
  if(f23!=NULL){
    fclose(f23);
  }
  if(f24!=NULL){
    fclose(f24);
  }
  if(f25!=NULL){
    fclose(f25);
  }
  if(f26!=NULL){
    fclose(f26);
  }
  if(f27!=NULL){
    fclose(f27);
  }
  if(f28!=NULL){
    fclose(f28);
  }
}

//returns MB
void print_used_mem(){
  FILE* proc = fopen("/proc/self/status", "r");
  string res;
  char line[128];

  while (fgets(line, 128, proc)){
      if (!strncmp(line, "VmRSS:", 6)){
	  printf("Resident RAM: %s",line);
          break;
      }
  }
  fclose(proc);
} 
