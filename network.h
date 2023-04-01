//Network connections and operations (e.g. LFP, homeostasis, distance)

#include "CellSyn.h" //Pair

//finds the maximum number of connections between each pair of cell types for normal and gap connections not including multiple connections between the same two cells
void find_max_connections(int Num_Types, Pair *cell_sizes,Cell_Info ****cells_info);

//uses the connections input files to figure the number of cells of each type we have
void set_cell_sizes(FILE* connections_file, int &Num_Types, Pair* &cell_sizes, int* &cell_numbers);

//warning this is not very error tolerant of the input file
//also sorry for it not being split into smaller funcs
//reads a file that specifies connections between cells
void create_connections_from_file(Cell_Info *****cells_info, Pair *cell_sizes, FILE *connections_file, int Num_Types);



class Homeostasis{
 public:
  //if and how homeostasis works
  int boost; // 1 for on 0 for off in input file
  double amp_boost; // how much to boost mini amp by
  double con_boost; // how much to boost cx*->cx* connection s by
  double fre_boost; //how much to boost mini fre by
  double target_f; // target frequency for firing
  int *boost_syn; //state variables if we should boost or nerf regions
  int num_regions; 
  int fre_window; //time steps for freqency tracking window 50 time steps per millisecond with tau at 0.02
  double *cut_locs; //boundaries of homeostatic regions
  int *total_region;
  double *frequencies;
  //this function implements the homeostatic mechanism
  void boost_activity(Pair *cell_sizes, CellSyn **group_cells, int num_cells);
  //allocate homeostatic stuff
  void allocate(){
    //allocate homeostatic stuff
    boost_syn = new int[num_regions];
    cut_locs = new double[num_regions+1];
    total_region = new int[num_regions];
    frequencies = new double[num_regions];
    double cut_inc = 1.0/num_regions;
    cut_locs[0] = 0.0;
    for(int iter = 0; iter < num_regions; iter++){
      cut_locs[iter+1] = cut_inc * (double)(iter+1); 
      boost_syn[iter] = 1;
    }//for
  }//allocate

  //calculates,prints,and broadcasts average spiking frequency of network
  void spike_fre_calc(Pair *cell_sizes, int *cell_numbers);
};//Homeostasis


//all stuff for local field effects
class LocalFieldPotential {
 public:
  //local field potential do not turn on unless using openmp
  int local_field_effect; //turns field on and off
  double lfp_scale; //how field is scaled
  int num_field_layers; //how many field layers to calc be careful with this 


  double **cx_base_v_SOMA;
  double **cxa_base_v_SOMA;
  double **cx5_local_field;
  double **cx6_local_field;
  double **cx5_soma_local_field;
  double **cx6_soma_local_field;
  double ***field_effect;
  FILE **field_file;
  //ideally all should be set via constructor
  void init(){
    field_file = new FILE*[num_field_layers];
  }
  void allocate_state_save(Pair *cell_sizes);
  //TODO !D only?
  //local field effect probably only works 1D. And certainly only works if CX and CXa and cx6 are the same size
  void apply_field(int ii,double time,double ttime,double t,Pair *cell_sizes);

}; //LocalFieldPotential

