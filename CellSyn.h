#ifndef CellSyn_h
#define CellSyn_h

#include "currents.h"    //header for all classes describing currents and basic cell types
#include <list>

using namespace std;

//it is important this comes before the other includes
#define I_HH    1     // 1 -- HH model; 0 -- for Maps;
// could be used with higher numbers for other models

// Removed this fix instead, modify the connection file for this fix 
//#define Nlayers 3  //quick fix to scale cortical exit conductance by number of layers. 

//------------Number of ODE for each cell -------------------------------

//TODO these numbers directly transfer to CellSyn.num_b_eq; maybe move it to constructor?
#define N_RE 7
#define N_TC 12 
#define N_GB 2
#define N_TCa 12 

#if (I_HH==1)
#define N_DEND   9
#define N_SOMA   4 //3
#else 
#define N_DEND   0
#define N_SOMA   0
#endif

#define N_CX     (N_DEND + N_SOMA)
#define N_CXa     (N_DEND + N_SOMA)
#define N_CX6     (N_DEND + N_SOMA)
#define N_IN     N_CX   //4
#define N_INa     N_CX 	// TB IN layer
#define N_IN6     N_CX 	

//other defines
#define TAU      0.02 //integration time step
#define TAU_Map  0.5  // time step for maps
#define TAUr     25   //(TAU_Map/TAU)


//information about each cell as read from network file. CellSyn contain members relater rather to computation
typedef struct{
  int *num_connects_norm; //number of incoming dendrites of each type
  int *num_connects_gap; //number of gap connections
  int *num_connects_total;

  //same as 3 above but does not include repeated connections involving the same two cells
  int *num_cell_cons; //includes every thing including short range and gap 
  int *num_cell_cons_gap; //does not include short range
  int *num_cell_cons_norm; //does not include short range

  int total_connects; //total numper of incoming 
  int cell_index; //the index identifing where in "cells" array this cell is 
  Syn_Info **syn_info;
}Cell_Info;

typedef struct{
  int x;
  int y;
}Pair;


/* bool ***current_spike=NULL; // Array which says if neuron spiked at previous time. */


//abtract class for the combination of a cell and its synapsis
//7 classes derive from it. TODO Fix remaining classes without I_HH switch
//The class includes a base_cell that can be one of 3 types CX,RE, or TC
class CellSyn{

 protected:
  int lhs_size; // how big y,yk,f and s are
  int num_b_eq; //used to decide how far into y an f you have to go to find the indexes for I_GB
  double  *yk; 
  double *s; 
  double cur_time;
  double xk;
  double tau; //time step size
  int current_step; //which of the 4 integration steps we are on

  Cell_Info *cell_info; //information about the cell sent from the neuro file

  double AMPA_Rsum; //Partial sums for R variables in synapse
  //double GABAA_Rsum;
  //list of giant synapses for this cell.
  //individual synapses are responsible for having proper pointer to this list
  list<GiantSyn*> giant_syns;
  //find or create new giant synapse in case it was not created yet
  GiantSyn* get_GiantSyn(enum Cell_Type cell, enum Syn_Type syn);
  //clear spikes from last time step
  void reset_GiantSyns();
  //call calcs for all giant syns and return sum of currents
  double calc_GiantSyns();

 public:

  Syn **syns; //synapses
  int num_syns; // number of synapses
  drand48_data *rand_buffer;
  double *y, *f;
  int m; // cell's location
  int n; // cell's location
  int flip; // if cell is spiking
  int flip_for_maps; // check if the cell is spiking for maps
  double old_v; // voltage before last calc of cell
  int type; // type of this cell
  int num_spikes; // how many spikes this window
  int ismap; // map or not


  BaseCell *base_cell; //Either an RE,CX, or TC cell defined in currents

  CellSyn(int i,int j,int b_size,int t_size,Cell_Info *cell_info,int type);

  //this function creates all the cells(CellSyns) and does some basic initialization
  //all arguments are inputs.
  static CellSyn* initialize(int type,int i, int j, Cell_Info ****cells_info, string output_location);

  //TODO no one using at the moment except of commented code
  double zero_one_rand(){
    double my_rand;
    //drand48_r(rand_buffer,&my_rand);
    my_rand = rand()/(RAND_MAX + 1.0);
    return my_rand;
  }
  int get_flip(){
    return this->flip;
  }
  int get_flip_for_maps(){
    return this->flip_for_maps;
  }

  //manages spiking event (signal & store needed info)
  void signal_spike();
  void reset_flip_maps();

  //RK computation, TODO is this fragmentation really needed?
  void zero_temps();
  void reset_y();
  void step();

  //step - 1 for first step of runge-kutta 0 for all other steps
  void calc(double x, double *y_ini, double *f_ini, int step);
  //create synapses and init base cell
  virtual void memory(Syn_Info **syn_info,int num_syns) = 0;
  //each type of neuron needs to define its own matrix of used Giant synapses
  //we ditched virtuality and use generic version for all from->to types right now
  //otherwise those need to be defined in particular neurons types and used in ::memory()
  GiantSyn* define_GiantSyn(enum Cell_Type cell, enum Syn_Type syn);
  //convenience function for searching
  GiantSyn* find_GiantSyn(enum Cell_Type cell, enum Syn_Type syn);
  //note that after calling this, syn_info parameter is already freed from memory
  Syn* initialize_synapses(Syn_Info *syn_info);
};





///////////////////////////////////////////////////////////////=================================
//we now define several classes that are derived from CellSyn.
class REsyn:public CellSyn{

 public:

 REsyn(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){
    base_cell = new RE();
    //type = E_RE;
    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    ((RE*)base_cell)->E_l = -77; 
    ((RE*)base_cell)->G_Ca = 2.2; //2.3; 
    ((RE*)base_cell)->G_kl = 0.012; //0.015; //0.005;
    this->AMPA_Rsum=0;
    this->ismap=0;
    //double R = 2.0 * zero_one_rand() - 1.0;
    //((RE*)base_cell)->G_kl = ((RE*)base_cell)->G_kl +R * 0.001; //0.015; //0.02
  }

  void memory(Syn_Info **syn_info,int num_syns){

    ((RE*)base_cell)->init(y);

    int i = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){
      syns[i] = initialize_synapses(syn_info[i]);
      syns[i]->strength = syns[i]->strength/((RE*)base_cell)->S_RE;

      if(syns[i]->type == E_GiantAMPAMap){
        ((GiantSynAMPA*) syns[i]->giant) -> alpha = 1;
      }
      if(syns[i]->to_type == E_RE){
        ((GABA_A*)syns[i])->E_GABA = -70;
      }
    } 

    //g_a[k]->E_GABA = -70;  //GABA-A from RE to TC has more neg.revers. 
    //g_a[k].Alpha = 20;
    //g_a[k].Beta = 0.162;

  }

};

//-------------------TC CELL core type and all Synapses from other cells-----------------------------------------------
class TCcore: public CellSyn {

  //double g_a1_I, g_b_I, a_cx_tc_I, a_cx6_tc_I;
  //AMPA   *a_cx_tc, *a_cx6_tc;
  //GABA_A *g_a1;
  //GB *g_b;

 public:

 TCcore(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){
    base_cell = new TC();
    //type = E_TC;
    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=0;
    this->AMPA_Rsum=0;

    ((TC*)base_cell)->G_A = 0;
    ((TC*)base_cell)->ginc = 2.0; // TB 1.5; //2.4; // 0.9; //2; //SWS-maxim
    //((TC*)base_cell)->E_l = -70;
    ((TC*)base_cell)->G_Ca = 2.5; //2.7; // TB 2.2; //2.7; //2.3; //2.; //SWS-maxim
   
    ((TC*)base_cell)->D = 2;
    ((TC*)base_cell)->pc = 0.007;
    ((TC*)base_cell)->k4 = 0.001;
    ((TC*)base_cell)->Vtr = -40;
    ((TC*)base_cell)->VtrK = -28;
    ((TC*)base_cell)->G_K = 12;
    ((TC*)base_cell)->G_h = 0.016;  
  
    ((TC*)base_cell)->G_kl = 0.024; //0.0142;  // TB 0.03; //0.0142;  SWS-maxim 
    ((TC*)base_cell)->DC = 0; //-0.05;
    
    //double R = 2.0 * zero_one_rand() - 1.0;
    //((TC*)base_cell)->G_kl = ((TC*)base_cell)->G_kl +R * 0.001; 
  }

  void memory(Syn_Info **syn_info,int num_syns){
    
    ((TC*)base_cell)->init(y);
    int i = 0;
    int k = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);

      if(syns[i]->type == E_AMPA){
        /* syns[i] = new AMPA(syn_info[i]); */
      }else if(syns[i]->type == E_GABA_A){
        /* syns[i] = new GABA_A(syn_info[i]); */
        ((GABA_A*)syns[i])->E_GABA = -83; //-85; //(J.Neur.1997,17(7),2348)
      }else if(syns[i]->type == E_GABA_B){
        /* syns[i] = new GB(syn_info[i]); */
        //}else if(syns[i]->type == E_GABA_B){
        ((GB*)syns[i])->K1 = 0.5; //0.09;
        ((GB*)syns[i])->K2 = 0.0012;
        ((GB*)syns[i])->K3 = 0.1; //0.18;
        ((GB*)syns[i])->K4 = 0.034;
        ((GB*)syns[i])->init(y,12+2*k);
        k = k + 1;
      }
      syns[i]->strength = syns[i]->strength/((TC*)base_cell)->S_TC;
    } 
  }
  
};  

//-------------------TC CELL matrix type (TCa) and all Synapses from other cells-----------------------------------------------
class TCmatrix: public CellSyn {

  //double g_a1_I, g_b_I, a_cx_tc_I;
  //AMPA   *a_cx_tc;
  //GABA_A *g_a1;
  //GB *g_b;

 public:

 TCmatrix(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){
    base_cell = new TC();
    //type = E_TCa;
    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=0;
    this->AMPA_Rsum=0;

    ((TC*)base_cell)->G_A = 0;
    ((TC*)base_cell)->ginc = 2; // TB 1.5; //2.4; // 0.9; //2; //SWS-maxim
    //((TC*)base_cell)->E_l = -70;
    ((TC*)base_cell)->G_Ca = 2.2; //2.7; // TB 2.2; //2.7; //2.3; //2.; //SWS-maxim
   
    ((TC*)base_cell)->D = 2;
    ((TC*)base_cell)->pc = 0.007;
    ((TC*)base_cell)->k4 = 0.001;
    ((TC*)base_cell)->Vtr = -40;
    ((TC*)base_cell)->VtrK = -28;
    ((TC*)base_cell)->G_K = 12;
    ((TC*)base_cell)->G_h = 0.015;  
  
    ((TC*)base_cell)->G_kl = 0.025; //0.0142;  // TB 0.03; //0.0142;  SWS-maxim 
    ((TC*)base_cell)->DC = 0; //-0.05;

    //double R = 2.0 * zero_one_rand() - 1.0;
    //((TC*)base_cell)->G_kl =((TC*)base_cell)->G_kl +R * 0.001; 
  }

  void memory(Syn_Info **syn_info,int num_syns){

    ((TC*)base_cell)->init(y);
    int i = 0;
    int k = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);

      if(syns[i]->type == E_GABA_A){
        /* syns[i] = new GABA_A(syn_info[i]); */
        ((GABA_A*)syns[i])->E_GABA = -83; //-85; //(J.Neur.1997,17(7),2348)
      }else if(syns[i]->type == E_GABA_B){
        /*syns[i] = new GB(syn_info[i]); */
        //}else if(syns[i]->type == E_GABA_B){
        ((GB*)syns[i])->K1 = 0.5; //0.09;
        ((GB*)syns[i])->K2 = 0.0012;
        ((GB*)syns[i])->K3 = 0.1; //0.18;
        ((GB*)syns[i])->K4 = 0.034;
        ((GB*)syns[i])->init(y,12+2*k);
        k = k + 1;
      }

      syns[i]->strength = syns[i]->strength/((TC*)base_cell)->S_TC;
    } 
  }
};  


//Below we describe first headers for specific cortical cell classes based on the HH type basic class and then we describe very similar headers for the same cortical classes but based on map type basic class. The switch I_HH select which headers are used. The functions associated with these specific cell classes are the same and are describe dbelow
//-------------------CX CELL and all Synapses from other cells-----------------------------------------------
class CXsyn: public CellSyn {

 public:
  
 CXsyn(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){


    base_cell = new CX();
    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=0;
    this->AMPA_Rsum=0;
    
    ((CX*)base_cell)->G_kl = 0.011; //0.002; //0.0025;
    ((CX*)base_cell)->E_l = -67.0; //-67;
    ((CX*)base_cell)->G_Km = 0.02; //0.020; // TB 0.01;      //0.02;    //SWS-maxim
    ((CX*)base_cell)->G_l = 0.011 + (((double) rand() / (RAND_MAX)) + 1) * 0.003;
    /* ((CX*)base_cell)->CX_DEND::G_Nap = 2.0; */
     //0.01;//     working --0.015; //1.0e3/30000; //0.06; // TB 1.0e3/30000; //0.06;  //SWS-maxim

  }
  
  void memory(Syn_Info **syn_info,int num_syns){
          
    ((CX*)base_cell)->init(y,N_DEND);

    int i = 0;
    int k = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);

      if(syns[i]->type == E_GABA_B){
        ((GB*)syns[i])->init(y,13+2*k);
        k = k + 1;
      }

      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend(); 
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();

    }
  }

  void print_c_stren(FILE *connection_record) { 
        
    double total = 0.0;
    int total_num = 0;
    int i = 0;

    for(i =0; i < num_syns; i++){
      if(syns[i]->to_type == E_CX && syns[i]->from_type == E_CX && syns[i]->type == E_AMPA_D2){
        total = total + ((AMPA_D2*)syns[i])->strength; // modify this for ampa d3 do not use strength
        total_num = total_num + 1;
      }
    }    
    double average = total / total_num;
    fprintf(connection_record,"%lf %lf %d %d",cur_time,average,m,n);
    
    for(i =0; i < num_syns; i++){
      if(syns[i]->to_type == E_CX && syns[i]->from_type == E_CX && syns[i]->type == E_AMPA_D2){ 
        fprintf(connection_record," %lf",((AMPA_D2*)syns[i])->strength); 
      }
    }    
    fprintf(connection_record,"\n");
  }
  

};  

//-------------------CXa CELL and all Synapses from other cells------------------
class CXasyn: public CellSyn {

 public:

 CXasyn(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){

    base_cell = new CX();
    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    ((CX*)base_cell)->G_kl = 0.01; //0.002; //0.0025;
    ((CX*)base_cell)->E_l = -68; //-67;
    ((CX*)base_cell)->G_Km = 0.02; //0.020; // TB 0.01;      //0.02;    //SWS-maxim
    ((CX*)base_cell)->G_l =  0.022; //1.0e3/30000; //0.06; // TB 1.0e3/30000; //0.06;  //SWS-maxim

    this->ismap=0;
    this->AMPA_Rsum=0;
  }
  
  void memory(Syn_Info **syn_info,int num_syns){
    
    ((CX*)base_cell)->init(y,N_DEND);
    
    int i = 0;
    int k = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);

      if(syns[i]->type == E_GABA_B){
        ((GB*)syns[i])->init(y,13+2*k);
        k = k + 1;
      }

      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend(); 
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();

 
    }
  }


};  

//CX6 layer
//-------------------Cx6 CELL and all Synapses from other cells-----------------------------------------------
class CXsyn6: public CellSyn {

 public:

 CXsyn6(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){


    base_cell = new CX();
    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);

    ((CX*)base_cell)->G_kl = 0.01; //0.002; //0.0025;
    ((CX*)base_cell)->E_l = -68; //-67;
    ((CX*)base_cell)->G_Km = 0.02; //0.020; // TB 0.01;      //0.02;    //SWS-maxim
    ((CX*)base_cell)->G_l =  0.022; //1.0e3/30000; //0.06; // TB 1.0e3/30000; //0.06;  //SWS-maxim
    this->ismap=0;
    this->AMPA_Rsum=0;
  }
  
  void memory(Syn_Info **syn_info,int num_syns){
    
    ((CX*)base_cell)->init(y,N_DEND);

    int i = 0;
    int k = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){
      syns[i] = initialize_synapses(syn_info[i]);

      if(syns[i]->type == E_GABA_B){
        ((GB*)syns[i])->init(y,13+2*k);
        k = k + 1;
      }

      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend(); 
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();


    }
  }
};  

// TB
//-------------------IN CELL core type and all Synapses from other cells-----------------------------------------------
class INsynCore: public CellSyn {
  
 public:  

 INsynCore(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){


    base_cell = new CX();
    ((CX*)base_cell)->rho = 50;         
    ((CX*)base_cell)->CX_DEND::G_Nap = 0.0;
    ((CX*)base_cell)->CX_SOMA::G_Nap = 0.0;
    ((CX*)base_cell)->CX_SOMA::G_Na = 2500;
    //((CX*)base_cell)->E_l = -75; //-67;
    //((CX*)base_cell)->G_Km = 0.02;  

    ((CX*)base_cell)->G_kl = 0.009; //0.002; //0.0025;
    ((CX*)base_cell)->G_l = 0.009 + (((double) rand() / (RAND_MAX)) + 1) * 0.003;

    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=0;
    this->AMPA_Rsum=0;
    
    //double R = 2.0 * zero_one_rand() - 1.0;
    //((CX*)base_cell)->CX_DEND::G_Na = ((CX*)base_cell)->CX_DEND::G_Na +R * 0.5; 
    //((CX*)base_cell)->CX_SOMA::G_Kv = ((CX*)base_cell)->CX_SOMA::G_Kv +R * 50; 
    //((CX*)base_cell)->CX_SOMA::G_Na = ((CX*)base_cell)->CX_SOMA::G_Na +R * 500; 
    //((CX*)base_cell)->E_l = ((CX*)base_cell)->E_l +R * 0.5; 


  }

  void memory(Syn_Info **syn_info,int num_syns){
    
    ((CX*)base_cell)->init(y,N_DEND);

    int i = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();
      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend();



    }
  }
};  

//-------------------IN CELL matrix type (INa) and all Synapses from other cells---------
class INsynMatrix: public CellSyn {

 public:

 INsynMatrix(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){

    base_cell = new CX();
    ((CX*)base_cell)->rho = 50;         
    ((CX*)base_cell)->CX_DEND::G_Nap = 0.0;
    ((CX*)base_cell)->CX_SOMA::G_Nap = 0.0;
    ((CX*)base_cell)->CX_SOMA::G_Na = 2500;

    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=0;
    this->AMPA_Rsum=0;
    
    //double R = 2.0 * zero_one_rand() - 1.0;
    //((CX*)base_cell)->CX_DEND::G_Na = ((CX*)base_cell)->CX_DEND::G_Na +R * 0.5; 
    //((CX*)base_cell)->CX_SOMA::G_Kv = ((CX*)base_cell)->CX_SOMA::G_Kv +R * 50; 
    //((CX*)base_cell)->CX_SOMA::G_Na = ((CX*)base_cell)->CX_SOMA::G_Na +R * 500; 
    //((CX*)base_cell)->E_l = ((CX*)base_cell)->E_l +R * 0.5; 

  }
  
  void memory(Syn_Info **syn_info,int num_syns){
    
    ((CX*)base_cell)->init(y,N_DEND);

    int i = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){
      syns[i] = initialize_synapses(syn_info[i]);
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();
      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend();

    }
  }
}; 

//-------------------IN CELL matrix type (INa) and all Synapses from other cells--------
class INsyn6: public CellSyn {
  
 public:

 INsyn6(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){


    base_cell = new CX();

    ((CX*)base_cell)->rho = 50;         
    ((CX*)base_cell)->CX_DEND::G_Nap = 0.0;
    ((CX*)base_cell)->CX_SOMA::G_Nap = 0.0;
    ((CX*)base_cell)->CX_SOMA::G_Na = 2500;
    //this line was not here unit matched with incore
    //((CX*)base_cell)->E_l = -75; //-67; 

    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=0;
    this->AMPA_Rsum=0;
    
    //double R = 2.0 * zero_one_rand() - 1.0;
    //((CX*)base_cell)->CX_DEND::G_Na = ((CX*)base_cell)->CX_DEND::G_Na +R * 0.5; 
    //((CX*)base_cell)->CX_SOMA::G_Kv = ((CX*)base_cell)->CX_SOMA::G_Kv +R * 50; 
    //((CX*)base_cell)->CX_SOMA::G_Na = ((CX*)base_cell)->CX_SOMA::G_Na +R * 500; 
    //((CX*)base_cell)->E_l = ((CX*)base_cell)->E_l +R * 0.5; 

  }
  
  void memory(Syn_Info **syn_info,int num_syns){
    
    ((CX*)base_cell)->init(y,N_DEND);
    
    int i = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();
      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend();

    }
  }
}; 


//------------------ Map cells---------------------------------------
class CXsyn_Map: public CellSyn {

 public:
   
 CXsyn_Map(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){

    base_cell = new RS();

    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=1;
		
    //double R = 2.0 * zero_one_rand() - 1.0;
    //((RS*)base_cell)->sigma = ((RS*)base_cell)->sigma + R * 0.02;

  }
  
  void memory(Syn_Info **syn_info,int num_syns){
          
    ((RS*)base_cell)->init();

    int i = 0;
    /* int k = 0; */
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();
      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend();

    }
  }

  void print_c_stren(FILE *connection_record){
        
    double total = 0.0;
    int total_num = 0;
    int i = 0;

    for(i =0; i < num_syns; i++){
      if(syns[i]->to_type == E_CX && syns[i]->from_type == E_CX && syns[i]->type == E_AMPA_D2){
        total = total + ((AMPA_D2*)syns[i])->strength; // modify this for ampa d3 do not use strength
        total_num = total_num + 1;
      }
    }    
    double average = total / total_num;
    fprintf(connection_record,"%lf %lf %d %d",cur_time,average,m,n);
    
    for(i =0; i < num_syns; i++){
      if(syns[i]->to_type == E_CX && syns[i]->from_type == E_CX && syns[i]->type == E_AMPA_D2){
        fprintf(connection_record," %lf",((AMPA_D2*)syns[i])->strength);
      }
    }    
    fprintf(connection_record,"\n");
  }
  
};  

class CXasyn_Map: public CellSyn {

 public:

 CXasyn_Map(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){

    base_cell = new RS();
    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=1;

    //double R = 2.0 * zero_one_rand() - 1.0;
    //((RS*)base_cell)->sigma = ((RS*)base_cell)->sigma + R * 0.02;
  }
  
  void memory(Syn_Info **syn_info,int num_syns){
    
    ((RS*)base_cell)->init();
    
    int i = 0;
    /* int k = 0; */
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();
      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend();

    }
  }

};  

class CXsyn6_Map: public CellSyn {

 public:

 CXsyn6_Map(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){

    base_cell = new RS();
    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=1;

    //double R = 2.0 * zero_one_rand() - 1.0;
    //((RS*)base_cell)->sigma = ((RS*)base_cell)->sigma + R * 0.02;

  }
  
  void memory(Syn_Info **syn_info,int num_syns){
    
    ((RS*)base_cell)->init();

    int i = 0;
    /* int k = 0; */
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){


      syns[i] = initialize_synapses(syn_info[i]);
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();
      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend();


    }
  }
};  

class INsynCore_Map: public CellSyn {
  
 public:  

 INsynCore_Map(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){

    base_cell = new FS1();
    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=1;
    
    //double R = 2.0 * zero_one_rand() - 1.0;
    //((FS1*)base_cell)->y = ((FS1*)base_cell)->y + R * 0.001;
  }

  void memory(Syn_Info **syn_info,int num_syns){
    
    ((FS1*)base_cell)->init();

    int i = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();
      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend();


    }
  }
};  

class INsynMatrix_Map: public CellSyn {

 public:

 INsynMatrix_Map(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){


    base_cell = new FS1();

    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=1;

    
    //double R = 2.0 * zero_one_rand() - 1.0;
    //((FS1*)base_cell)->y = ((FS1*)base_cell)->y + R * 0.001;

  }
  
  void memory(Syn_Info **syn_info,int num_syns){
    
    ((FS1*)base_cell)->init();

    int i = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();
      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend();

    }
  }
}; 

class INsyn6_Map: public CellSyn {
  
 public:

 INsyn6_Map(int i, int j,int b_size, int t_size,Cell_Info *cell_info,int type):CellSyn(i,j,b_size,t_size,cell_info,type){

    base_cell = new FS1();

    this->old_v = base_cell->get_v_soma();
    memory(cell_info->syn_info,num_syns);
    this->ismap=1;

    //double R = 2.0 * zero_one_rand() - 1.0;
    //((FS1*)base_cell)->y = ((FS1*)base_cell)->y + R * 0.001;

  }
  
  void memory(Syn_Info **syn_info,int num_syns){
    
    ((FS1*)base_cell)->init();
    
    int i = 0;
    syns = new Syn*[num_syns];
    for(i =0; i < num_syns; i++){

      syns[i] = initialize_synapses(syn_info[i]);  
      syns[i]->strength = syns[i]->strength / ((OtherBase*)(this->base_cell))->get_s_cx_dend();
      syns[i]->mini_s = syns[i]->mini_s / ((OtherBase*)(this->base_cell))->get_s_cx_dend();
    
    }
  }
}; 


//peters funcs
//Peter - fix all of these to make them match also clean up arguments
int get_cell_index(int type, int m, int n);
void apply_field(double ***field_effect, double **cx5_local_field, double **cx6_local_field,double **cx5_soma_local_field, double **cx6_soma_local_field,int ii,double time,FILE **field_file);
void allocate_state_save(double ***cx_base_v_SOMA, double ***cxa_base_v_SOMA, double ***cx5_local_field,double ***cx6_local_field,double ***cx5_soma_local_field,double ***cx6_soma_local_field, double ****field_effect);
void spike_fre_calc(int *total_region,double *frequencies);
void boost_activity();
void start_critical();
void end_critical();
void root_critical();
int receive_spike(int m, int n, enum Cell_Type type);
int receive_spike_for_maps(int m, int n, enum Cell_Type type);
double receive_dend(int my_type,int m, int n, enum Cell_Type type);
extern int **max_connect;
extern int **max_connect_gap;

#endif //CellSyn_h
