#include "network.h"
#include <vector>
using std::vector;

//TODO SynCall touching this variable directly in calc!
//max number of connections for a cell of each type from a cell of each type not including if their are multipul connections between the same two cells
int **max_connect =NULL;
int **max_connect_gap =NULL;

void find_max_connections(int Num_Types, Pair *cell_sizes,Cell_Info ****cells_info){
  int i,j,k,l;
  max_connect = new int*[Num_Types];
  max_connect_gap = new int*[Num_Types];
  for(i=0; i< Num_Types; i++){
    max_connect[i] = new int[Num_Types];
    max_connect_gap[i] = new int[Num_Types];
    for(j=0; j< Num_Types; j++){
      max_connect[i][j] = 0;
      max_connect_gap[i][j] = 0;
      for(k = 0; k<cell_sizes[i].x; k++){
		for(l = 0; l<cell_sizes[i].y; l++){
		  if(max_connect[i][j] < cells_info[i][k][l]->num_cell_cons_norm[j]){
			max_connect[i][j] = cells_info[i][k][l]->num_cell_cons_norm[j];
		  }
		  if(max_connect_gap[i][j] < cells_info[i][k][l]->num_cell_cons_gap[j]){
			max_connect_gap[i][j] = cells_info[i][k][l]->num_cell_cons_gap[j];
		  }
		}
      }
      // printf("max connect[%d][%d]:%d \n", i,j,max_connect[i][j]);
    }
  }

  

}


void set_cell_sizes(FILE* connections_file,int &Num_Types, Pair* &cell_sizes, int* &cell_numbers){

  int x_size;
  int y_size;
  int iter = 0;

  if(fscanf(connections_file,"%d\n",&Num_Types) != 1){
    printf("bad input file could not get number of types\n");
    exit(1);
  }

  cell_sizes = new Pair[Num_Types];
  cell_numbers = new int[Num_Types];

  while(fscanf(connections_file,"%d %d\n",&x_size,&y_size) == 2){
    cell_sizes[iter].x = x_size;
    cell_sizes[iter].y = y_size;
    cell_numbers[iter] = x_size*y_size;
    iter = iter+1;
    if(iter > Num_Types){
      printf("error setting cell sizes. probably a bad input file \n");
      exit(1);
    }
  }
  return;
}


void create_connections_from_file(Cell_Info *****cells_info,Pair* cell_sizes,FILE *connections_file,int Num_Types){

  vector<vector<vector<int> > > conn_number;

  int i =0;
  int j =0;
  int k =0;
  int l = 0;

  //for debugging keeps track of how many connection between the same two cells their are in a row in an input file these repeats are ignored for purposes of calculating max connections       
  int num_repeats = 0;

  printf("Cells no of type:%d \n", Num_Types);

  conn_number.resize(Num_Types);
  for (int i = 0; i < Num_Types; ++i) {
    conn_number[i].resize(Num_Types);

    for (int j = 0; j < Num_Types; ++j)
      conn_number[i][j].resize(30);
  }
 
  //first we allocate everything
  *cells_info = new Cell_Info***[Num_Types];
  for(i=0; i<Num_Types; i++){
    (*cells_info)[i] = new Cell_Info**[cell_sizes[i].x];
    for(j=0; j<cell_sizes[i].x; j++){
      (*cells_info)[i][j] = new Cell_Info*[cell_sizes[i].y];
      for(k=0; k<cell_sizes[i].y; k++){
		(*cells_info)[i][j][k]=new Cell_Info; 
		(*cells_info)[i][j][k]->num_connects_norm = new int[Num_Types];
		(*cells_info)[i][j][k]->num_connects_gap = new int[Num_Types];
		(*cells_info)[i][j][k]->num_connects_total = new int[Num_Types];

		(*cells_info)[i][j][k]->num_cell_cons = new int[Num_Types];
		(*cells_info)[i][j][k]->num_cell_cons_gap = new int[Num_Types];
		(*cells_info)[i][j][k]->num_cell_cons_norm = new int[Num_Types];
		(*cells_info)[i][j][k]->total_connects = 0;
      }
    }
	printf("Cells in type:%d -> %d , %d \n", i,cell_sizes[i].x, cell_sizes[i].y);
  }

  printf("Cells are set \n");
  int type;
  int x_loc;
  int y_loc;

  //now we go through the connections file and fill in the info we need
  while(1){

    //get information for a cell
    // if(fscanf(connections_file,"In: %d %d %d\n",&type,&x_loc,&y_loc) != 3){
	// Incoming connections for cell: type: 0 x: 683 y: 1
    // if(fscanf(connections_file,"Incoming connections for cell: type: %d x: %d y: %d\n",&type,&x_loc,&y_loc) != 3){
 
	if(fscanf(connections_file,"In: %d %d %d\n",&type,&x_loc,&y_loc) != 3){
      break;
    }
	// Correction for 0 vs 1
	// x_loc=x_loc-1; y_loc=y_loc-1;

    fpos_t position;
    //saves where the information about this cell is located in the file because we have to make two passes through it
    fgetpos(connections_file,&position);
    int in_sizes[Num_Types];
    int i = 0;
    for(i=0; i<Num_Types; i++){
      in_sizes[i] = 0;
    }
    int in_type;
    int in_x_loc;
    int in_y_loc;

    //first pass going through all the connections to the cell
    int total_con = 0;
    while(1){
      char syn_type[250];
      double strength = 0.0;
      double mini_s = 0.0;
      double mini_fre = 0.0;
      int short_range = 0;
      //type x y Syntype strength mini_s mini_fre short_range

	  //type: 0 x: 479 y: 1 Syntype: GABA_A range: 1
      // if(fscanf(connections_file,"type: %d x: %d y: %d Syntype: %s range: %d \n",&in_type,&in_x_loc,&in_y_loc,syn_type,&short_range) != 5){
      // if(fscanf(connections_file,"%d %d %d %s \n",&in_type,&in_x_loc,&in_y_loc,syn_type) != 4){
      if(fscanf(connections_file,"%d %d %d %s %lf %lf %lf %d\n",&in_type,&in_x_loc,&in_y_loc,syn_type,&strength,&mini_s,&mini_fre,&short_range) != 8){
		// ,&strength,&mini_s,&mini_fre		
		// Range is not set 
		break;
      }
	  
	  // short_range = 1;
	  // Correction for 0 vs 1
	  // in_x_loc=in_x_loc-1; in_y_loc=in_y_loc-1;

      //we are finding out how many connections of each type it has
      in_sizes[in_type] = in_sizes[in_type] + 1;
      total_con = total_con + 1;
    }
    
	// printf("Setup for %d, %d , %d \n", type,x_loc,y_loc);

    (*cells_info)[type][x_loc][y_loc]->syn_info = new Syn_Info*[total_con];
    //allocate space for the connections of each type
    for(l=0; l<Num_Types; l++){
      //num connects set to 0 here even though we know the answer because we use it as an index for how far we have gone later in this function
      (*cells_info)[type][x_loc][y_loc]->num_connects_norm[l] = 0; 
      (*cells_info)[type][x_loc][y_loc]->num_connects_gap[l] = 0; 
      (*cells_info)[type][x_loc][y_loc]->num_connects_total[l] = 0; 

      (*cells_info)[type][x_loc][y_loc]->num_cell_cons[l] = 0; 
      (*cells_info)[type][x_loc][y_loc]->num_cell_cons_gap[l] = 0; 
      (*cells_info)[type][x_loc][y_loc]->num_cell_cons_norm[l] = 0; 
      
    }
	

    //resets the file pointer to the start of the connections for the cell we are working on
    fsetpos(connections_file,&position);

    int in_type_old = -1;
    int in_x_loc_old = -1;
    int in_y_loc_old= -1;
    //second pass going through all the connections to the cell
    while(1){
      char syn_type[250];
      double strength = 0.0;
      double mini_s = 0.0;
      double mini_fre = 0.0;
      int gap = 0;
      int short_range = 0;

      // if(fscanf(connections_file,"type: %d x: %d y: %d Syntype: %s range: %d \n",&in_type,&in_x_loc,&in_y_loc,syn_type,&short_range) != 5){
      // if(fscanf(connections_file,"%d %d %d %s\n",&in_type,&in_x_loc,&in_y_loc,syn_type) != 4){
      if(fscanf(connections_file,"%d %d %d %s %lf %lf %lf %d\n",&in_type,&in_x_loc,&in_y_loc,syn_type,&strength,&mini_s,&mini_fre,&short_range) != 8){
		break;
      }
      string syn_type_lit(syn_type);

	  // Correction for 0 vs 1
	  // in_x_loc=in_x_loc-1; in_y_loc=in_y_loc-1;

	  // Match the synapse type to set mini's and strength
	  // cell types - 0 - re, 1- tc, 2- tca, 3- cx, 4- cxa, 5- cx6, 6- in, 7- ina, 8- in6

	  

      if(syn_type_lit == "GAP"){
		gap = 1;
		// printf("type:%d , x:%d y:%d ; IN: type:%d , x:%d y:%d range: %d \n", type, x_loc, y_loc, in_type, in_x_loc, in_y_loc, short_range);
      }
      
      if(in_x_loc >= cell_sizes[in_type].x || in_y_loc >= cell_sizes[in_type].y){
		printf("probably messed up connections input file \n");
		printf("type:%d , x:%d y:%d ; IN: type:%d , x:%d y:%d range: %d \n", type, x_loc, y_loc, in_type, in_x_loc, in_y_loc, short_range);
		printf("Size: type:%d , x:%d y:%d ; IN: type:%d , x:%d y:%d range: %d \n", type, x_loc, y_loc, in_type, cell_sizes[in_type].x, cell_sizes[in_type].y, short_range);
		exit(1);
      }
      
      (*cells_info)[type][x_loc][y_loc]->num_connects_total[in_type]++;
      if(gap ==1){
		(*cells_info)[type][x_loc][y_loc]->num_connects_gap[in_type]++;
      }else{
		(*cells_info)[type][x_loc][y_loc]->num_connects_norm[in_type]++;
      }
      
      if(in_type == in_type_old && in_x_loc == in_x_loc_old && in_y_loc == in_y_loc_old){
		num_repeats = num_repeats + 1;
		//printf("num_repeats: %d \n",num_repeats);
      }else{
		num_repeats = 0;
		//printf("norm\n");
		(*cells_info)[type][x_loc][y_loc]->num_cell_cons[in_type]++;

		if(short_range == 1){ //only count short range connections for connection scaling
		  if(gap ==1){
			(*cells_info)[type][x_loc][y_loc]->num_cell_cons_gap[in_type]++;
		  }else{
			(*cells_info)[type][x_loc][y_loc]->num_cell_cons_norm[in_type]++;
		  }
		}
		in_type_old = in_type;
		in_x_loc_old = in_x_loc;
		in_y_loc_old = in_y_loc;
      }

      int total_connects = (*cells_info)[type][x_loc][y_loc]->total_connects;
      (*cells_info)[type][x_loc][y_loc]->syn_info[total_connects] = new Syn_Info;
      (*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->strength = strength;
      (*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->mini_s = mini_s;
      (*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->mini_fre = mini_fre;
      (*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->to_type = (enum Cell_Type)type;
      (*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->from_type = (enum Cell_Type)in_type;
      (*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->to_x = x_loc;
      (*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->to_y = y_loc;
      (*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->from_x = in_x_loc;
      (*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->from_y = in_y_loc;

      if(syn_type_lit == "AMPA"){
		(*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->type = E_AMPA;
      }else if(syn_type_lit == "GABA_A"){
		(*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->type = E_GABA_A;
      }else if(syn_type_lit == "GABA_B"){
		(*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->type = E_GABA_B;
      }else if(syn_type_lit == "AMPA_D2"){
		(*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->type = E_AMPA_D2;
      }else if(syn_type_lit == "NMDA_D1"){
		(*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->type = E_NMDA_D1;
      }else if(syn_type_lit == "GABA_A_D2"){
		(*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->type = E_GABA_A_D2;
      }else if(syn_type_lit == "AMPA_D1"){
		(*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->type = E_AMPA_D1;
      }else if(syn_type_lit == "AMPA_D3"){
		(*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->type = E_AMPA_D3;
      }else if(syn_type_lit == "GAP"){ 
		(*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->type = E_GAP;
      }else if(syn_type_lit == "AMPAMap_D1"){ 
		(*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->type = E_AMPAMap_D1;
      }else if(syn_type_lit == "GiantAMPA"){
		(*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->type = E_GiantAMPA;
      }else if(syn_type_lit == "GiantAMPAMap"){
		(*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->type = E_GiantAMPAMap;
      }else if(syn_type_lit == "NMDAMap_D1"){
		(*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->type = E_NMDAMap_D1;
      }else if(syn_type_lit == "GABAMap_A_D1"){
		(*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->type = E_GABAMap_A_D2;
      }else{
		printf("could not find syn type %s requested \n", syn_type_lit.c_str());
		exit(1);
      }

      (*cells_info)[type][x_loc][y_loc]->total_connects++;
	  conn_number[type][in_type][(*cells_info)[type][x_loc][y_loc]->syn_info[total_connects]->type]++;
    }
  }
  fclose(connections_file);
  find_max_connections(Num_Types,cell_sizes,*cells_info);

  printf("All Connections are set\n");
  
  for(unsigned int i=0; i<conn_number.size(); i++) {
      for(unsigned int j=0; j<conn_number[i].size(); j++) {
          for(unsigned int k=0; k<conn_number[i][j].size(); k++){
              if (conn_number[i][j][k]>0) 
              {
                  printf("No of connections for %d from %d type %d = %d \n", i,j,k,conn_number[i][j][k]);
              }
          }
      }
  }

}


void Homeostasis::boost_activity(Pair *cell_sizes, CellSyn **cells, int num_cells){

  if (boost!=1) return;
  int i = 0;
  int j = 0;
  //printf("boosted boost_syn_left:%d boost_syn_right:%d\n",boost_syn_left,boost_syn_right); 
  for(i=0; i<num_cells; i++){
	if(cells[i]->type == E_CX || cells[i]->type == E_CXa || cells[i]->type == E_CX6 ){
	  for(j = 0; j<num_regions; j++){
		if( (cells[i]->m > int(cell_sizes[cells[i]->type].x * cut_locs[j])) && (cells[i]->m < int(cell_sizes[cells[i]->type].x * cut_locs[j+1]))  ){
		  int k = 0;
		  for(k=0; k < cells[i]->num_syns; ++k){
			if(cells[i]->syns[k]->from_type ==  cells[i]->syns[k]->to_type && (cells[i]->syns[k]->type == E_AMPA_D2 || cells[i]->syns[k]->type == E_AMPA_D3)){
			  cells[i]->syns[k]->strength= cells[i]->syns[k]->strength + cells[i]->syns[k]->strength*con_boost*boost_syn[j];
			  cells[i]->syns[k]->mini_s = cells[i]->syns[k]->mini_s+ cells[i]->syns[k]->mini_s*amp_boost * boost_syn[j];
			  cells[i]->syns[k]->mini_fre = cells[i]->syns[k]->mini_fre- cells[i]->syns[k]->mini_fre*fre_boost * boost_syn[j];
			}
		  }
		}
	  }
	}
  }
}

int fre_receive(int m , int n, enum Cell_Type type);

void Homeostasis::spike_fre_calc(Pair *cell_sizes, int *cell_numbers){
  printf("calculating spike frequency \n");
  int total = 0;  
  int i = 0;
  int j = 0;
  int k = 0;


  for(k =0; k<num_regions; k++){
    total_region[k] = 0;
  }
  
  printf("adding up spikes \n");
  for(i = 0; i < cell_sizes[E_CX].x; ++i){
    for(j = 0; j < cell_sizes[E_CX].y; ++j){
      int cur_fre = fre_receive(i,j,E_CX);
      for(k=0; k<num_regions; k++){
		if((i  > int(cell_sizes[E_CX].x * cut_locs[k])) && (i < int(cell_sizes[E_CX].x * cut_locs[k+1])) ){  
		  total_region[k] = cur_fre + total_region[k];
		}
      }
      total = total + cur_fre;
    }
  }
    
  for(i = 0; i < cell_sizes[E_CXa].x; ++i){
    for(j = 0; j < cell_sizes[E_CXa].y; ++j){
      int cur_fre = fre_receive(i,j,E_CXa);
      for(k=0; k<num_regions; k++){
		if((i  > int(cell_sizes[E_CXa].x * cut_locs[k])) && ( i < int(cell_sizes[E_CXa].x * cut_locs[k+1])) ){  
		  total_region[k] = cur_fre + total_region[k];
		}
      }
      total = total + cur_fre;
    }
  }

  for(i = 0; i < cell_sizes[E_CX6].x; ++i){
    for(j = 0; j < cell_sizes[E_CX6].y; ++j){
      int cur_fre = fre_receive(i,j,E_CX6);
      for(k=0; k<num_regions; k++){
		if( (i > int(cell_sizes[E_CX6].x * cut_locs[k])) && (i < int(cell_sizes[E_CX6].x * cut_locs[k+1])) ){  
		  total_region[k] = cur_fre + total_region[k];
		}
      }
      total = total + cur_fre;
    }
  }


  printf("finding frequency \n");
  int number_of_cx_cells = cell_numbers[E_CX] + cell_numbers[E_CXa] + cell_numbers[E_CX6];
  double frequency = ((float)total) / number_of_cx_cells / (((float)fre_window)/50.0);
  for(k=0; k<num_regions; k++){
    frequencies[k] =  ((float)total_region[k]) / (number_of_cx_cells/num_regions) / (((float)fre_window)/50.0);
    printf("total_freqency: %lf  region:%d frequency:%lf total:%d \n", frequency,k,frequencies[k],total_region[k]);
  }

  printf("setting boosts \n");
  for(k=0; k<num_regions; k++){
    if(frequencies[k] <target_f){
      boost_syn[k] = 1;
    }else{
      boost_syn[k] = -1;
    }
  }
}


void LocalFieldPotential::allocate_state_save(Pair * const cell_sizes){

  //initialization should be put else where needed for print
  int i =0;
  int j = 0;
  int k = 0;
  cx_base_v_SOMA = new double*[cell_sizes[E_CX].x];
  for(i=0; i <cell_sizes[E_CX].x; i++){
    cx_base_v_SOMA[i] = new double[cell_sizes[E_CX].y];
    for(j=0; j<cell_sizes[E_CX].y; j++){
      cx_base_v_SOMA[i][j] = 0;
    }
  }

  cxa_base_v_SOMA = new double*[cell_sizes[E_CXa].x];
  for(i=0; i <cell_sizes[E_CXa].x; i++){
    cxa_base_v_SOMA[i] = new double[cell_sizes[E_CXa].y];
    for(j=0; j<cell_sizes[E_CXa].y; j++){
      cxa_base_v_SOMA[i][j] = 0;
    }
  }

  field_effect = new double**[num_field_layers];
  for(i=0; i< num_field_layers; i++){
    field_effect[i] =  new double*[cell_sizes[E_CX].x];
    for(j = 0; j < cell_sizes[E_CX].x; j++){
      field_effect[i][j] =  new double[cell_sizes[E_CX].y];
      for(k = 0; k < cell_sizes[E_CX].y; k++){
        field_effect[i][j][k] = 0;
      }
    }
  }   

  cx5_local_field = new double*[cell_sizes[E_CX].x];
  cx6_local_field = new double*[cell_sizes[E_CX].x]; 
  cx5_soma_local_field = new double*[cell_sizes[E_CX].x];
  cx6_soma_local_field = new double*[cell_sizes[E_CX].x];
  
  for(i = 0; i < cell_sizes[E_CX].x; i++){
    cx5_local_field[i] = new double[cell_sizes[E_CX].y];
    cx6_local_field[i] = new double[cell_sizes[E_CX].y];
    cx5_soma_local_field[i] = new double[cell_sizes[E_CX].y];
    cx6_soma_local_field[i] = new double[cell_sizes[E_CX].y];
    for(j = 0; j < cell_sizes[E_CX].y; j++){
      cx5_local_field[i][j] = 0;
      cx6_local_field[i][j] = 0;
      cx5_soma_local_field[i][j] = 0;
      cx6_soma_local_field[i][j] = 0;
    }
  }
}


//calculates distance between 2 points in 3 space. Not this is non standard as in quadrouples distances in the z direction because layers are 4 times further apart than cells maybe
double distance(double x, double y, double z, double i, double j, double k){
  return sqrt( ((x-i)*(x-i)) + ((y-j)*(y-j)) + (((z-k)*4)*((z-k)*4)) );   
}

//TODO kill this extern, perhaps new class needed
extern CellSyn **cells;
int get_cell_index(int type, int m, int n);

//dendrite part
double get_cx5_local_field(int i, int j){

  double message = 0;
  int group_index = get_cell_index(E_CXa,i,j);
  message =  -1.0*((CX*)((CXasyn*)cells[group_index])->base_cell)->axil_current;
  return message;
}

//soma part
double get_cx5_soma_local_field(int i, int j){
  
  double message = 0;
  int group_index = get_cell_index(E_CXa,i,j);
  message =  1.0*((CX*)((CXasyn*)cells[group_index])->base_cell)->axil_current;
  return message;
}

double get_cx6_soma_local_field(int i, int j){

  double message = 0;  
  int group_index = get_cell_index(E_CX6,i,j);
  message =  1.0*((CX*)((CXsyn6*)cells[group_index])->base_cell)->axil_current;
  return message;
}

double get_cx6_local_field(int i, int j){

  double message = 0;
  int group_index = get_cell_index(E_CX6,i,j);
  message =  -1.0*((CX*)((CXsyn6*)cells[group_index])->base_cell)->axil_current;	
  return message;
}


//local field effect  probably only works 1D. And certainly only works if CX and CXa and cx6 are the same size
void LocalFieldPotential::apply_field(int ii,double time,double ttime, double t, Pair *cell_sizes){
  int i,j,k,q,z;
  
  //gets information about the field cells are generating.(this is currently current:)
  for(i = 0; i < cell_sizes[E_CX].x; ++i){
    for(j = 0; j < cell_sizes[E_CX].y; ++j){
      cx6_local_field[i][j] = get_cx6_local_field(i,j);
      cx5_local_field[i][j] = get_cx5_local_field(i,j);
      cx6_soma_local_field[i][j] = get_cx6_soma_local_field(i,j);
      cx5_soma_local_field[i][j] = get_cx5_soma_local_field(i,j);
    }
  }
  
  //relative location of each layer in space
  //cx6 dend 4 soma 6
  //cx5 dend 1 soma 5
  //cx  dend 4 soma 4

  for(z=0; z<num_field_layers; z++){
    for(k=0; k<cell_sizes[E_CX].x; k++){
      for(q=0; q<cell_sizes[E_CX].y; q++){
        for(i = 0; i < cell_sizes[E_CX].x; i++){
		  for(j = 0; j < cell_sizes[E_CX].y; j++){
            if(i!=k || j!=q || z!= 4){
              field_effect[z][k][q] = field_effect[z][k][q] + (cx6_local_field[i][j] / distance(k,q,z,i,j,4));
            }
            if(i!=k || j!=q || z!= 1){
              field_effect[z][k][q] = field_effect[z][k][q] + (cx5_local_field[i][j] / distance(k,q,z,i,j,1));
            }
            if(i!=k || j!=q || z!= 6){
              field_effect[z][k][q] = field_effect[z][k][q] + (cx6_soma_local_field[i][j] / distance(k,q,z,i,j,6));
            }
            if(i!=k || j!=q || z!= 5){
              field_effect[z][k][q] = field_effect[z][k][q] + (cx5_soma_local_field[i][j] / distance(k,q,z,i,j,5));
            }
          }
        }
      }
    }
  }
  
  //apply the calculated field effects
  for(i = 0; i < cell_sizes[E_CX].x; ++i){  
    for(j = 0; j < cell_sizes[E_CX].y; ++j){
	  int group_index = get_cell_index(E_CX,i,j);
	  ((CX*)((CXsyn*)cells[group_index])->base_cell)->field_effect_dend =  field_effect[4][i][j]*lfp_scale;
	  ((CX*)((CXsyn*)cells[group_index])->base_cell)->field_effect_soma =  field_effect[4][i][j]*lfp_scale;

	  group_index = get_cell_index(E_CXa,i,j);
	  ((CX*)((CXasyn*)cells[group_index])->base_cell)->field_effect_dend = field_effect[1][i][j]*lfp_scale;
	  ((CX*)((CXasyn*)cells[group_index])->base_cell)->field_effect_soma = field_effect[5][i][j]*lfp_scale;

	  group_index = get_cell_index(E_CX6,i,j);
	  ((CX*)((CXasyn*)cells[group_index])->base_cell)->field_effect_dend = field_effect[4][i][j]*lfp_scale;
	  ((CX*)((CXasyn*)cells[group_index])->base_cell)->field_effect_soma = field_effect[6][i][j]*lfp_scale;
    }
  }
  
  //print all the field effect information on multipules of 50
  if((t > ttime) && ((ii/(50))*(50) == ii)){ 
    for(i=0; i<num_field_layers; i++){
      fprintf(field_file[i], "%lf ", t);
      double total = 0.0;

      for(k=0; k<cell_sizes[E_CX].x; k++){
        for(q=0; q<cell_sizes[E_CX].y; q++){
		  total = total + field_effect[i][k][q];
        }
      }

      double average = total/ (cell_sizes[E_CX].x * cell_sizes[E_CX].y);
      fprintf(field_file[i],"%.15le ", average*lfp_scale);

      for(k=0; k<cell_sizes[E_CX].x; k++){
        for(q=0; q<cell_sizes[E_CX].y; q++){
          fprintf(field_file[i],"%.15le ",field_effect[i][k][q]*lfp_scale);
        }
      }

      fprintf(field_file[i],"\n");
    }
  }
  
  //reset counters for field effects for next round(probably should be local variables so we do not have to worry about this)

  for(i=0; i<num_field_layers; i++){
    for(k=0; k<cell_sizes[E_CX].x; k++){
      for(q=0; q<cell_sizes[E_CX].y; q++){
        field_effect[i][k][q] = 0.0;  
      }
    }
  }      
}

