#include <string>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm>

using namespace std;
ofstream * ConnSummaryFile;

typedef unsigned char BYTE;
//default_random_engine generator;



#define MAXTYPES 33
class CellTypes {
 public:
  struct size { int x, y; string name;};
  int n; 		//number of all types we have
  size cell[MAXTYPES]; 	//info for each type
  int get_type(const string &t){  //map name into index
    for (int i=0; i<n; i++) if (cell[i].name==t) return i;
    assert(0); //no such neuron type exists
  }//get_type
  void dump(){ for (int i=0; i<n; i++) {dump(i);cerr<<"\n";} }
  void dump(int i){
     assert(i<n);
     cerr<<cell[i].name<<" "<<cell[i].x<<" "<<cell[i].y;
  }
  void output(){ //fixed output format used by simulation
    cout<<n<<"\n";
    for (int i=0; i<n; i++)
     cout<<cell[i].x<<" "<<cell[i].y<<"\n";
  }
  CellTypes():n(0){};
};//CellTypes

CellTypes CT;

//connection types (not particular connections between neurons)
class Connections{
 public:
  //connection types between layers
  struct conn { 
    int from, to; string synapse;
      double radius_min, radius_max, cs; //radius; column size
      double probab, probab_oc; // probability within  and outside of column
    string distribution; 
    double strength,
    mini_strength,   //minimal strength of miniature synaptic potentials
    mini_freq;       //frequency of spontaneous minis (miniature synaptic potentials)
    int range;       //long/short range; this value has different meaning for generated (1==short range) and MRI network
   };//conn
  int n; //number of connections
  conn C[MAXTYPES*MAXTYPES]; //connections
  
  int find(int from,int to,string synapse,int range_type){
    for (int i=0; i<n; i++)
      if (C[i].from == from && C[i].to == to && C[i].synapse == synapse && range_type == C[i].range)
        return i;
   cerr<<"Missing connection "<<from<<"->"<<to<<" "<<synapse<<", range "<<range_type<<"\n";
   //exit(0); //connection type not entered in config file
   return 100;
  }//find 
  void dump(){ for (int i=0; i<n; i++) dump(i); }//dump
  void dump(int i){
     assert(i<n);
     cerr<<C[i].from<<" "<<C[i].to<<" "<<C[i].synapse<<" "<<
     C[i].radius_min<<" "<<C[i].radius_max<<" "<<C[i].probab<<"\n";
  }
  int stats[MAXTYPES*MAXTYPES]; //connection statistics for debug purposes
  void dump_stats(){
    for (int c=0; c<n; c++) {
      cerr<<c+1<<" of "<<n<<" , edges: "<<stats[c]<<" \ttype: (";
      CT.dump(C[c].from);cerr<<") -> (";CT.dump(C[c].to);cerr<<") ";
      dump(c);
    }
  }//dump stats
  Connections():n(0){ for (int i=0;i<MAXTYPES*MAXTYPES;i++) stats[i]=0;};
};//Connections

Connections CN;

//C++11 has std::round
double round(double x){
  return (x<0.0) ? ceil(x-0.5f) : floor(x+0.5f);
}

#define sqr(x) ((x)*(x))
//euclidean distance; we could use manhattan for 2D networks actually
//check carefuly comments in generate_connections routine in case
//you want to change this metric.
inline double dist(double x1, double y1, double x2, double y2){
  return sqrt(sqr(x1-x2)+sqr(y1-y2));
}

inline int is_column(int x1, int y1, int x2, int y2, int CS){
    if ((x1/CS == x2/CS) && (y1/CS==y2/CS)) return 1;
    else return 0;
}


//TODO use mersenne twister instead of this weak guy
inline double rand01(){return ((double) rand() / (RAND_MAX));}

string get_input_line(FILE *f,string &line){
  char buf[4096];
  string tok1;

  while(1) {
    if (!fgets(buf,sizeof(buf),f)) {return line="";}; //EOF
    line=buf;
    if (line.empty()) continue;
    stringstream parser(buf);
    if (!(parser>>tok1)) continue;    //whitespace line
    if (tok1[0]=='#') continue;       //just comment
    break;
  }
  return tok1;
}

void parse_config(char const *file){
  string line,tok;
  double x,y,ratio_x,ratio_y;
  FILE *f=fopen(file,"r");
  assert(f);

  // Types
  tok=get_input_line(f,line);
  assert(tok=="[Types]");
  //ratios
  tok=get_input_line(f,line);
  stringstream parser(line);
  assert(parser >> ratio_x); assert(parser >> ratio_y);
  //types enumeration
  tok=get_input_line(f,line);
  while (tok!="[Connections]"){
    assert(CT.n<MAXTYPES);
    parser.str(line);
    assert(parser >> CT.cell[CT.n].name >> x >> y);

    CT.cell[CT.n].x=round(x*ratio_x); CT.cell[CT.n].y=round(y*ratio_y);
    CT.n++;

    tok=get_input_line(f,line);
  }//types

  //Short range connections
  string from, to, syn, dist; double rad,cs,p,poc,str,mf,ms;
  assert(tok=="[Connections]");
  tok=get_input_line(f,line);
  while (tok!="[Longrange]"){
    assert(CN.n<MAXTYPES*MAXTYPES);
    parser.str(line);
    assert(parser >> from >> to >> syn >> rad >> cs >> p >> poc >> dist >> str >> ms >> mf);

    Connections::conn &c=CN.C[CN.n];	 //current connection
    c.from=CT.get_type(from); c.to=CT.get_type(to); c.synapse=syn;
    c.radius_min=0; c.radius_max=rad; c.cs=cs; c.probab=p; c.probab_oc=poc;
    c.distribution=dist; c.strength=str; c.mini_strength=ms; c.mini_freq=mf;
    c.range=1;
    CN.n++; 
    
    tok=get_input_line(f,line);
  }//connections

  //Long range connections
  double rad_min,rad_max;
  assert(tok=="[Longrange]");
  tok=get_input_line(f,line);
  while (!feof(f)){
    assert(CN.n<MAXTYPES*MAXTYPES);
    parser.str(line);
    assert(parser >> from >> to >> syn >> rad_min >> rad_max >> p >> dist >> str >> ms >> mf);

    Connections::conn &c=CN.C[CN.n];	 //current connection
    c.from=CT.get_type(from); c.to=CT.get_type(to); c.synapse=syn;
    c.radius_min=rad_min; c.radius_max=rad_max; c.probab=p;
    c.range=0;
    CN.n++; 

    //TODO we should change rad_min for small (1D) networks to be bigger or use specific function
    
    tok=get_input_line(f,line);
  }//connections

  fclose(f);
}

//representing single synaptic connection when generating network
class edge {
 public:
  int from_t, from_x, from_y, to_t, to_x, to_y; //out==from, in==to
  Connections::conn *syntype; 
  edge(int ft,int fx,int fy,int tt,int tx,int ty,Connections::conn *s):
    from_t(ft),from_x(fx),from_y(fy),to_t(tt),to_x(tx),to_y(ty),syntype(s){};
  static int ct,cx,cy; //cache for last to_ printed -> we want to print "In: " line sometimes
  int output(){
     int new_neuron=0; //did we enter new neurons while printing connections; for total neurons count.
     if (to_t!=ct || to_x!=cx || to_y!=cy) {
       cout<<"In: "<<to_t<<" "<<to_x<<" "<<to_y<<"\n";
       ct=to_t; cx=to_x; cy=to_y;
       new_neuron=1;
     }
     cout<<from_t<<" "<<from_x<<" "<<from_y<<" "<<syntype->synapse<<" "<<syntype->strength
         <<" "<<syntype->mini_strength<<" "<<syntype->mini_freq<<" "<<syntype->range<<"\n";
     return new_neuron;
  }
  void dump(){
     cerr<<from_t<<" "<<from_x<<" "<<from_y<< " -> " <<to_t<<" "<<to_x<<" "<<to_y<<" "<<syntype->synapse<<"\n";
  }
  inline static bool cmp(const edge &a, const edge &b){ //for sorting
     if (a.to_t!=b.to_t) return (a.to_t<b.to_t);
     if (a.to_y!=b.to_y) return (a.to_y<b.to_y);
     if (a.to_x!=b.to_x) return (a.to_x<b.to_x);
     if (a.from_t!=b.from_t) return (a.from_t<b.from_t);
     if (a.from_y!=b.from_y) return (a.from_y<b.from_y);
     if (a.from_x!=b.from_x) return (a.from_x<b.from_x);
     return false;
  }//cmp
};//edge
int edge::ct=-1; int edge::cx=-1; int edge::cy=-1;
class Edges{
 public:
  vector<edge> edges;
  //this is rough version, performance can be substantially improved in case we still need it by
  //1) testing/enumerating neurons only in radius - DONE
  //2) Using 3D fixed array for sorting results instead of 1D vector
  void generate_connections(){ //generates geometry topology
    cerr<<"No input file given, generating network on our own\n";

    for (int c=0; c<CN.n; c++){ //layer connections
      int FX = CT.cell[CN.C[c].from].x; int FY = CT.cell[CN.C[c].from].y; //input layer dimension
      int TX = CT.cell[CN.C[c].to].x; int TY = CT.cell[CN.C[c].to].y; //output layer dimension
      int ec=0;	//edges counter for statistics
      for (int fx=0; fx<FX; fx++)
        for (int fy=0; fy<FY; fy++) {
          //get scaled position of target neuron in TO layer
          int sc_x = (double)TX/FX * fx; int sc_y = (double)TY/FY * fy;
          //restrict the local neighbourhood where we look for possible connections
	  //by given radius. this speedup trick can fail in case someone switch to
	  //non-eucledian metric.
          //[max(sc_x-radius,0)..min(sc_x+radius,TX)]
	  //[max(sc_y-radius,0)..min(sc_y+radius,TY)]
          int radius = CN.C[c].radius_max + 1;

          for (int tx=max(sc_x-radius,0); tx<min(sc_x+radius,TX); tx++) //particular connections between neurons
            for (int ty=max(sc_y-radius,0); ty<min(sc_y+radius,TY); ty++) { 

             if (sc_x==tx && sc_y==ty && CT.cell[CN.C[c].from].name==CT.cell[CN.C[c].to].name) continue; //never connect neuron to itself
             //distance is measured between 'From' neuron projected into 'TO' layer
             double d=dist(sc_x,sc_y, tx,ty);

             if ( d>=CN.C[c].radius_min && d<=CN.C[c].radius_max ){
                  // check if neurons are in same or different column
                 if (CN.C[c].range==1) {
                     if (is_column(sc_x,sc_y,tx,ty,CN.C[c].cs)) {
                       if (rand01() <= CN.C[c].probab) 
                         edges.push_back(edge(CN.C[c].from,fx,fy,CN.C[c].to,tx,ty,&CN.C[c])),ec++;
                    } else {  //out of column
                        if (rand01() <= CN.C[c].probab_oc) 
                          edges.push_back(edge(CN.C[c].from,fx,fy,CN.C[c].to,tx,ty,&CN.C[c])),ec++;
                      } 
                 } else {  //longrange
                     if (rand01() <= CN.C[c].probab) 
                       edges.push_back(edge(CN.C[c].from,fx,fy,CN.C[c].to,tx,ty,&CN.C[c])),ec++;
                   }
             }//radius
            }//ty
          }//fy

      //dbg
      cerr<<c+1<<" of "<<CN.n<<" , edges: "<<ec<<" \ttype: (";
      CT.dump(CN.C[c].from);cerr<<") -> (";CT.dump(CN.C[c].to);cerr<<") ";
      CN.dump(c);
    }//c
  }//generate_connections
  void generate_3D_connections();
  void write_connections(){
   cerr<<"Writing output file\n";
   //header
   CT.output();
   int te=0,tn=0;	//total edges, total neurons
   vector<edge>::iterator it=edges.begin();
   vector<edge>::iterator end=edges.end();
   for (;it!=end;it++){
    te++; 
    tn += (*it).output(); 
    //(*it).dump();
   }
   cerr<<"Total counts: "<<tn<<" neurons, "<<te<<" edges.\n";
  }//write_connections

  //assign connection attributes for each edge in network (this can be dynamic) 
  //this routine is separated from topology generation since we load some topology-only networks 
  void apply_synaptic_types(){

   cerr<<"Sorting connections\n";
   sort(edges.begin(),edges.end(),edge::cmp);

   cerr<<"Adjusting synaptic types\n";
   vector<edge>::iterator it=edges.begin();
   vector<edge>::iterator end=edges.end();
   for (;it!=end;it++){
    edge &c=(*it);
    if (c.syntype->distribution=="fixed") 
      continue; //syntype->already contains values we will use in output

    //TODO decide what to do with gaussian or normal distribution
//    if (c.syntype->distribution=="uniform") {}// ..
//    if (c.syntype->distribution=="gauss") {}// ..
   }    
 }//apply

 //dump all output connection for a given neuron (we usually write down the opposite)
 void dump_from_neuron_edges(int type, int x, int y, int to_type){
   vector<edge>::iterator it=edges.begin();
   vector<edge>::iterator end=edges.end();
   for (;it!=end;it++){
    edge &c=(*it);
    if (c.from_t==type && c.from_x == x && c.from_y == y && (to_type==-1 || c.to_t == to_type))
      cerr<<type<<" "<<x<<","<<y<<" -> "<<c.to_t<<" "<<c.to_x<<","<<c.to_y<<"\n";
   }    
 }//dump
 void dump_from_neuron_edges(int type, int x, int y){ dump_from_neuron_edges(type, x, y, -1); }
 void load_MRI_network(const char *file);
};//Edges

Edges CE;

//Tables for mapping specific synapses to MAP types for the hybrid model
//Tables are derived from the loading code, but the whole thing can be explained
//much simpler at the end: CX & IN are Map based neurons and _any_ synapse on them
//must be converted into MAP type; lazy to simplify this now.
//The mechanism needs to be checked in case short X long range connection differ
//in synapses used.
//All to all connections, e.g. CX* -> CX*
const string rule_all[][4] =
                       { {"CX","CX",  "AMPA_D2",  "AMPAMap_D1"},
                         {"CX","CX",  "NMDA_D1",  "NMDAMap_D1"},
                         {"TC","CX",  "AMPA_D1",  "AMPAMap_D1"},
                         {"TC","IN",  "AMPA_D1",  "AMPAMap_D1"} };

//Paired connections, eg. CX/a/6 -> IN/a/6 (not CX6->INa)
const string rule_paired[][4]=
                       { {"CX","IN",  "AMPA_D2",  "AMPAMap_D1"},
                         {"CX","IN",  "NMDA_D1",  "NMDAMap_D1"},
                         {"IN","CX","GABA_A_D2","GABAMap_A_D1"},
                         {"CX","CX",  "NMDA_D1",  "NMDAMap_D1"}}; //these should be actually part of the previous

const int rules_all=sizeof(rule_all)/sizeof(rule_all[0]);
const int rules_paired=sizeof(rule_paired)/sizeof(rule_paired[0]);

//does the rule and string match in the beginning in both cases
int match(string const from, string const rule_from, string const to, string const rule_to){
  return !from.compare(0,rule_from.length(),rule_from) && !to.compare(0,rule_to.length(),rule_to);
}

//does the rule and string match in the beginning in both cases, does the suffix matches (e.g. 6 in  CX6)
int match_paired(string const from, string const rule_from, string const to, string const rule_to){ 
   assert(from.length() && to.length());
   if (from==rule_from && to==rule_to) //without suffix
     return true;
   return !from.compare(0,rule_from.length(),rule_from) &&
          !to.compare(0,rule_to.length(),rule_to) &&
	  from[from.length()-1] == to[to.length()-1];          
}

//load network from the file (fixed structure) we
//get from ucsd folks approximating MRI topology
void Edges::load_MRI_network(const char *file){
  cerr<<"Loading network\n";

  //cleanup x,y header information, we will fill it with what we read from here.
  for (int ii=0; ii<CT.n; ii++) CT.cell[ii].x=CT.cell[ii].y=0;

  FILE *f=fopen(file,"r"); assert(f);
  //now we go through the connections file and fill in the info we need
  while(1){
    int type,x_loc,y_loc;
    if(fscanf(f,"Incoming connections for cell: type: %d x: %d y: %d\n",&type,&x_loc,&y_loc) != 3)
      break;

    CT.cell[type].x = max(x_loc,CT.cell[type].x); CT.cell[type].y = max(y_loc,CT.cell[type].y);

    //second pass going through all the connections to the cell
    while(1){
      char *syn_type = new char[256];
      double strength = 0.0;
      double mini_s = 0.0;
      double mini_fre = 0.0;
      int mri_range = 0;
      int in_type,in_x_loc,in_y_loc;

      if(fscanf(f,"type: %d x: %d y: %d Syntype: %s range: %d \n",&in_type,&in_x_loc,&in_y_loc,syn_type,&mri_range) != 5)
	break;
      
      //transform synaptic type from given file to hybrid model with MAP synapse type currently used.
      string syn_type_lit(syn_type);
      string from(CT.cell[in_type].name);
      string to(CT.cell[type].name);

      //Short/long range connection conversion:
      //3 used for inter-hemispheric "short" range connections, 0 for intrahemispheric long range connections, 1 forshort range
      //3 is currently ignored, we don't have connections defined
      if (mri_range == 3) continue;
      int short_range = mri_range; //except for 3 the same meaning as in normal generation

      //transform synapse names according to the tables defined above
      for (int i=0; i<rules_all; i++){
        if (rule_all[i][2]!=syn_type_lit) continue;
	if (match(from,rule_all[i][0],to,rule_all[i][1]))
	  syn_type_lit=rule_all[i][3];
      }
      for (int i=0; i<rules_paired; i++){
        if (rule_paired[i][2]!=syn_type_lit) continue;
	if (match_paired(from,rule_paired[i][0],to,rule_paired[i][1]))
	  syn_type_lit=rule_paired[i][3];
      }

      int conn_type = CN.find(in_type,type,syn_type_lit,short_range);
      edges.push_back(edge(in_type,in_x_loc,in_y_loc,type,x_loc,y_loc,&CN.C[conn_type]));
      CN.stats[conn_type]++; //for stats printing
    }//while 2nd pass
   }//while file

  CN.dump_stats();
  fclose(f);
}//load_MRI_network

#define MRI_POINTS 12500
float Distance3D[MRI_POINTS][MRI_POINTS]; //distance lookup table given for MRI data

void load_3D_distances(const char *file){
  cerr<<"Loading 3D network distances\n";

  for (int i=0; i<MRI_POINTS; i++)
    for (int ii=0; ii<MRI_POINTS; ii++)
      Distance3D[i][ii]=-1;

  FILE *f=fopen(file,"r"); assert(f);

  //go through the distance matrix
  int i=0;
  while(1){
    int from,to;
    float dist;
    if(fscanf(f,"%d %d %f\n",&from,&to,&dist)!=3)
      break;
    from--; to--; //MATLAB vs C indexing
    assert(from<MRI_POINTS); assert(to<MRI_POINTS);
    Distance3D[from][to]=dist;
  }//while

  fclose(f);
}

#define MAX_SUBNET_POINTS 12500	//The actual number of smaller net is to be read out from network.cfg (TC cells)
float Subnet[MAX_SUBNET_POINTS]; //distance lookup table given for MRI data

//load subset of CX neurons to be flagged as indexes for TC,RE & IN
void load_3D_subnet(const char *file){
  cerr<<"Loading 3D sub network mapping\n";

  for (int i=0; i<MAX_SUBNET_POINTS; i++)
    Subnet[i]=-1;

  FILE *f=fopen(file,"r"); assert(f);

  //go through the distance matrix
  int i=0;
  while(1){
    int id,map;
    if(fscanf(f,"%d %d\n",&id,&map)!=2)
      break;
    id--; map--; //MATLAB vs C indexing
//    cerr<<id<<" "<<map<<" " << CT.cell[CT.get_type("TC")].x<<"\n";
    assert(id<MRI_POINTS); assert(map<MAX_SUBNET_POINTS); assert(id< CT.cell[CT.get_type("TC")].x);
    Subnet[id]=map;
  }//while

  fclose(f);
}


inline double dist3D(int id_from, int id_to){
 assert(Distance3D[id_from][id_to]!=-1);
 return Distance3D[id_from][id_to];
}

void Edges::generate_3D_connections(){ //generates geometry topology
    cerr<<"Generating connection file from pre-set 3D MRI data\n";

    int small = CT.cell[CT.get_type("TC")].x; int large = CT.cell[CT.get_type("CX")].x;
    cerr<<"Full set(CX): "<<large<<", Subnet(TC,RE,IN...): "<<small<<"\n";



    for (int c=0; c<CN.n; c++){ //layer connections
      int FX = CT.cell[CN.C[c].from].x; int FY = CT.cell[CN.C[c].from].y; //input layer dimension
      int TX = CT.cell[CN.C[c].to].x; int TY = CT.cell[CN.C[c].to].y; //output layer dimension
      assert(FY==1); assert(TY==1); //MRI data are stored as 1D indexes

      int ec=0;	//edges counter for statistics

      for (int fxi=0, fx= (FX==large) ? 0:Subnet[0]; (FX==large) ? fx<FX: Subnet[fxi]!=-1; (FX==large)? fx++: fx=Subnet[++fxi]) {
          //scaling is actually hidden in Subnet array, identity in large case, map in small case
	  int fy = 0;
          int sc_x = fx; int sc_y = fy;

          for (int txi=0, tx= (TX==large) ? 0:Subnet[0]; (TX==large) ? tx<TX: Subnet[txi]!=-1; (TX==large)? tx++: tx=Subnet[++txi]) { //particular connections between neurons
	    int ty=0;
	    if (sc_x==tx && sc_y==ty && CT.cell[CN.C[c].from].name==CT.cell[CN.C[c].to].name) continue; //never connect neuron to itself
             //distance is measured between 'From' neuron projected into 'TO' layer
             double d=dist3D(sc_x, tx);

             if ( d>=CN.C[c].radius_min && d<=CN.C[c].radius_max ){
                // column is no more checked
               if (CN.C[c].range==1) {
                 if (rand01() <= CN.C[c].probab_oc) 
                   edges.push_back(edge(CN.C[c].from,(FX==large)? fx:fxi,fy,CN.C[c].to,(TX==large)? tx:txi,ty,&CN.C[c])),ec++;
               } else {  //longrange
                   if (rand01() <= CN.C[c].probab) 
                       edges.push_back(edge(CN.C[c].from,(FX==large)? fx:fxi,fy,CN.C[c].to,(TX==large)? tx:txi,ty,&CN.C[c])),ec++;
                 } //longrange
             }//radius
            }//tx
          }//fx

      //dbg
      cerr<<c+1<<" of "<<CN.n<<" , edges: "<<ec<<" \ttype: (";
      CT.dump(CN.C[c].from);cerr<<") -> (";CT.dump(CN.C[c].to);cerr<<") ";
      CN.dump(c);
    }//c
  }//generate_3Dconnections



int main(int argc, char *argv[]){
	if (argc==1){
		throw std::invalid_argument("At least network configuration file should be provided");
	}
	char const *file  = argv[1];
	parse_config(file);
	//any commandline parameter taken as network connections file
	if (argc==3)
		CE.load_MRI_network(argv[2]); 		//connection file we get verbatim from ucsd group
	else if (argc==4) {				//partly processed data from ucsd folks - only subnet and distances given
		load_3D_subnet(argv[3]);
		load_3D_distances(argv[2]);
		CE.generate_3D_connections();
	} else
		CE.generate_connections(); 		//our own random connectivity

	CE.apply_synaptic_types();	//needed in case we need sample parameters for strength/mini_X parameters
	//CE.dump_from_neuron_edges(4,20,20);
	//CE.dump_from_neuron_edges(0,0,0,2);
	CE.write_connections();
	//CT.dump();
	//CN.dump();
	return 0;
}
