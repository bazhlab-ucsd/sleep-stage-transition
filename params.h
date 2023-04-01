/*v 0.1
simple parameters parsing helper


Example of code usage:

  double param1; add_double_param(param1);		//add "param1" to be recognized as input parameter
  int param2;    add_int_param(param2);
  ...
  assert(load_parameters("input_parameters.txt"));	//read parameters from input file
  assert(cmdline_parameters(argc,argv));		//can be overriden by commandline, eg. ./simul [any_pars_without_minus] -param1 value1 -param3 value3
  print_parameters();					//dump table of parameters values


Configuration file example (input_parameters.txt):
------
#any comment starts with '#'

 #any numbers of blanks around

param1 3.14
#what is the meaning of life?
param2 42
-----

*/

#ifndef PARAMSHELP
#define PARAMSHELP

#include <string>
#include <sstream>
#include <assert.h>
#include <stdio.h>

using std::string;
using std::stringstream;

#define MAX_PARAMS 100

//buffer for parameter names and pointers to
//variables used for them in the main program
double *double_pars[MAX_PARAMS];
string double_par_name[MAX_PARAMS];
int double_count=0;
#define add_double_param(par) \
assert(double_count!=MAX_PARAMS);\
double_pars[double_count]=&par;\
/*(*double_pars[double_count])=0;*/\
double_par_name[double_count++]=#par; 

int *int_pars[MAX_PARAMS];
string int_par_name[MAX_PARAMS];
int int_count=0;
#define add_int_param(par) \
assert(int_count!=MAX_PARAMS);\
int_pars[int_count]=&par;\
/*(*int_pars[int_count])=0;*/\
int_par_name[int_count++]=#par; 

//parse one parameter on the given config line; returns 1 if line format is recognized
//result is put directly into global parameter
int parse_param(string line){
  if (line.empty()) return 1;

  string par1,par2;		//line should be pair of inputs
  stringstream parser(line);
  
  if (!(parser>>par1)) return 1; 	//whitespace line
  if (par1[0]=='#') return 1; 		//just comment
  if (par1[0]=='-') par1=par1.substr(1); //allow "-param" in commandline 
  if (!(parser>>par2)) return 0; 	//missing value

  int i;
  for (i=0; i<double_count; i++) //double parameters
    if (par1==double_par_name[i]){
      if (!(stringstream(par2)>>*(double_pars[i]))) return 0;
      return 1;
    }

  for (i=0; i<int_count; i++) //integer parameters
    if (par1==int_par_name[i]){
      if (!(stringstream(par2)>>*(int_pars[i]))) return 0;
      return 1;
    }

  fprintf(stderr, "Unknown input parameter found.\n");
  return 0; //no recognized parameter found
}//parse

//load parameters from file, returns 0 on error
int load_parameters(char const *file){
  FILE *f=fopen(file,"r");
  if (f==NULL) { fprintf(stderr, "File with input parameters not found.\n"); return 0; }

  char line[4096];
  while ( fgets(line,sizeof(line),f) != NULL )
    if (!parse_param(line)) return 0;

  fclose(f);  
  return 1;
}//load_parameters

//load parameters from commandline ("-param value" vector)
//expecting standard argc, argv arguments
//we allow different parameters on cmdline when they don't use "-" (eg. ./simul input.txt -param val)
int cmdline_parameters(int argc, char *argv[]){
  argc--;
  int i=0;
  while (argc>++i) {
    string s(argv[i]);
    if (s.empty()) continue;
    if (s[0]!='-') continue;
    if (!parse_param(s+" "+string(argv[i+1]))) return 0;
    i++; 	//jump over the value (can be negative, ie "-")
  }
 return 1;
}//cmdline_parameters

void print_parameters(){
  for (int i=0; i<int_count;i++)
    fprintf(stderr,"%s:\t%d\n",int_par_name[i].c_str(),*(int_pars[i]));

  for (int i=0; i<double_count;i++)
    fprintf(stderr,"%s:\t%lf\n",double_par_name[i].c_str(),*(double_pars[i]));
}//print_parameters

#endif
