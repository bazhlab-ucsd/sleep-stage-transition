//Created by M.Bazhenov 1997-2009
//All rights reserved
//

#include "currents.h"
using namespace std;


//-------------------RS CELL (main class to describe regular spiking neuron)---------------

void RS::calc(double t, double I, double *y1, double *f1){

  beta_n = beta_e * I + beta_dc * Idc;
  sigma_n = sigma_e * I + sigma_dc * Idc;

  
	//     if (beta_n < -0.0001) beta_n = -0.0001;
  if (beta_n < -1.0) beta_n = -1.0;
  if (beta_n > 1.0) beta_n = 1.0;

  if(xp <= 0.0) { 
    x = alpha / (1.0 - xp) + y +beta_n;
    spike = 0;
  }else{
    if(xp <= alpha + y +beta_n && xpp <= 0.0) { 
      x = alpha + y + beta_n;
      spike = 1;
    }else {
      x = -1;
      spike = 0;
    }
  }

  y = y - mu* (xp +1.0) + mu * sigma + mu * sigma_n;
  if (y>yLimit){
    y = yLimit;
  }
  
  xpp = xp;
  xp = x;
  //     y=1;
  v_DEND = v_SOMA = x*50-15;

}



//-------------------FS1 CELL (main class to describe fast spiking interneuron)-----------------
void FS1::calc(double t, double I, double *y1, double *f1){

  if(spike > 0.1){  
    ii =  0.60 * ii - gg;
  }else{
    ii = 0.60 * ii;
  }       

  beta_n = beta_e * I + beta_dc * Idc; 
  sigma_n = sigma_e * I + sigma_dc * Idc; 
  //if (beta_n < -0.0001) beta_n = -0.0001;
  if (beta_n < -1.0) beta_n = -1.0;
  if (beta_n > 1.0) beta_n = 1.0;
  beta_n = beta_n + beta_ii *ii;
  sigma_n = sigma_n + sigma_ii *ii;
  if(xp <= 0.0) {  
    x = alpha / (1.0 - xp) + y +beta_n; 
    spike = 0; 
  }else{ 
    if(xp <= alpha + y +beta_n && xpp <= 0.0) {  
      x = alpha + y + beta_n; 
      spike = 1; 
    }else{ 
      x = -1; 
      spike = 0; 
    } 
  }

  //y = y - mu* (xp +1.0) + mu * sigma + mu * sigma_n;
  y = -2.90;
  xpp = xp;
  xp = x; 
  v_DEND = v_SOMA = x*50-15;

} 

//---------------Low-threshold Ca2+ current (RE cell)---------------------
double IT_RE::Shift = 2, IT_RE::Ca_0 = 2, IT_RE::Cels = 36;

void IT_RE::calc(double m, double h, double &fm, double &fh,
                 double v, double cai, double x) {
  ratio = Ca_0/cai;
  if(ratio <= 0.)printf("\n LOG ERROR: RE: cai=%lf ratio=%lf",cai,ratio);
  eca = eca0 * log(ratio);
  iT = G_Ca*m*m*h*(v - eca);                              
  m_inf = 1/(1 + exp(-(v + 52)/7.4));
  tau_m = (3 + 1/(exp((v + 27)/10) + exp(-(v + 102)/15)))/Phi_m;
  h_inf = 1/(1 + exp((v + 80)/5));
  tau_h = (85 + 1/(exp((v + 48)/4) + exp(-(v + 407)/50)))/Phi_h;
  fm = -(m - m_inf)/tau_m;                                  
  fh = -(h - h_inf)/tau_h;                                  
}

//--------------fast Na and K current (RE and TC cells)------------------
double INaK::E_K = -95, INaK::E_Na = 50, INaK::Cels = 36; 

void INaK::calc(double m, double h, double n, double &fm, double &fh, double &fn,
                double v, double x){
  v2 = v - Vtr;
  v2K = v - VtrK;
  iNa = G_Na*m*m*m*h*(v - E_Na);
  Alpha1 = 0.32*(13 - v2)/(exp((13 - v2)/4) - 1);
  Beta1 = 0.28*(v2 - 40)/(exp((v2 - 40)/5) - 1);
  tau_m = 1/(Alpha1 + Beta1) / Phi;
  m_inf = Alpha1/(Alpha1 + Beta1);

  Alpha2 = 0.128*exp((17 - v2)/18);
  Beta2 = 4/(exp((40 - v2)/5) + 1);
  tau_h = 1/(Alpha2 + Beta2) / Phi;
  h_inf = Alpha2/(Alpha2 + Beta2);

  fm = -(m - m_inf)/tau_m;                 
  fh = -(h - h_inf)/tau_h;                 

  iK = G_K* n*n*n*n*(v - E_K);    
  Alpha3 = 0.032*(15 - v2K)/(exp((15 - v2K)/5) - 1);
  Beta3 = 0.5*exp((10 - v2K)/40);
  tau_n = 1/(Alpha3 + Beta3) / Phi;
  n_inf = Alpha3/(Alpha3 + Beta3);

  fn  = -(n - n_inf)/tau_n;                 
}

//------------------Ca-dynamics------------------------------------
double ICa::Ca_inf = 2.4e-4;
double ICa::K_T = 0.0001, ICa::K_d = 0.0001;

void ICa::calc(double cai, double &fcai, double iT, double x) {
  drive = -drive0 * iT / D;
  if(drive < 0) drive = 0;
  fcai = drive + (Ca_inf - cai)/Taur; // - K_T*cai/(cai + K_d);
}

//------------------Low-theshold Ca2+ current (TC cell)-----------------
double IT_TC::Shift = 2, IT_TC::Ca_0 = 2, IT_TC::Cels = 36;
double IT_TC::Qm = 3.55, IT_TC::Qh = 3; //2.8;

void IT_TC::calc(double m, double h, double &fm, double &fh,
                 double v, double cai, double x) {
  ratio = Ca_0/cai;
  if(ratio <= 0.)printf("\n LOG ERROR: TC: cai=%lf ratio=%lf",cai,ratio);
  eca = eca0 * log(ratio);
  iT = G_Ca*m*m*h*(v - eca); 

  m_inf = 1 / (1+exp(-(v+59)/6.2));//////////////////////////////////////
  h_inf = 1 / (1+exp((v+83)/4.)); /////////////////////////////////////

  tau_m = (1/(exp(-(v+131.6)/16.7)+exp((v+16.8)/18.2)) + 0.612) / Phi_m;
  tau_h = (30.8 + (211.4 + exp((v + Shift + 113.2)/5))/
           (1+exp((v + Shift + 84)/3.2))) / Phi_h;

  fm = -(1/tau_m)*(m - m_inf);                                
  fh = -(1/tau_h)*(h - h_inf);
}

//----------------- h-current (TC cell) -----------------------------------
double Ih_TC::E_h = -40; 
double Ih_TC::Shift = 0, Ih_TC::Cels = 36;
double Ih_TC::k2 = 0.0004;
double Ih_TC::nca = 4, Ih_TC::nexp = 1, Ih_TC::taum = 20;

void Ih_TC::calc(double o1, double p1, double o2, double &fo1, double &fp1, 
                 double &fo2, double v, double cai, double x) {
  ih = G_h*(o1 + ginc * o2)*(v - E_h);
  h_inf = 1/(1 + exp((v + 75 - Shift + fac_gh_TC)/5.5));
  tau_s = (taum + 1000 / (exp((v + 71.5 - Shift + fac_gh_TC)/14.2) + 
                          exp(-(v + 89 - Shift + fac_gh_TC)/11.6))) / Phi;
  alpha = h_inf/tau_s;
  beta = (1 - h_inf)/tau_s;
  k1ca = k2 * pow((cai/cac),nca);
  k3p = k4 * pow((p1/pc),nexp);
  fo1 = alpha * (1-o1-o2) - beta * o1; // + k4 * o2 - k3p * o1;
  fp1 = k1ca * (1-p1) - k2 * p1; // + k4 * o2 - k3p * o1;
  fo2 = k3p * o1 - k4 * o2;
}

//----------------------Potassium A-current (TC cell)------------------------
double IA_TC::Cels = 36, IA_TC::E_K = -95;

void IA_TC::calc(double m, double h, double &fm, double &fh, double v, double x){
  iA = G_A*m*m*m*m*h*(v - E_K);
  tau_m = (1.0/( exp((v+35.82)/19.69)+exp(-(v+79.69)/12.7) ) +0.37) / Tad;
  m_inf = 1.0 / (1+exp(-(v+60)/8.5));
  tau_h = 1.0/((exp((v+46.05)/5)+exp(-(v+238.4)/37.45))) / Tad;
  if(v >= -63) tau_h = 19.0/Tad;
  h_inf = 1.0/(1+exp((v+78)/6));
  fm = -(1/tau_m)*(m - m_inf);       
  fh = -(1/tau_h)*(h - h_inf);
}

//---------------------Hight-threshold Ca2+ current (CX cell)----------------
double IHVA_CX::Shift = 0, IHVA_CX::Ca_0 = 2, IHVA_CX::E_Ca = 140; 
double IHVA_CX::Qm = 2.3, IHVA_CX::Qh = 2.3, IHVA_CX::Cels = 36;

void IHVA_CX::calc(double m, double h, double &fm, double &fh,
                   double v, double cai, double x) {
  //------------ECa is fixed (=140mV) instead of using Nerst eq.-----------------

  iHVA = Phi_m * G_HVA * m*m*h * (v - E_Ca);
  vm = v + Shift;

  a = 0.055*(-27 - vm)/(exp((-27-vm)/3.8) - 1);
  b = 0.94*exp((-75-vm)/17);
  tau_m = (1/(a+b))/Phi_m;
  m_inf = a/(a+b);

  a = 0.000457*exp((-13-vm)/50);
  b = 0.0065/(exp((-vm-15)/28) + 1);
  tau_h = (1/(a+b))/Phi_h;
  h_inf = a/(a+b);

  fm = -(m - m_inf)/tau_m;                                  
  fh = -(h - h_inf)/tau_h;
}

//--------------Ca-dependent potassium current (CX cell)-----------------------
double IKCa_CX::E_KCa = -90, IKCa_CX::Q = 2.3, IKCa_CX::caix = 1;
double IKCa_CX::Ra = 0.01, IKCa_CX::Rb = 0.02, IKCa_CX::Cels = 36;   

void IKCa_CX::calc(double m, double &fm, double v, double cai, double x){
  iKCa = Tad * G_KCa * m * (v - E_KCa);                         

  //a = Ra * pow(cai,caix);
  a = Ra * cai;  
  b = Rb;
  tau_m = (1/(a+b))/Tad;
  m_inf = a/(a+b);
  fm = -(1/tau_m)*(m - m_inf);                   
}

//--------------------Potassium M-current (CX cell)-------------------------
double IKm_CX::E_Km = -90, IKm_CX::Q = 2.3;
double IKm_CX::tha = -30, IKm_CX::qa = 9;
double IKm_CX::Ra = 0.001, IKm_CX::Rb = 0.001, IKm_CX::Cels = 36;   

void IKm_CX::calc(double m, double &fm, double v, double x){
  iKm = Tad * G_Km * m * (v - E_Km);
  a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa));
  b = -Rb * (v - tha) / (1 - exp((v - tha)/qa));
  tau_m = (1/(a+b))/Tad;
  m_inf = a/(a+b);
  fm = -(1/tau_m)*(m - m_inf);                   
}

/* MAX MODIF */
// NEED TO ADD LEAK POTASSIUM CURRENT HERE
// FOR NEUROMODULATORS PURPOSES
// --> Preliminary changes with increase in M-current conductance
/* MAX MODIF END */

//--------------------Fast potassium current (CX cell)-------------------
double IKv_CX::E_Kv = -90, IKv_CX::Q = 2.3;
double IKv_CX::tha = 25 /*25*/, IKv_CX::qa = 9;
double IKv_CX::Ra = 0.02, IKv_CX::Rb = 0.002, IKv_CX::Cels = 36;   

void IKv_CX::calc(double m, double &fm, double v, double x){
  g_Kv = Tad * G_Kv * m;
  iKv = g_Kv * (v - E_Kv);                         
  a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa));
  b = -Rb * (v - tha) / (1 - exp((v - tha)/qa));
  tau_m = (1/(a+b))/Tad;
  m_inf = a/(a+b);
  fm = -(1/tau_m)*(m - m_inf);                   
}

//----------------Fast sodium current (CX cell)--------------------------
double INa_CX::Shift = -10, INa_CX::E_Na = 50; 
double INa_CX::Qm = 2.3, INa_CX::Qh = 2.3, INa_CX::Cels = 36;
double INa_CX::tha = -35, INa_CX::qa = 9;
double INa_CX::Ra = 0.182,INa_CX::Rb = 0.124; 
double INa_CX::thi1 = -50, INa_CX::thi2 = -75, INa_CX::qi = 5;
double INa_CX::thinf = -65, INa_CX::qinf = 6.2;
double INa_CX::Rg = 0.0091, INa_CX::Rd = 0.024; 

void INa_CX::calc(double m, double h, double &fm, double &fh,
                  double v, double x) {

  g_Na = Phi_m * G_Na * m*m*m*h;
  iNa = g_Na * (v - E_Na);
  vm = v + Shift;

  a = trap0(vm,tha,Ra,qa);
  b = trap0(-vm,-tha,Rb,qa);
  tau_m = (1/(a+b))/Phi_m;
  m_inf = a/(a+b);

  //"h" inactivation 
  a = trap0(vm,thi1,Rd,qi);
  b = trap0(-vm,-thi2,Rg,qi);
  tau_h = (1/(a+b))/Phi_h;
  h_inf = 1/(1+exp((vm-thinf)/qinf));

  fm = -(m - m_inf)/tau_m;                                  
  fh = -(h - h_inf)/tau_h;
}

//-------------------Persist Sodium current (CX cell)-------------------
void INap_CX::calc(double m, double &fm, double v, double x){
  g_Nap = G_Nap * m;
  iNap = g_Nap * (v - INa_CX::E_Na);
  m_inf = f/(1 + exp(-(v - Tet)/Sig));
  fm = -(m - m_inf)/tau_m; 
}
double INap_CX::cond(){
  return(m_inf);
}

// ----------------------------------------------------------------------
//---------SYNAPCES DESCRIPTION-------------------------------------------
//---------first order kinet model for GABA-A synapse---------------------
double GABA_A::Cdur = 0.3, GABA_A::Cmax = 0.5, GABA_A::Deadtime = 1;
double GABA_A::Prethresh = 0;
double GABA_A::calc(double x, double y_pre, int y_post) { 

	q = ((x - lastrelease) - Cdur);         
	if (q > Deadtime) {
	  if (y_post) {        
      C = Cmax;                
      R0 = R;
      lastrelease = x;
	  }
	}else if (q < 0) { 
	}else if (C == Cmax) {                  
	  R1 = R;
	  C = 0.;
	}
	if (C > 0) {                            
	  R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
	}else{                              
	  R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
	}

  I = strength * fac_GABA_TC * R * (y_pre - E_GABA);
  return I;
}

//------second order kinet model (including G-proteins) for GABA-B synapse----
double GABA_B::E_GABA = -95, GABA_B::Cmax = 0.5, GABA_B::Deadtime = 1;
double GABA_B::Prethresh = 0;
double GABA_B::Kd = 100, GABA_B::n = 4; 

//peter
//y_pre is itself
//y_post is who it is getting it from
//this seems backwards
//if there is a message we reverse weather it is spiking or not
void GABA_B::calc_gaba_b(double r, double g, double &fr, double &fg, 
                         double g_GABA_B, double x, double y_pre, int y_post) {


  //if(y_post){
  //if(old_v == 1){
  //  old_v = 0;
  //}else{
  //  old_v = 1;
  //}
  //}

  Gn = pow(g,n); 
  Gn1 = Gn/(Gn + Kd);
  I = fac_GABA_TC * g_GABA_B * Gn1 * (y_pre - E_GABA);

  q = ((x - lastrelease) - Cdur);         
  if (q > Deadtime) {
    if (y_post) {        
      C = Cmax;                
      lastrelease = x; 
    }
  }else if (q < 0) {                     
  }else if (C == Cmax){
    C = 0;  
  }
  fr = K1 * C * (1 - r) - r * K2;
  fg = K3 * r - K4 * g;
}

double GB::calc_gb(double x, double *y, double *f, double y_pre, int y_post){
  GABA_B::calc_gaba_b(y[0], y[1], f[0], f[1], strength, x, y_pre, y_post); 
  return I;
}; 

//------------first order kiner model for AMPA synapse---------------------
double AMPA::E_AMPA = 0, AMPA::Cdur = 0.3, AMPA::Cmax = 0.5, AMPA::Deadtime = 1;
double AMPA::Cdel = 0; //synaptic delay to get some latency between depolarization
//of the presynaptic membrain and transmitter release
double AMPA::Prethresh = 0, AMPA::Alpha = 1.1, AMPA::Beta = 0.19;
double AMPA::de = 0.9981;//exp(-Beta * (HHH/2));
double AMPA::Cconst=(exp(Beta*Cdur)-exp(-Alpha*Cmax*Cdur))*(Alpha*Cmax/(Alpha*Cmax+Beta)); // R(spike_time)-R_inf; Will produce higher values for t<Cdur
double AMPA::calc(double x, double y_pre, int y_post) {

  q = ((x - lastrelease) - Cdur); 
  
  if(q > Deadtime) {
    if(y_post) {
      if( (x - lastspike) > (Cdel + Cdur) ){
        lastspike = x;
        s = 1; 
      } 
    }  //the flag that spike was but wasn't utilized yet

    if((s == 1) && ((x - lastspike) > Cdel)) {
      s = 0; //spike was utilized
      C = Cmax;                
      R0 = R;
      lastrelease = x;  
    }
  }else if (q < 0){     // x-lr<0.3
  }else if (C == Cmax){ // q>0 && q<1              
    R1 = R;
    C = 0.;
  }

  if (C > 0) {                            
    R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
  } else {                              
    R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
  }

  // I = strength * R * (y_pre - E_AMPA);

  I = strength * fac_AMPA_TC * R * (y_pre - E_AMPA);
  return I;
}

// void AMPA::reset_R(double x) {
//   R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
// }


// double* generate_decay(double *decay,int size, double initial,double de){
//   decay[0]=initial;
//   for(int i=1; i<size; i++){
// 	decay[i]=decay[i-1]*de;
//   }
//   return decay;
// }


//------------first order kiner model for AMPA synapse WITH depression------
double AMPA_D1::E_AMPA = 0, AMPA_D1::Cdur = 0.3, AMPA_D1::Cmax = 0.5;
double AMPA_D1::Deadtime = 1;
double AMPA_D1::Cdel = 0; 
double AMPA_D1::Prethresh = 0, AMPA_D1::Alpha = 1.1, AMPA_D1::Beta = 0.19;

double AMPA_D1::calc(double x, double y_pre, int y_post) {
  q = ((x - lastrelease) - Cdur); 
  
  if(q > Deadtime) {
    if(y_post) {
      if( (x - lastspike) > (Cdel + Cdur) ){
        lastspike = x;
        s = 1; 
      } 
    }  //the flag that spike was but wasn't utilized yet

    if((s == 1) && ((x - lastspike) > Cdel)) {
      s = 0; //spike was utilized
      C = Cmax;                
      R0 = R;
      lastrelease = x;  
    }
  }else if (q < 0){     // x-lr<0.3
  }else if (C == Cmax){ // q>0 && q<1              
    R1 = R;
    C = 0.;
  }

  if (C > 0) {                            
    R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
  } else {                              
    R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
  }

  // I = strength * R * (y_pre - E_AMPA);

  I = strength * fac_AMPA_D2 * R * (y_pre - E_AMPA);
  return I;

}

//--first order kiner model for AMPA synapse WITH depression & spont releases------
double AMPA_D2::E_AMPA = 0, AMPA_D2::Cdur = 0.3;
double AMPA_D2::Cmax = 0.5, AMPA_D2::Deadtime = 1;
double AMPA_D2::Cdel = 0;
double AMPA_D2::Prethresh = 0, AMPA_D2::Alpha = 1.1, AMPA_D2::Beta = 0.19;

//0.30303,0.151515,time,own volts,input volts~-50 to pos 50  mostly -50s
//strength,mini amplitude,current time,my_voltage,pre_spike?,mini_frequency
double AMPA_D2::calc(double x, double y_pre, int y_post){

	q = ((x - lastrelease) - Cdur);     //lastrelease is time of last event either spike or mini
	q1 = ((x - lastrelease1) - Cdur);   //lastrelease1 is time of last spike; it is used to calculate depression length
	if(q > Deadtime){ //if it has been long enough seince last spike or mini
	  if(y_post){ //if the pre cell spiked
      //Use = 0.073; //0.08;
      C = Cmax;                
      R0 = R;
      E = 1 - (1 - E*(1-Use)) * exptable(-q1/Tr);
      lastrelease = x; // this equas when we spiked
      lastrelease1 = x;
      g1 = strength*E;
      if(strength <= 0){
        printf("ampa d2 neg or 0 strength \n");
        exit(1);
      }
	  }else if( ((x - lastrelease1) > 100.0) && ((x - lastrelease) > newrelease) ) {//if it has been long enough since last spike and mini

      //set when the next mini should go off based on the time since the last spike
      //SS = (2.0/(1.0+exp(-(x-lastrelease1)/mini_fre))-1)/100.0;      
      //S = rand()/(RAND_MAX + 1.0); // other way of doing rand
      //SS = log((x - lastrelease1 +Tau)/Tau)/mini_fre;
      SS = (2.0/(1.0+exp(-(x-lastrelease1)/mini_fre))-1.0)/250.0;
      drand48_r(rand_buffer,&S); //uniform random number between 1 and 0

      if(S < 0.000001) S = 0.000001;
      newrelease = -(log(S))/SS;

      g1 = mini_s;
      //Use = 0; //0.01;
      C = Cmax;                
      R0 = R;
      //E = 1 - (1 - E*(1-Use)) * exptable(-q1/Tr);
      lastrelease = x; //This means we released transmitter
	  }
	}else if (q < 0){                     
	}else if (C == Cmax) {                  
	  R1 = R;
	  C = 0.;
	}

	if (C > 0) {                            
	  R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
	} else {                              
	  R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
	}

  I = g1 * fac_AMPA_D2 * R * (y_pre - E_AMPA);

  return I;
}


//------------------D3--------------------- 
//ampa connection with built in stdp

double AMPA_D3::E_AMPA = 0, AMPA_D3::Cdur = 0.3;
double AMPA_D3::Cmax = 0.5, AMPA_D3::Deadtime = 1;
double AMPA_D3::Cdel = 0; //synaptic delay to get some latency between depolarization
//of the presynaptic membrain and transmitter reliase
double AMPA_D3::Prethresh = 0, AMPA_D3::Prethresh1 = -30, AMPA_D3::Alpha = /*0.94*/1.1, AMPA_D3::Beta = 0.19;/*0.18*/ //changed to make it like D2//
//--------------------------------
double AMPA_D3::calc(double time, double y_post, int y_pre) {


  //g_ampa0 is what is actually used
  //g_ampa is the starting value
  // g_ampa00 is probably useless

  if(time<100){
    g_AMPA0=strength; //for first 100 milliseconds use starting strength
    if(strength == 0){
      MinFactor=1;
    }else{
      MinFactor=mini_s/strength;
    }
  }

  q = ((time - lastrelease) - Cdur);
  q1 = ((time - lastrelease1) - Cdur);
  
  //plasticity part
  if(PS) { // if plasticity is on 
    if( ((time-spike_post) > 1.0) && (y_post > Prethresh1) ){ //legit post spike has happend
      if(keyP==1) { //some pre spike has also happend
        //majic number is a mean from some undefined input file we use
        g_AMPA00 = g_AMPA0 + cP*0.1 *exp(-(time-spike_pre)/TauL); //strenthen based on time since pre spike //0.1 //0.8
        // printf("\n weight change UP: t= %lf gOLD= %lf gNEW= %lf", x, g_AMPA0, g_AMPA00);
        if(g_AMPA00 < 3*strength) { //only works if less than 3 times starting strength
          g_AMPA0 = g_AMPA00;
        }else{
          g_AMPA0 = 3*strength;
          printf("weight violation1: time= %lf g_ampa0= %lf  g_ampa00=%lf cp=%lf g_ampa=%lf spike_pre =%lf Taul=%lf \n",time, g_AMPA0,g_AMPA00,cP, strength,spike_pre,TauL);
        }
        keyP=0; // can not have another + stdp even until another prespike
      }
      spike_post=time;
      keyD=1; //a post spike has happend marker
    }

    if( ((time-spike_pre) > 0) && (y_pre > Prethresh) ){ //legit pre-spike has happend
      if(keyD==1) { // a post spike has happened before
        //g_AMPA00 = g_AMPA0 - cD*strength*exp(-(time-spike_post)/TauL); //weaken based on time since post spike
        g_AMPA00 = g_AMPA0 - cD*g_AMPA0*exp(-(time-spike_post)/TauL); //weaken based on time since post spike
        //printf("\n weight change DN: t= %lf gOLD= %lf gNEW= %lf", x, g_AMPA0, g_AMPA00);
        if(g_AMPA00 > 0){
          g_AMPA0 = g_AMPA00;
        }else{
          g_AMPA0 = 0;
          printf("weight violation2: time=%lf g_ampa0=%lf cD=%lf spike_post=%lf TauL=%lf f_ampa00=%lf g_ampa0= %lf  \n", time, g_AMPA0,cD,spike_post,TauL,g_AMPA00,g_AMPA0);
        }
        keyD=0; 
      }
      spike_pre=time;
      keyP=1;
    }

  }

  //spiking part
  if(q > Deadtime) { // if it has been long enough since last spike
    if(y_pre > Prethresh) { //if we get enough input to spike
      //spike
      g1 = E*g_AMPA0;
      factor = 1;
      Use = /*0.1*/0.07; //0.077; //0.08;//changed to make it like D2//
      C = Cmax;                
      R0 = R;
      E = 1 - (1 - E*(1-Use)) * exptable(-q1/Tr);
      lastrelease = time;
      lastrelease1 = time;
    }else if( ((time - lastrelease1) > 70.0) && ((time - lastrelease) > newrelease) ) { //mini spike
      //SS = 1/(1+exp(-(x - lastrelease1 - Period)/Tau)); 

      SS = log((time - lastrelease1 +Tau)/Tau)/mini_fre;
      drand48_r(rand_buffer,&S);
      //S = rand()/(RAND_MAX + 1.0);y 0.30 
      if(S < 0.000001) S = 0.000001;
      newrelease = -(log(S))/SS; //decide when 
      g1 = g_AMPA0*MinFactor; //g_AMPAmin;
      factor = 2;
      //Use = 0; //0.01;
      C = Cmax;                
      R0 = R;
      //E = 1 - (1 - E*(1-Use)) * exptable(-q1/Tr);
      lastrelease = time;
    }
  }else if (q < 0) {                     
  }else if (C == Cmax) {                  
    R1 = R;
    C = 0.0;
  }

  //figure out decay
  if (C > 0) {                            
    R = Rinf + (R0 - Rinf) * exptable (-(time - lastrelease) / Rtau);
  } else {                              
    R = R1 * exptable (-Beta * (time - (lastrelease + Cdur)));
  }

  I = g1 * R * (y_post - E_AMPA);

  return I;
  //       I = (g_AMPA/factor) * R * (y_post - E_AMPA);
}


//------------first order kiner model for NMDA synapse WITH depression------
double NMDA_D1::E_NMDA = 0, NMDA_D1::Cdur = 0.3, NMDA_D1::Cmax = 0.5;
double NMDA_D1::Deadtime = 1;
double NMDA_D1::Cdel = 0; 
double NMDA_D1::Prethresh = 0, NMDA_D1::Alpha = 1, NMDA_D1::Beta = 0.0067;

double NMDA_D1::calc(double x, double y_pre, int y_post) {

	q = ((x - lastrelease) - Cdur); 

	if(q > Deadtime) {
	  if(y_post) {
      if( (x - lastspike) > (Cdel + Cdur) ){
        lastspike = x;
        s = 1; //the flag that spike was but wasn't utilized
      } 
	  }  

	  if((s == 1) && ((x - lastspike) > Cdel)) {
      s = 0; //spike was utilized
      C = Cmax;                
      R0 = R;
      E = 1 - (1 - E*(1-Use)) * exptable(-q/Tr);
      lastrelease = x;  
	  }
	} else if (q < 0) {                     
	} else if (C == Cmax) {                  
	  R1 = R;
	  C = 0.;
	}

	if (C > 0) {                            
	  R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
	} else {                              
	  R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
	}

  fn = 1/(1+exp(-(y_pre - (-25))/12.5));
  I = strength * R * fn * E * (y_pre - E_NMDA);
  return I;
}

//---first order kiner model for GABA-A synapse with DEPRESSION & spont IPSPs--
double GABA_A_D2::Cdur = 0.3, GABA_A_D2::Cmax = 0.5, GABA_A_D2::Deadtime = 1;
double GABA_A_D2::Prethresh = 0, GABA_A_D2::Alpha = 0.53, GABA_A_D2::Beta =0.18; 
double GABA_A_D2::calc(double x, double y_pre, int y_post) {

	q = ((x - lastrelease) - Cdur);
	q1 = ((x - lastrelease1) - Cdur);
	if (q > Deadtime) {
	  if (y_post) {        
      //Use = 0.07;
      C = Cmax;                
      R0 = R;
      E = 1 - (1 - E*(1-Use)) * exptable(-q1/Tr);
      lastrelease = x;
      lastrelease1 = x;
      g1 = strength*E;
	  }else if( ((x - lastrelease1) > 70.0) && ((x - lastrelease) > newrelease) ){
      // SS = log((x - lastrelease1 +Tau)/Tau)/400;
      // //S = rand()/(RAND_MAX + 1.0);
      // drand48_r(rand_buffer,&S);
      // if(S < 0.000001) S = 0.000001;
      // newrelease = -(log(S))/SS;

      SS = (2.0/(1.0+exp(-(x-lastrelease1)/mini_fre))-1.0)/250.0;
      drand48_r(rand_buffer,&S); //uniform random number between 1 and 0

      if(S < 0.000001) S = 0.000001;
      newrelease = -(log(S))/SS;

      //Use = 0; //0.01;
      C = Cmax;                
      R0 = R;
      //E = 1 - (1 - E*(1-Use)) * exptable(-q1/Tr);
      lastrelease = x;
      g1 = mini_s;

	  }
	} else if (q < 0) { } 
	else if (C == Cmax) {                  
	  R1 = R;
	  C = 0.;
	}
	if (C > 0) {                            
	  R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
	} else {                              
	  R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
	}

  I = g1 * fac_GABA_D2 * R * (y_pre - E_GABA);

	// int st=10000; 
	// int et=15000;
	// double fac=2.0;

	// if (x>=et) {       
	//   I = g1 * fac * R * (y_pre - E_GABA);	}  
	// else{ 
	//   if ((x>st) & (x<et)) {       
	//    	I = g1 * (1 + (((et-x)/(et-st)) * (fac-1))) * R * (y_pre - E_GABA);
	//    }
	//    else{
	// 	I = g1 * R * (y_pre - E_GABA);
	//    }
	// }

  return I;

}


// ----------------------------------------------------------------------
//---------MAP SYNAPCES-------------------------------------------

double AMPAmapD1::E_AMPA = 0; //, AMPAmapD1::gamma = 0.98; //0.6;
double AMPAmapD1::d_dep=0.2, AMPAmapD1::d_rec=0.00005;
// double AMPAmapD1::calc(double g_AMPA, double g_AMPAmin, double x, double x_post, double x_pre,double mini_fre) {
//   if(x>Tcr){
//     if(x_pre > 0.1){
//       g = gamma * g + d*g_AMPA/15;   //factor is to match EPSC amplitude to the kinetics model
//       d = (1.0 - d_dep)*d;
//       lastrelease = x;   //any release
//       lastrelease1 = x;  //only spike  
//     }else if( ((x - lastrelease1) > 70.0) && ((x - lastrelease) > newrelease) ) {
//       SS = log((x - lastrelease1 +Tau)/Tau)/mini_fre;
//       drand48_r(rand_buffer,&S);
//       //S = rand()/(RAND_MAX + 1.0);
//       if(S < 0.000001) S = 0.000001;
//       newrelease = -(log(S))/SS;
//       g = gamma * g + g_AMPAmin/5;      //factor is to match EPSC amplitude to the kinetics model
//       lastrelease = x;
//     }else{
//       g = gamma * g;
//       d = 1.0 - (1.0 - d_rec) * (1.0 - d);
//     }
//     I = g * (x_post - E_AMPA);
//     return I;
//     Tcr=x+h;
//   }
//   printf("this functions should not be used in its current state\n");
//   exit(1);
//   return I;
// }


/// From Emily
// double AMPAmapD1::calc(double g_AMPA, double g_AMPAmin, double x, double x_post,  double x_pre) {
double AMPAmapD1::calc(double x, double y_pre, int y_post){
  // calc(double g_AMPA, double g_AMPAmin, double x, double x_post, double x_pre,double mini_fre) {
  // g_AMPA=strength;
	if(y_post > 0.1){
	  g = gamma * g + d*strength/15;   //factor is to match EPSC amplitude to the kinetics model
	  d = (1.0 - d_dep)*d;
	  lastrelease = x;   //any release
	  lastrelease1 = x;  //only spike  
	}
	else if( ((x - lastrelease1) > 70.0) && ((x - lastrelease) > newrelease) ) {   //70.0 is replaced by 70000000.0


	  // S = rand()/(RAND_MAX + 1.0);
    drand48_r(rand_buffer,&S);
	  if(S < 0.000001) S = 0.000001;

    SS = log((x - lastrelease1 +Tau)/Tau)/mini_fre;
	  newrelease = -(log(S))/SS;


	  // g = gamma * g + g_AMPAmin/5;      //factor is to match EPSC amplitude to the kinetics model
	  /// **** REQUIRED TO BE FIXED *****
	  
	  /* //----- In previous code 
       g = gamma * g + strength/5;      //factor is to match EPSC amplitude to the kinetics model
	  */

	  g = gamma * g + mini_s/5; 

	  // (strength/15)/100; // Make mini's some constant * strength as AMPA strength
	  // this scaling is based on testing --- need to be made as a function of mini strength -- OR to think about the correct way to do

	  lastrelease = x;
	}
	else{
	  g = gamma * g;
	  d = 1.0 - (1.0 - d_rec) * (1.0 - d);
	}

	I = g * (y_pre - E_AMPA);
	return I;
	//Tcr=x+h;

}


//=------------------------------
double AMPAmapD::E_AMPA = 0; //, AMPAmapD::gamma = 0.98; //0.6;
double AMPAmapD::d_dep=0.1, AMPAmapD::d_rec=0.00005;
// double AMPAmapD::calc(double g_AMPA, double x, double x_post, double x_pre) {
//   if(x>Tcr){
//     if(x_pre > 0.1){
//       g = gamma * g + d*g_AMPA/15;   //factor is to match EPSC amplitude to the kinetics model
//       d = (1.0 - d_dep)*d;  
//     }else{
//       g = gamma * g;
//       d = 1.0 - (1.0 - d_rec) * (1.0 - d);
//     }
//     I = g * (x_post - E_AMPA);
//     return I;
//     Tcr=x+h;
//   }
//   printf("this functions should not be used in its current state\n");
//   exit(1);
//   return I;
// }

/// From Emily
// double AMPAmapD::calc(double g_AMPA, double x_post, double x_pre) {
double AMPAmapD::calc(double x, double y_pre, int y_post){
  // (double g_AMPA, double x, double x_post, double x_pre) {
  // g_AMPA=strength;
  if(y_post){
    I = gamma * I - d*strength * (y_pre - E_AMPA);
    d = (1.0 - d_dep)*d;
  }
  else{
    I = gamma * I;
    d = 1.0 - (1.0 - d_rec) * (1.0 - d);
  }
  return I;
}


//=------------------------------
double AMPAmap::E_AMPA = 0; //, AMPAmap::gamma = 0.98; //0.6;
double AMPAmap::calc(double x, double y_pre, int y_post){
  // calc(double g_AMPA, double x, double x_post, double x_pre) {
  if(x>Tcr){
    if(y_post){
      g = gamma * g + strength/15;   //factor is to match EPSC amplitude to the kinetics model  
    }else{
      g = gamma * g;
    }
    I = g * (y_pre - E_AMPA);
    return I;
    Tcr=x+h;
  }
  printf("this functions should not be used in its current state\n");
  exit(1);
  return I;
}


//=------------------------------
double NMDAmapD1::E_NMDA = 0; //, NMDAmapD1::gamma = 0.99925; //0.6;
double NMDAmapD1::d_dep=0.0, NMDAmapD1::d_rec=0.00005;
// double NMDAmapD1::calc(double g_NMDA, double x, double x_post, double x_pre) {
//   if(x>Tcr){
//     if(x_pre > 0.1){
//       g = gamma * g + d*g_NMDA/15;   //factor is to match EPSC amplitude to the kinetics model
//       d = (1.0 - d_dep)*d;
//     }else{
//       g = gamma * g;
//       d = 1.0 - (1.0 - d_rec) * (1.0 - d);
//     }
//     fn = 1/(1+exp(-(x_post - (-25))/12.5));
//     I = g * fn * (x_post - E_NMDA);
//     return I;
//     Tcr=x+h;
//   }
//   printf("this functions should not be used in its current state\n");
//   exit(1);
//   return I;
// }

double NMDAmapD1::calc(double x, double y_pre, int y_post){
  // calc(double g_NMDA, double x, double x_post, double x_pre) {
  // g_NMDA=strength;
  if(y_post){
    g = gamma * g + d*strength/15;   //factor is to match EPSC amplitude to the kinetics model
    d = (1.0 - d_dep)*d;
  }

  else{
    g = gamma * g;
    d = 1.0 - (1.0 - d_rec) * (1.0 - d);
  }
  fn = 1/(1+exp(-(y_pre - (-25))/12.5));
  I = g * fn * (y_pre - E_NMDA);
  return I;

}

//=------------------------------
double GABAAmapD1::E_GABAA = -70; //, GABAAmapD1::gamma = 0.9;
//double GABAAmapD1::d_dep=0.5, GABAAmapD1::d_rec=0.00005;
// double GABAAmapD1::calc(double g_GABA, double x, double x_post, double x_pre) {
//   if(x>Tcr){
//     if(x_pre > 0.1){
//       g = gamma * g + g_GABA/15;
//       //d = (1.0 - d_dep)*d;
//     }else{
//       g = gamma * g;
//       //d = 1.0 - (1.0 - d_rec) * (1.0 - d);
//     }
//     I = g * (x_post - E_GABAA);
//     return I;
//     Tcr=x+h;
//   }
//   printf("this functions should not be used in its current state\n");
//   exit(1);
//   return I;
// }

double GABAAmapD1::calc(double x, double y_pre, int y_post) {
  // calc(double g_GABA, double x, double x_post, double x_pre) {
  // No depression is used ???
  // g_GABA=strength;
  if(y_post){
    g = gamma * g + strength/15;
    //           d = (1.0 - d_dep)*d;
  }

  else{
    g = gamma * g;
    //          d = 1.0 - (1.0 - d_rec) * (1.0 - d);
  }
  I = g * (y_pre - E_GABAA);
  return I;
}


//------------------------------------------------------------------------------------------
//-----first order kiner model for AMPA synapse used for external stimulation----
double Extern_ampa::Cdur = 0.3, Extern_ampa::Cmax = 0.5, Extern_ampa::Deadtime = 1;
double Extern_ampa::Prethresh = 0;
void Extern_ampa::calc(double g_Extern_ampa, double x) {

	q = ((x - lastrelease) - Cdur);         
	if (q > Deadtime) {
	  if ((x - lastrelease) > TR) {        
      C = Cmax;                
      R0 = R;
      lastrelease = x;
	  }
	} else if (q < 0) {                     
	} else if (C == Cmax) {                  
	  R1 = R;
	  C = 0.;
	}
	if (C > 0) {                            
	  R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau);
	} else {                              
	  R = R1 * exptable (-Beta * (x - (lastrelease + Cdur)));
	}
	g = g_Extern_ampa * R;

}

//-------------------basic RE CELL-----------------------------------------------
double RE::V0 = -61, RE::Cai0 = 0.0001;

void RE::calc(double x, double I, double *y, double *f){
  IT_RE::calc(y[2], y[3], f[2], f[3], y[0], y[1], x);
  INaK::calc(y[4], y[5], y[6], f[4], f[5], f[6], y[0], x);
  ICa::calc(y[1], f[1], iT, x);

  f[0] = -G_l * (y[0] - E_l) - iT - iNa - iK - fac_gkl_RE * G_kl * (y[0] - INaK::E_K) + I;

  v_DEND = v_SOMA = y[0];

}

//-------------------basic TC CELL-----------------------------------------------
double TC::Cai0 = 0.0001;
double TC::G_l = 0.01;//0.05; //0.01;///
double TC::V0 = -68;

void TC::calc(double x, double I, double *y, double *f){
  IT_TC::calc(y[2], y[3], f[2], f[3], y[0], y[1], x);
  Ih_TC::calc(y[4], y[5], y[6], f[4], f[5], f[6], y[0], y[1], x);
  INaK::calc(y[7], y[8], y[9], f[7], f[8], f[9], y[0], x); 
  IA_TC::calc(y[10], y[11], f[10], f[11], y[0], x);
  ICa::calc(y[1], f[1], iT, x);

  f[0] = -G_l * (y[0] - E_l)  - iNa - iK - ih - 2.0*iT - iA - 1. * fac_gkl_TC * G_kl * (y[0] - INaK::E_K) + I +DC + DC_TC;
  // f[0] = -G_l * (y[0] - E_l) - iT -  ih - iNa - iK - iA -  G_kl * (y[0] - INaK::E_K) + I +DC;
  v_DEND = v_SOMA = y[0];

}

//-------------------basic CX CELL------------------------------------------------
//-------------------CX CELL (DENDRITE)-------------------------------------
void CX_DEND::calc(double x, double *y, double *f){
	IHVA_CX::calc(y[2], y[3], f[2], f[3], y[0], y[1], x);
	IKCa_CX::calc(y[4], f[4], y[0], y[1], x);
	IKm_CX::calc(y[5], f[5], y[0], x);
	INa_CX::calc(y[6], y[7], f[6], f[7], y[0], x);
	ICa::calc(y[1], f[1], iHVA, x);
	INap_CX::calc(y[8], f[8], y[0], x);

	iDEND =  -G_l * (y[0] - E_l) - iHVA -  iKCa - fac_gkm_cx * iKm - iNa -iNap + I_Stim1 - 1.25* fac_gkl * G_kl * (y[0] - INaK::E_K);
}

//-------------------CX CELL (SOMA)-------------------------------------
void CX_SOMA::calc(double x, double *y, double *f){
  IKv_CX::calc(y[0], f[0], v_SOMA, x);
  INa_CX::calc(y[1], y[2], f[1], f[2], v_SOMA, x);
  INap_CX::calc(y[3], f[3], v_SOMA, x);

  g1_SOMA = g_Na + fac_gkv_cx * g_Kv + g_Nap;
  g2_SOMA = g_Na * INa_CX::E_Na + fac_gkv_cx * g_Kv * IKv_CX::E_Kv + g_Nap * INa_CX::E_Na + I_Stim2 + 6.74172; //
  iSOMA =  - iNa - iKv - iNap;     
}

//------------CX CELL (connect DENDRITE and SOMA)---------------------------
double CX::Cai0 = 0.0001, CX::V0 = -68;
double CX::C = 0.75;   // uF/cm^2

void CX::calc(double x, double I, double *y, double *f){
    
  y[0] = y[0] + field_effect_dend;
  v_SOMA = v_SOMA + field_effect_soma;
  //y[N_DEND] = y[N_DEND] + field_effect_soma;
  CX_SOMA::calc(x, y+N_DEND, f+N_DEND);
  v_SOMA = (y[0] + kappa * S_CX_SOMA * g2_SOMA) / (1 + kappa*S_CX_SOMA * g1_SOMA);
  v_DEND = y[0];
  CX_DEND::calc(x, y, f);
  axil_current = 1.0/(kappa*S_CX_DEND) * (v_SOMA - y[0]);
  f[0] = (1.0/C) * ( iDEND + 1.0/(kappa*S_CX_DEND) * (v_SOMA - y[0])) + I;
  
}
//-------------------------------------------
