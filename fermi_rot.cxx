#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"


using namespace std;

//----------------------------------------------------

 void fermi_rot(Float_t &E_beam_fermi, Float_t &theta_rot2,Float_t E_beam,TLorentzVector  P4_E_prime, TLorentzVector  &P4_E_prime_boosted) {
 
//This subroutine performs the Lab-->quasiLab transformation, where quasiLab is the system, where the initial proton is at rest, while the initial electron moves along the Z-axis. The transformation is performed via the sequence of auxiliary systems: Lab-->System 1-->System 2-->quasiLab(System 3).  
//As an output it gives the energy of the incoming electron in the quasiLab (E_beam_fermi), the four-momentum of the scattered electron in the quasiLab (P4_E_prime_boosted), and the rotation angle of the System 2-->quasiLab transformation (theta_rot2).
 
Float_t theta_fermi, phi_fermi;
TLorentzVector P4_EL, P4_in_Prot;

//The four-momenta of the incoming electron and the target proton in the Lab frame. 
P4_EL.SetXYZT(0.,0.,E_beam,E_beam);
P4_in_Prot.SetXYZT(px_fermi,py_fermi,pz_fermi,sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));

//Transformation Lab --> System 1
theta_fermi = acos(pz_fermi/sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
phi_fermi = acos(fabs(px_fermi)/sqrt(px_fermi*px_fermi+py_fermi*py_fermi));

if ((px_fermi < 0.)&&(py_fermi > 0.)) phi_fermi = M_PI-phi_fermi;
if ((px_fermi < 0.)&&(py_fermi < 0.)) phi_fermi = phi_fermi + M_PI;
if ((px_fermi > 0.)&&(py_fermi < 0.)) phi_fermi = 2.*M_PI - phi_fermi;

P4_EL.RotateZ(-phi_fermi);
P4_EL.RotateY(-theta_fermi);
 
P4_E_prime.RotateZ(-phi_fermi);
P4_E_prime.RotateY(-theta_fermi);
 
P4_in_Prot.RotateZ(-phi_fermi);
P4_in_Prot.RotateY(-theta_fermi);


//Transformation System 1 --> System 2 
Float_t beta;
 
beta = sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi)/sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi);

P4_in_Prot.Boost(0,0,-beta);
P4_EL.Boost(0,0,-beta);
P4_E_prime.Boost(0,0,-beta);


//Transformation System 2 --> quasiLab (Syatem 3)
theta_rot2 = acos(P4_EL[2]/sqrt(P4_EL[0]*P4_EL[0]+P4_EL[1]*P4_EL[1]+P4_EL[2]*P4_EL[2]));


P4_EL.RotateY(theta_rot2);
P4_E_prime.RotateY(theta_rot2);
P4_in_Prot.RotateY(theta_rot2);

P4_E_prime_boosted=P4_E_prime;
E_beam_fermi = P4_EL[3]; 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%TEST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//Test whether the target proton is at rest in the quasiLab.  If not --> smth is wrong!
 if ((fabs(P4_in_Prot[0])>0.001)||(fabs(P4_in_Prot[1])>0.001)||(fabs(P4_in_Prot[2])>0.001)||(fabs(P4_in_Prot[3]-MP)>0.0001))  cout << "ALARM! Wrong Lab-->quasiLab transformation! Proton is not at rest in the quasiLab! \n";
 
//Test whether the incoming electron moves along Z-axis in the quasiLab.  If not --> smth is wrong! 
 if ((fabs(P4_EL[0])>0.001)||(fabs(P4_EL[1])>0.001)||(fabs(P4_EL[2]-P4_EL[3])>0.001))  cout << "ALARM! Wrong Lab-->quasiLab transformation! Electron does not move along Z-axis in the quasiLab! \n";

};
 
 
 
void get_EpsL_Ebeam_ferm( Float_t E_beam,Float_t Ep, Float_t phi_e,Float_t theta_e, Float_t &eps_l,Float_t &eps_t, Float_t &E_beam_fermi,Float_t &E_p_fermi) { 

//This subroutine also performs the Lab-->quasiLab transformation, where quasiLab is the system, where the initial proton is at rest, while the initial electron moves along the Z-axis. The transformation is performed via the sequence of auxiliary systems: Lab-->System 1-->System 2-->quasiLab(System 3).  
//As an output it gives the dergees of virtual photon polarization in the quasiLab (eps_l, eps_t) as well as the energies of the incoming and scattered electrons in the quasiLab (E_beam_fermi, E_p_fermi).
//This subroutine is similiar to fermi_rot() and needed for the simulation of the radiative effects.
 
Float_t Q2, theta_fermi, phi_fermi,theta_rot2;
Float_t nu, E_E_prime, Theta_e_prime;
TLorentzVector P4_EL, P4_in_Prot, P4_gamma,P4_E_prime;

//The four-momenta of the incoming electron, target proton, scattered electron, and virtual photon in the Lab frame.  
P4_EL.SetXYZT(0.,0.,E_beam,E_beam);
P4_in_Prot.SetXYZT(px_fermi,py_fermi,pz_fermi,sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
P4_E_prime.SetXYZT(Ep*cos(phi_e)*sin(theta_e),Ep*sin(phi_e)*sin(theta_e),Ep*cos(theta_e),Ep);
P4_gamma = P4_EL -  P4_E_prime;

//Transformation Lab --> System 1
theta_fermi = acos(pz_fermi/sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
phi_fermi = acos(fabs(px_fermi)/sqrt(px_fermi*px_fermi+py_fermi*py_fermi));
 
if ((px_fermi < 0.)&&(py_fermi > 0.)) phi_fermi = M_PI-phi_fermi;
if ((px_fermi < 0.)&&(py_fermi < 0.)) phi_fermi = phi_fermi + M_PI;
if ((px_fermi > 0.)&&(py_fermi < 0.)) phi_fermi = 2.*M_PI - phi_fermi;


P4_EL.RotateZ(-phi_fermi);
P4_EL.RotateY(-theta_fermi);
  
P4_gamma.RotateZ(-phi_fermi);
P4_gamma.RotateY(-theta_fermi);
   
P4_E_prime.RotateZ(-phi_fermi);
P4_E_prime.RotateY(-theta_fermi);
   
P4_in_Prot.RotateZ(-phi_fermi);
P4_in_Prot.RotateY(-theta_fermi);


//Transformation System 1 --> System 2
Float_t beta;
 
beta = sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi)/sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi);

P4_in_Prot.Boost(0,0,-beta);
P4_EL.Boost(0,0,-beta);
P4_E_prime.Boost(0,0,-beta);
P4_gamma.Boost(0,0,-beta);


//Transformation System 2 --> quasiLab (Syatem 3)
theta_rot2 = acos(P4_EL[2]/sqrt(P4_EL[0]*P4_EL[0]+P4_EL[1]*P4_EL[1]+P4_EL[2]*P4_EL[2]));

P4_EL.RotateY(theta_rot2);
P4_gamma.RotateY(theta_rot2);
P4_E_prime.RotateY(theta_rot2);
P4_in_Prot.RotateY(theta_rot2);
  
//----------------------------

E_beam_fermi = P4_EL[3];
E_p_fermi = P4_E_prime[3];
 
Q2 = -P4_gamma.Mag2();

eps_t = 1./(1.+Q2*((P4_gamma.Vect()).Mag2())/2./(((P4_EL.Vect()).Cross(P4_E_prime.Vect())).Mag2()));
eps_l = Q2*eps_t/P4_gamma[3]/P4_gamma[3];


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%TEST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//Test whether the target proton is at rest in the quasiLab.  If not --> smth is wrong!
 if ((fabs(P4_in_Prot[0])>0.001)||(fabs(P4_in_Prot[1])>0.001)||(fabs(P4_in_Prot[2])>0.001)||(fabs(P4_in_Prot[3]-MP)>0.0001))  cout << "ALARM! Wrong Lab-->quasiLab transformation! Proton is not at rest in the quasiLab! \n";

//Test whether the incoming electron moves along Z-axis in the quasiLab.  If not --> smth is wrong! 
 if ((fabs(P4_EL[0])>0.001)||(fabs(P4_EL[1])>0.001)||(fabs(P4_EL[2]-P4_EL[3])>0.001))  cout << "ALARM! Wrong Lab-->quasiLab transformation! Electron does not move along Z-axis in the quasiLab! \n";

};
 
 
void get_rot2( Float_t E_beam,Float_t Ep, Float_t phi_e,Float_t theta_e,Float_t &theta_rot2, Float_t &ph_e_ferm) { 
 
//This subroutine also performs the Lab-->quasiLab transformation, where quasiLab is the system, where the initial proton is at rest, while the initial electron moves along the Z-axis. The transformation is performed via the sequence of auxiliary systems: Lab-->System 1-->System 2-->quasiLab(System 3).  
//As an output it gives the rotation angle of the System 2-->quasiLab transformation (theta_rot2) as well as the azimuthal angle of the scattered electron in quasiLab (ph_e_ferm). 
//This subroutine is similiar to fermi_rot() and needed for the simulation of the radiative effects.
 
Float_t theta_fermi, phi_fermi;
Float_t nu, E_E_prime, Theta_e_prime;
TLorentzVector P4_EL, P4_in_Prot,P4_E_prime;

//The four-momenta of the incoming electron, target proton, and scattered electron in the Lab frame.   
P4_EL.SetXYZT(0.,0.,E_beam,E_beam);
P4_in_Prot.SetXYZT(px_fermi,py_fermi,pz_fermi,sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
P4_E_prime.SetXYZT(Ep*cos(phi_e)*sin(theta_e),Ep*sin(phi_e)*sin(theta_e),Ep*cos(theta_e),Ep);
      
//Transformation Lab --> System 1
theta_fermi = acos(pz_fermi/sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
phi_fermi = acos(fabs(px_fermi)/sqrt(px_fermi*px_fermi+py_fermi*py_fermi));
 
if ((px_fermi < 0.)&&(py_fermi > 0.)) phi_fermi = M_PI-phi_fermi;
if ((px_fermi < 0.)&&(py_fermi < 0.)) phi_fermi = phi_fermi + M_PI;
if ((px_fermi > 0.)&&(py_fermi < 0.)) phi_fermi = 2.*M_PI - phi_fermi;
 

P4_EL.RotateZ(-phi_fermi);
P4_EL.RotateY(-theta_fermi);
  
P4_E_prime.RotateZ(-phi_fermi);
P4_E_prime.RotateY(-theta_fermi);
   
P4_in_Prot.RotateZ(-phi_fermi);
P4_in_Prot.RotateY(-theta_fermi);
  

//Transformation System 1 --> System 2 
Float_t beta;
 
beta = sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi)/sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi);

P4_in_Prot.Boost(0,0,-beta);
P4_EL.Boost(0,0,-beta);
P4_E_prime.Boost(0,0,-beta);


//Transformation System 2 --> quasiLab (Syatem 3)
theta_rot2 = acos(P4_EL[2]/sqrt(P4_EL[0]*P4_EL[0]+P4_EL[1]*P4_EL[1]+P4_EL[2]*P4_EL[2]));

P4_EL.RotateY(theta_rot2);
P4_E_prime.RotateY(theta_rot2);
P4_in_Prot.RotateY(theta_rot2);

ph_e_ferm = P4_E_prime.Phi();

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%TEST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//Test whether the target proton is at rest in the quasiLab.  If not --> smth is wrong!
 if ((fabs(P4_in_Prot[0])>0.001)||(fabs(P4_in_Prot[1])>0.001)||(fabs(P4_in_Prot[2])>0.001)||(fabs(P4_in_Prot[3]-MP)>0.0001))  cout << "ALARM! Wrong Lab-->quasiLab transformation! Proton is not at rest in the quasiLab! \n";

//Test whether the incoming electron moves along Z-axis in the quasiLab.  If not --> smth is wrong! 
 if ((fabs(P4_EL[0])>0.001)||(fabs(P4_EL[1])>0.001)||(fabs(P4_EL[2]-P4_EL[3])>0.001))  cout << "ALARM! Wrong Lab-->quasiLab transformation! Electron does not move along Z-axis in the quasiLab! \n"; 

 }; 
//-------------------------------------
 
 
void from_qulab_to_lab(Float_t theta_rot2,Float_t theta_e_ferm, Float_t ph_e_ferm,Float_t Es,Float_t Ep,Float_t &Es_lab,Float_t &Ep_lab) {

//This subroutine performs the quasiLab-->Lab transformation, where quasiLab is the system, where the initial proton is at rest, while the initial electron moves along the Z-axis. The transformation is performed via the sequence of auxiliary systems: quasiLab(System 3)-->System 2-->System 1-->Lab.
//As an output it gives the energies of the incoming and scattered electrons in the Lab (Es_lab, Ep_lab).
//This subroutine is similiar to fermi_antirot() and needed for the simulation of the radiative effects.

Float_t theta_fermi, phi_fermi;
TLorentzVector P4_Eini_qualab,P4_E_prime_qualab;

//The four-momenta of the incoming electron and scattered electron in the quasiLab frame.  
P4_Eini_qualab.SetXYZT(0.,0.,Es,Es); P4_E_prime_qualab.SetXYZT(Ep*cos(ph_e_ferm)*sin(theta_e_ferm),Ep*sin(ph_e_ferm)*sin(theta_e_ferm),Ep*cos(theta_e_ferm),Ep); 

//Defining of the spatial angles of the Fermi momentum  
theta_fermi = acos(pz_fermi/sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
phi_fermi = acos(fabs(px_fermi)/sqrt(px_fermi*px_fermi+py_fermi*py_fermi));
 
if ((px_fermi < 0.)&&(py_fermi > 0.)) phi_fermi = M_PI-phi_fermi;
if ((px_fermi < 0.)&&(py_fermi < 0.)) phi_fermi = phi_fermi + M_PI;
if ((px_fermi > 0.)&&(py_fermi < 0.)) phi_fermi = 2.*M_PI - phi_fermi;
 
//Defining beta of the boost 
Float_t beta;

beta = sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi)/sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi);
 

//Transformation quasiLab-->System 2 
P4_Eini_qualab.RotateY(-theta_rot2);
P4_E_prime_qualab.RotateY(-theta_rot2);

//Transformation System 2--> System 1  
P4_Eini_qualab.Boost(0,0,beta);
P4_E_prime_qualab.Boost(0,0,beta);

//Transformation System 1--> Lab
P4_Eini_qualab.RotateY(theta_fermi);
P4_Eini_qualab.RotateZ(phi_fermi);

P4_E_prime_qualab.RotateY(theta_fermi);
P4_E_prime_qualab.RotateZ(phi_fermi);
  

Ep_lab = P4_E_prime_qualab[3];
Es_lab = P4_Eini_qualab[3];
   


}; 
 
 
 
