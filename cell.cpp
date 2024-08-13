#include<iostream> //standard header file containing cout, cin etc
#include<algorithm> //for using max function
#include<cstdlib> //for using atof function
#include<cmath> //for using math functions
#include<limits> //for using numeric limits
#include<vector> //for using vectors
#include<fstream> //for using files
#include<sstream> //for using strings
#include<omp.h>
using namespace std;

int N;
double ep, I, Rm, dr, l_B, E_Rc, delta_E_Rc, alpha, chi, theta, a;
vector <double> delta_psi(2, 0.), delta_y(2,0.);

//define function y
double eta(double y)
{
  if (y >= 0.)
    return 1./(1. + exp(-y));
  else
    return exp(y)/(1. + exp(y));
}

//defining Gradient_psi: derivative of omega w. r. to psi_s

double Gradient_psi(int s, const vector<double> &psi, double y, double e_rc)
{
    double sum = 0.0;
    if(s == 0)
    {
        double dpsi = psi[s+1] - psi[s];
        
        double t1_2_c = ((fabs(dpsi)>0.001)? (dpsi-sinh(dpsi))/(dpsi*dpsi) : -dpsi/6.);
        double t1_2_s = ((fabs(dpsi)>0.001)? (1.-cosh(dpsi))/(dpsi*dpsi) : -0.5);
        double t2_2_c = ((fabs(dpsi)>0.01)? (dpsi*sinh(dpsi)-2.*cosh(dpsi)+2.)/(dpsi*dpsi*dpsi) : dpsi/12.);
        double t2_2_s = ((fabs(dpsi)>0.001)? (dpsi*cosh(dpsi)-2.*sinh(dpsi)+dpsi)/(dpsi*dpsi*dpsi) : 1./6.);
        double t3_2_c = ((fabs(dpsi)>0.02)? ((dpsi*dpsi+6.)*sinh(dpsi)-4.*dpsi*cosh(dpsi)-2.*dpsi)/(dpsi*dpsi*dpsi*dpsi) : dpsi/20.);
        double t3_2_s = ((fabs(dpsi)>0.02)? ((dpsi*dpsi+6.)*cosh(dpsi)-4.*dpsi*sinh(dpsi)-6.)/(dpsi*dpsi*dpsi*dpsi) : 1./12.);
        
        double T1_2 = -6.022e-4*8*M_PI*I*dr*(Rm+s*dr)*(Rm+s*dr)*(t1_2_c*cosh(psi[s])+t1_2_s*sinh(psi[s]));
        double T2_2 = 6.022e-4*16*M_PI*I*dr*dr*(Rm+s*dr)*(t2_2_c*cosh(psi[s])+t2_2_s*sinh(psi[s]));
        double T3_2 = 6.022e-4*8*M_PI*I*dr*dr*dr*(t3_2_c*cosh(psi[s])+t3_2_s*sinh(psi[s]));
        double T4_2 = -(ep/l_B)*(Rm*Rm/dr+Rm*(2*s+1)+s*(s+1)+dr/3.)*dpsi;
        double T5 = -4*M_PI*Rm*Rm*(eta(y)-1./theta)/(a*a);
        
        sum+=T1_2+T2_2+T3_2+T4_2+T5;
    }
    else if(s == N)
    {
        double dpsi = psi[s] - psi[s-1];
        
        double t1_1_c = ((fabs(dpsi)>0.001)? (dpsi-sinh(dpsi))/(dpsi*dpsi) : -dpsi/6.);
        double t1_1_s = ((fabs(dpsi)>0.001)? (-1.+cosh(dpsi))/(dpsi*dpsi) : 0.5);
        double t2_1_c = ((fabs(dpsi)>0.01)? (dpsi*dpsi+2.-2.*cosh(dpsi))/(dpsi*dpsi*dpsi) : -dpsi/12.);
        double t2_1_s = ((fabs(dpsi)>0.001)? (-2.*dpsi+2.*sinh(dpsi))/(dpsi*dpsi*dpsi) : 1./3.);
        double t3_1_c = ((fabs(dpsi)>0.02)? (dpsi*dpsi*dpsi+6.*dpsi-6.*sinh(dpsi))/(dpsi*dpsi*dpsi*dpsi) : -dpsi/20.);
        double t3_1_s = ((fabs(dpsi)>0.02)? (-3.*dpsi*dpsi-6.+6.*cosh(dpsi))/(dpsi*dpsi*dpsi*dpsi) : 0.25);
        
        double T1_1 = 6.022e-4*8*M_PI*I*dr*(Rm+(s-1)*dr)*(Rm+(s-1)*dr)*(t1_1_c*cosh(psi[s])+t1_1_s*sinh(psi[s]));
        double T2_1 = 6.022e-4*16*M_PI*I*dr*dr*(Rm+(s-1)*dr)*(t2_1_c*cosh(psi[s])+t2_1_s*sinh(psi[s]));
        double T3_1 = 6.022e-4*8*M_PI*I*dr*dr*dr*(t3_1_c*cosh(psi[s])+t3_1_s*sinh(psi[s]));
        double T4_1 = (ep/l_B)*(Rm*Rm/dr+Rm*(2*s-1)+s*(s-1)+dr/3.)*dpsi;
        double T5 = ep*e_rc;
        
        sum+=T1_1+T2_1+T3_1+T4_1+T5;
    }
    else
    {
        double dpsi = psi[s] - psi[s-1];
        
        double t1_1_c = ((fabs(dpsi)>0.001)? (dpsi-sinh(dpsi))/(dpsi*dpsi) : -dpsi/6.);
        double t1_1_s = ((fabs(dpsi)>0.001)? (-1.+cosh(dpsi))/(dpsi*dpsi) : 0.5);
        double t2_1_c = ((fabs(dpsi)>0.01)? (dpsi*dpsi+2.-2.*cosh(dpsi))/(dpsi*dpsi*dpsi) : -dpsi/12.);
        double t2_1_s = ((fabs(dpsi)>0.001)? (-2.*dpsi+2.*sinh(dpsi))/(dpsi*dpsi*dpsi) : 1./3.);
        double t3_1_c = ((fabs(dpsi)>0.02)? (dpsi*dpsi*dpsi+6.*dpsi-6.*sinh(dpsi))/(dpsi*dpsi*dpsi*dpsi) : -dpsi/20.);
        double t3_1_s = ((fabs(dpsi)>0.02)? (-3.*dpsi*dpsi-6.+6.*cosh(dpsi))/(dpsi*dpsi*dpsi*dpsi) : 0.25);
        
        double T1_1 = 6.022e-4*8*M_PI*I*dr*(Rm+(s-1)*dr)*(Rm+(s-1)*dr)*(t1_1_c*cosh(psi[s])+t1_1_s*sinh(psi[s]));
        double T2_1 = 6.022e-4*16*M_PI*I*dr*dr*(Rm+(s-1)*dr)*(t2_1_c*cosh(psi[s])+t2_1_s*sinh(psi[s]));
        double T3_1 = 6.022e-4*8*M_PI*I*dr*dr*dr*(t3_1_c*cosh(psi[s])+t3_1_s*sinh(psi[s]));
        double T4_1 = (ep/l_B)*(Rm*Rm/dr+Rm*(2*s-1)+s*(s-1)+dr/3.)*dpsi;
        
        dpsi = psi[s+1] - psi[s];
        
        double t1_2_c = ((fabs(dpsi)>0.001)? (dpsi-sinh(dpsi))/(dpsi*dpsi) : -dpsi/6.);
        double t1_2_s = ((fabs(dpsi)>0.001)? (1.-cosh(dpsi))/(dpsi*dpsi) : -0.5);
        double t2_2_c = ((fabs(dpsi)>0.01)? (dpsi*sinh(dpsi)-2.*cosh(dpsi)+2.)/(dpsi*dpsi*dpsi) : dpsi/12.);
        double t2_2_s = ((fabs(dpsi)>0.001)? (dpsi*cosh(dpsi)-2.*sinh(dpsi)+dpsi)/(dpsi*dpsi*dpsi) : 1./6.);
        double t3_2_c = ((fabs(dpsi)>0.02)? ((dpsi*dpsi+6.)*sinh(dpsi)-4.*dpsi*cosh(dpsi)-2.*dpsi)/(dpsi*dpsi*dpsi*dpsi) : dpsi/20.);
        double t3_2_s = ((fabs(dpsi)>0.02)? ((dpsi*dpsi+6.)*cosh(dpsi)-4.*dpsi*sinh(dpsi)-6.)/(dpsi*dpsi*dpsi*dpsi) : 1./12.);
        
        double T1_2 = -6.022e-4*8*M_PI*I*dr*(Rm+s*dr)*(Rm+s*dr)*(t1_2_c*cosh(psi[s])+t1_2_s*sinh(psi[s]));
        double T2_2 = 6.022e-4*16*M_PI*I*dr*dr*(Rm+s*dr)*(t2_2_c*cosh(psi[s])+t2_2_s*sinh(psi[s]));
        double T3_2 = 6.022e-4*8*M_PI*I*dr*dr*dr*(t3_2_c*cosh(psi[s])+t3_2_s*sinh(psi[s]));
        double T4_2 = -(ep/l_B)*(Rm*Rm/dr+Rm*(2*s+1)+s*(s+1)+dr/3.)*dpsi;
        
        sum+=T1_1+T2_1+T3_1+T4_1+T1_2+T2_2+T3_2+T4_2;
    }
    
    return sum;
}

//defining Gradient_y: derivative of omega w. r. to y

double Gradient_y(const vector<double> &psi, double y)
{
    double sum = 0.0;
    
    double eta_y = eta(y);
    
    double T1 = -4*M_PI*Rm*Rm*psi[0]/(a*a);
    double T2 = 4*M_PI*Rm*Rm*(alpha + chi*eta_y - y)/(a*a);
    
    sum+=T1+T2;
    
    return sum*eta_y*(1.-eta_y); //converting derivative w. r. to eta to that w. r. to y
}

//defining Gradient_E_Rc: derivative of omega w. r. to E_Rc

double Gradient_E_Rc(const vector<double> &psi)
{
    return ep*psi[N];
}

int main (int args, char *arg[])
{
    if(args!=22)
    {
        cerr << "Program <ep> <I> <Rm> <dr> <N> <l_B> <tol_psi> <tol_y> <tol_E_Rc> <iter_psi_max> <iter_y_max> <iter_E_Rc_max> <tau_psi> <tau_y> <tau_E_Rc> <alpha> <chi> <theta> <a> <eta_left_init> <eta_rightt_init>" << endl;
        return 1;
    }
    
    ep = atof(arg[1]);
    I = atof(arg[2]);
    Rm = atof(arg[3]);
    dr = atof(arg[4]);
    N = atoi(arg[5]);
    l_B = atof(arg[6]);
    double tol_psi = atof(arg[7]);
    double tol_y = atof(arg[8]);
    double tol_E_Rc = atof(arg[9]);
    int iter_psi_max = atoi(arg[10]);
    int iter_y_max = atoi(arg[11]);
    int iter_E_Rc_max = atoi(arg[12]);
    double tau_psi = atof(arg[13]);
    double tau_y = atof(arg[14]);
    double tau_E_Rc = atof(arg[15]);
    alpha = atof(arg[16]);
    chi = atof(arg[17]);
    theta = atof(arg[18]);
    a = atof(arg[19]);
    double eta_left_init = atof(arg[20]);
    double eta_right_init = atof(arg[21]);

    double E_Rc_0 = 0.;
  
    double E_Rc_new = 0.;
  
    vector<double> y(2,0.), y_new(2,0.);
    y[0] = log(eta_left_init)-log1p(-eta_left_init);
    y[1] = log(eta_right_init)-log1p(-eta_right_init);

    vector < vector <double> > psi(2, vector <double> (N+1,0.)), psi_new(2, vector <double> (N+1,0.));
  
    for(int iter_E_Rc = 0; iter_E_Rc <= iter_E_Rc_max; ++iter_E_Rc)
    {
        #pragma omp parallel num_threads(2)
        {
            #pragma omp for
            for(int i = 0; i <= 1; ++i)
            {
                for(int iter_y = 0; iter_y <= iter_y_max; ++iter_y)
                {
                    for(int iter_psi=0; iter_psi<=iter_psi_max; ++iter_psi)
                    {
                        delta_psi[i] = 0.0;
                        for(int s = 0; s <= N; ++s)
                        {
                            psi_new[i][s] = psi[i][s]-tau_psi*Gradient_psi(s,psi[i],y[i],((i==0)? E_Rc : -E_Rc));
                            delta_psi[i] = max(delta_psi[i],fabs(psi_new[i][s]-psi[i][s]));
                        }
                      
                        psi[i] = psi_new[i];
                      
                        if(delta_psi[i] <= tol_psi)
                            break;
                    }
                  
                    y_new[i] = y[i]+tau_y*Gradient_y(psi[i],y[i]);
                    delta_y[i] = fabs(y_new[i]-y[i]);
                  
                    y[i] = y_new[i];
                  
                    if(delta_y[i] <= tol_y)
                        break;
                }
            }
        }
      
        E_Rc_new = E_Rc+tau_E_Rc*(Gradient_E_Rc(psi[0])-Gradient_E_Rc(psi[1]));
        delta_E_Rc = fabs(E_Rc_new-E_Rc);
      
        E_Rc = E_Rc_new;
      
        if(delta_E_Rc <= tol_E_Rc)
            break;
    }

//calculating Omega [compared to notes, the extra minus sign in each term comes due to conversion from fem to true potential]
     
    double Omega = 0.;
    
    for(int v = 0; v <= N-1; ++v)
    {
        double dpsi = psi[0][v+1] - psi[0][v];
        
        double t1_c = ((fabs(dpsi)>0.0000003)? sinh(dpsi)/dpsi : 1.);
        double t1_s = ((fabs(dpsi)>0.001)? (cosh(dpsi)-1.)/dpsi : 0.5*dpsi);
        double t2_c = ((fabs(dpsi)>0.001)? (dpsi*sinh(dpsi)-cosh(dpsi)+1.)/(dpsi*dpsi) : 0.5);
        double t2_s = ((fabs(dpsi)>0.001)? (dpsi*cosh(dpsi)-sinh(dpsi))/(dpsi*dpsi) : dpsi/3.);
        double t3_c = ((fabs(dpsi)>0.001)? ((dpsi*dpsi+2.)*sinh(dpsi)-2.*dpsi*cosh(dpsi))/(dpsi*dpsi*dpsi) : 1./3.);
        double t3_s = ((fabs(dpsi)>0.01)? ((dpsi*dpsi+2.)*cosh(dpsi)-2.*dpsi*sinh(dpsi)-2.)/(dpsi*dpsi*dpsi) : 0.25*dpsi);
        
        double T1 = -6.022e-4*8*M_PI*I*dr*(Rm+v*dr)*(Rm+v*dr)*(t1_c*cosh(psi[0][v])+t1_s*sinh(psi[0][v]));
        Omega += T1;
        double T2 = -6.022e-4*16*M_PI*I*dr*dr*(Rm+v*dr)*(t2_c*cosh(psi[0][v])+t2_s*sinh(psi[0][v]));
        Omega += T2;
        double T3 = -6.022e-4*8*M_PI*I*dr*dr*dr*(t3_c*cosh(psi[0][v])+t3_s*sinh(psi[0][v]));
        Omega += T3;
        double T4 = -(ep/(2.*l_B))*(Rm*Rm/dr+Rm*(2*v+1)+v*(v+1)+dr/3.)*dpsi*dpsi;
        Omega += T4;
    }
    
    for(int v = 0; v <= N-1; ++v)
    {
        double dpsi = psi[1][v+1] - psi[1][v];
        
        double t1_c = ((fabs(dpsi)>0.0000003)? sinh(dpsi)/dpsi : 1.);
        double t1_s = ((fabs(dpsi)>0.001)? (cosh(dpsi)-1.)/dpsi : 0.5*dpsi);
        double t2_c = ((fabs(dpsi)>0.001)? (dpsi*sinh(dpsi)-cosh(dpsi)+1.)/(dpsi*dpsi) : 0.5);
        double t2_s = ((fabs(dpsi)>0.001)? (dpsi*cosh(dpsi)-sinh(dpsi))/(dpsi*dpsi) : dpsi/3.);
        double t3_c = ((fabs(dpsi)>0.001)? ((dpsi*dpsi+2.)*sinh(dpsi)-2.*dpsi*cosh(dpsi))/(dpsi*dpsi*dpsi) : 1./3.);
        double t3_s = ((fabs(dpsi)>0.01)? ((dpsi*dpsi+2.)*cosh(dpsi)-2.*dpsi*sinh(dpsi)-2.)/(dpsi*dpsi*dpsi) : 0.25*dpsi);
        
        double T1 = -6.022e-4*8*M_PI*I*dr*(Rm+v*dr)*(Rm+v*dr)*(t1_c*cosh(psi[1][v])+t1_s*sinh(psi[1][v]));
        Omega += T1;
        double T2 = -6.022e-4*16*M_PI*I*dr*dr*(Rm+v*dr)*(t2_c*cosh(psi[1][v])+t2_s*sinh(psi[1][v]));
        Omega += T2;
        double T3 = -6.022e-4*8*M_PI*I*dr*dr*dr*(t3_c*cosh(psi[1][v])+t3_s*sinh(psi[1][v]));
        Omega += T3;
        double T4 = -(ep/(2.*l_B))*(Rm*Rm/dr+Rm*(2*v+1)+v*(v+1)+dr/3.)*dpsi*dpsi;
        Omega += T4;
    }
    
    double T5_l = 4*M_PI*(eta(y[0])-1/theta)*Rm*Rm*psi[0][0]/(a*a);
    Omega += T5_l;
    double T5_r = 4*M_PI*(eta(y[1])-1/theta)*Rm*Rm*psi[1][0]/(a*a);
    Omega += T5_r;
    
    double T6_l = -ep*E_Rc*psi[0][N];
    Omega += T6_l;
    double T6_r = ep*E_Rc*psi[1][N];
    Omega += T6_r;

    double T7_l = -4*M_PI*Rm*Rm*(alpha*eta(y[0]) + 0.5*chi*eta(y[0])*eta(y[0]) - eta(y[0])*log(eta(y[0])) - (1-eta(y[0]))*log(1-eta(y[0])))/(a*a);
    Omega += T7_l;
    double T7_r = -4*M_PI*Rm*Rm*(alpha*eta(y[1]) + 0.5*chi*eta(y[1])*eta(y[1]) - eta(y[1])*log(eta(y[1])) - (1-eta(y[1]))*log(1-eta(y[1])))/(a*a);
    Omega += T7_r;
    
//printing results
  
    cout.precision(numeric_limits<double>::digits10);
    cout << "# ep = " << ep << endl;
    cout << "# I = " << I << endl;
    cout << "# Rm = " << Rm << endl;
    cout << "# dr = " << dr << endl;
    cout << "# N = " << N << endl;
    cout << "# l_B = " << l_B << endl;
    cout << "# tol_psi = " << tol_psi << endl;
    cout << "# iter_psi_max = " << iter_psi_max << endl;
    cout << "# tau_psi = " << tau_psi << endl;
    cout << "# delta_psi_left = " << delta_psi[0] << endl;
    cout << "# delta_psi_right = " << delta_psi[1] << endl;
    cout << "# tol_y = " << tol_y << endl;
    cout << "# iter_y_max = " << iter_y_max << endl;
    cout << "# tau_y = " << tau_y << endl;
    cout << "# delta_y_left = " << delta_y[0] << endl;
    cout << "# delta_y_right = " << delta_y[1] << endl;
    cout << "# tol_E_Rc = " << tol_E_Rc << endl;
    cout << "# iter_E_Rc_max = " << iter_E_Rc_max << endl;
    cout << "# tau_E_Rc = " << tau_E_Rc << endl;
    cout << "# delta_E_Rc = " << delta_E_Rc << endl;
    cout << "# alpha = " << alpha << endl;
    cout << "# chi = " << chi << endl;
    cout << "# theta = " << theta << endl;
    cout << "# eta_left_initial  = " << eta_left_init << endl;
    cout << "# eta_right_initial = " << eta_right_init << endl;
    cout << "#" << endl;
    cout << "# eta_left  = " << eta(y[0]) << endl;
    cout << "# eta_right = " << eta(y[1]) << endl;
    cout << "# E(R_c) = " << E_Rc << endl;
    cout << "# Omega = " << Omega << endl;
    cout << "#" << endl;
    cout << "# s psi_left(s) psi_right(s)" << endl;

    for(int s=0; s<=N; ++s)
    {
        cout << s << " " << psi[0][s] << " " << psi[1][s] << endl;
    }
    return 0;
}
