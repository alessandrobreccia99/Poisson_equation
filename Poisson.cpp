//
//  main.cpp
//  ProbSfere
//
//  Created by Alessandro Breccia on 11/03/21.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <ccomplex>

using namespace std;

void ResetPiastre (double *Ex,double *Ey,int colc1,int colc2,int N,int l);
void ResetVecInt (double *Ex,double *Ey,int colc1,int colc2,int N,int l,double dcol);
void NormVec (double *Ex,double *Ey,double maxx,double maxy,double maxt,int N,double dcol);
void CampoEx (double *Ex,double *V,int N,double me);
void CampoEy (double *Ey,double *V,int N,double me);
void MesTer1 (double *V,int col,int rig,int N,double r2);
void MesTer2 (double *V,int col1,int rig1,int N,double r2);
void MesTerC (double *V,int colc1,int colc2,int N);
void Pote (double *V,double *V1,double *f,int N,double mt,double tempq,double temp,double diff,double thdiff);
void RiempPiastre(double *f,int N,double ro0,double ro1,int colc1,int colc2,int l,double mll);
void RiempSfere(double *f,int N,double ro0,double ro1,int col,int col1,int rig,int rig1,double r2_1,double r2_2,double mll,char intsf1,char intsf2);
void StampaPiastre(int N,double *V,double *Ex,double *Ey);
void StampaSfere(int N,double *V,double *Ex,double *Ey);
void InizM1 (double *V,double *f,int N);
void InizF (double *f,int N);

int main() {
    
    int N=100,ro0,ro1,col,rig,col1,rig1,colc1,colc2,l;
    int n=N*N;
    double *Ex=new double[n];
    double *Ey=new double[n];
    double *V=new double[n];
    double *V1=new double[n];
    double *f=new double[n];
    double h=1;
    double mud=1./pow(h,2.);
    double mll=-2.0*(1./pow(h,2)+1./pow(h,2));
    double mt=-mud/mll;
    double diff=10,temp=0,tempq=0,thdiff=1;
    double me=1/(2*h);
    double r2_1t=0,r2_2t=0,r2_1=0,r2_2=0;
    double maxx=0,maxy=0,maxt=0;
    char intsf1,intsf2;
    char tyPr;
    cout<<"Condensatori a piastre (c) o Sfere (s): "<<'\n';
    cin>>tyPr;
    
    if (tyPr=='s') {
        cout<<"Inserisci colonna e riga 1a sfera [0,99]:"<<'\n';
        cin>>col;
        cin>>rig;
        cout<<"Inserisci raggio 1a sfera:"<<'\n';
        cin>>r2_1t;
        cout<<"1a Sfera internamente carica (c) o messa a terra (t):"<<'\n';
        cin>>intsf1;
        cout<<"Carica esterna 1a sfera (-100C : +100C)"<<'\n';
        cin>>ro0;
        cout<<"Inserisci colonna e riga 2a sfera [0,99]:"<<'\n';
        cin>>col1;
        cin>>rig1;
        cout<<"Inserisci raggio 2a sfera:"<<'\n';
        cin>>r2_2t;
        cout<<"2a Sfera internamente carica (c) o messa a terra (t):"<<'\n';
        cin>>intsf2;
        cout<<"Carica esterna 2a sfera (-100C : +100C)"<<'\n';
        cin>>ro1;
    }
    else if (tyPr=='c') {
        cout<<"Inserisci colonna 1° condensatore [0,99]: "<<'\n';
        cin>>colc1;
        cout<<"Inserisci colonna 2° condensatore [0,99]: "<<'\n';
        cin>>colc2;
        cout<<"Inserisci lunghezza piastre [0:99]: "<<'\n';
        cin>>l;
        cout<<"Carica 1a piastra (-100C : +100C)"<<'\n';
        cin>>ro0;
        cout<<"Carica 2a piastra (-100C : +100C)"<<'\n';
        cin>>ro1;
    }
    double dcol=colc2-colc1;
    r2_1=r2_1t*r2_1t;
    r2_2=r2_2t*r2_2t;
// Matrice f
    InizF(f, N);
    // Riempimento Sfere
   if (tyPr=='s') {
       RiempSfere(f,N,ro0,ro1,col,col1,rig,rig1,r2_1,r2_2,mll,intsf1,intsf2);
   }
    // Riempimento Piastre
   else if(tyPr=='c'){
       RiempPiastre(f, N, ro0, ro1, colc1, colc2, l,mll);
   }

   // Matrice M
        InizM1(V, f, N);
   // Matrice Potenziale
        Pote(V, V1, f, N, mt, tempq, temp, diff,thdiff);
    
    if (intsf1=='t') {
        MesTer1(V, col, rig, N, r2_1);
    }
    if (intsf2=='t') {
        MesTer2(V, col1, rig1, N, r2_2);
    }
   // if (tyPr=='c'){
   //     MesTerC (V,colc1,colc2,N);
   // }
    
    
// Matrice Campo
       CampoEx(Ex, V, N, me);
    
       CampoEy(Ey, V, N, me);
    
       NormVec(Ex, Ey, maxx, maxy, maxt, N, dcol);
   if (tyPr=='c') {
       ResetVecInt(Ex, Ey, colc1, colc2, N, l, dcol);
    
       ResetPiastre(Ex, Ey, colc1, colc2, N, l);
   }
    cout<<"ok"<<endl;
   // Stampa Risultati
    // V
    if (tyPr=='s') {
        StampaSfere(N, V, Ex, Ey);
    }
        else if(tyPr=='c'){
            StampaPiastre( N, V, Ex, Ey);
            }
    return 0;
}


void ResetPiastre (double *Ex,double *Ey,int colc1,int colc2,int N,int l){
    for (int i=(N-l)/2; i<=N-(N-l)/2; i++) {
        Ex[i*N+colc1]=0;
        Ex[i*N+colc2]=0;
        Ey[i*N+colc1]=0;
        Ey[i*N+colc2]=0;
    }
}
void ResetVecInt (double *Ex,double *Ey,int colc1,int colc2,int N,int l,double dcol){
    for (int i=0; i<N; i++) {
       for (int j=colc1; j<=colc2; j++) {
           if ( sqrt(pow(Ex[i*N+j],2)+pow(Ey[i*N+j],2))>dcol) {
               Ex[i*N+j]=Ex[i*N+j]/abs(Ex[i*N+j])*dcol/3;
           }
       }
    }
}

void NormVec (double *Ex,double *Ey,double maxx,double maxy,double maxt,int N,double dcol){
    for (int i=0; i<N; i++) {
       for (int j=0; j<N; j++) {
           if (abs(Ex[i*N+j])<abs(Ex[i*N+j+1])) {
               maxx=abs(Ex[i*N+j+1]);
           }
           if (abs(Ey[i*N+j])<abs(Ey[i*N+j+1])) {
               maxy=abs(Ey[i*N+j+1]);
           }
       }
    }
    
    maxt=max(maxx, maxy);
    for (int i=0; i<N; i++) {
       for (int j=0; j<N; j++) {
           Ey[i*N+j]=Ey[i*N+j]/maxt;
           Ex[i*N+j]=Ex[i*N+j]/maxt;
       }
    }
}
void CampoEx (double *Ex,double *V,int N,double me){
    for (int i=0; i<N; i++) {
       for (int j=0; j<N; j++) {
               if (i==0 || i==(N-1) || j==0 || j==(N-1)) {
                   Ex[i*N+j]=0;
               }
              else Ex[i*N+j]=(V[i*N+j+1]-V[i*N+j-1])*me;
       }
    }
}
void CampoEy (double *Ey,double *V,int N,double me){
    for (int i=0; i<N; i++) {
       for (int j=0; j<N; j++) {
               if (i==0 || i==(N-1) || j==0 || j==(N-1)) {
                   Ey[i*N+j]=0;
               }
              else Ey[i*N+j]=(V[(i+1)*N+j]-V[(i-1)*N+j])*me;
       }
    }
}
void MesTer1 (double *V,int col,int rig,int N,double r2_1){
    for (int i=0; i<N; i++) {
      for (int j=0; j<N; j++) {
             if ( floor((i-rig)*(i-rig)+(j-col)*(j-col))<floor(r2_1-36) ) {
                 V[i*N+j]=0;
             }
      }
    }
}
void MesTer2 (double *V,int col1,int rig1,int N,double r2_2){
    for (int i=0; i<N; i++) {
      for (int j=0; j<N; j++) {
             if ( floor((i-rig1)*(i-rig1)+(j-col1)*(j-col1))<floor(r2_2-36) ) {
                 V[i*N+j]=0;
             }
      }
    }
}
void MesTerC (double *V,int colc1,int colc2,int N){
    for (int i=0; i<N; i++) {
      for (int j=0; j<N; j++) {
             if (j<colc1 || j>colc2 ){
                 V[i*N+j]=0;
             }
      }
    }
}
void Pote (double *V,double *V1,double *f,int N,double mt,double tempq,double temp,double diff,double thdiff){
    while (diff>thdiff) {
    for (int i=0; i<N; i++) {
       for (int j=0; j<N; j++) {
           if (i==0 || i==(N-1) || j==0 || j==(N-1)) {
               V1[i*N+j]=V[i*N+j];
           }
           else {
                   V1[i*N+j]=V[i*N+j-1]*mt+V[i*N+j+1]*mt+V[(i-1)*N+j]*mt+V[(i+1)*N+j]*mt+f[i*N+j];
           }
       }
     }
   
       for (int i=0; i<N; i++) {
           for (int j=0; j<N; j++) {
               tempq=pow(V1[i*N+j]-V[i*N+j],2);
               temp+=sqrt(tempq);
           }
       }
       diff=temp;
       temp=0;
        for (int i=0; i<N; i++) {
           for (int j=0; j<N; j++) {
               V[i*N+j]=V1[i*N+j];
           }
       }
    }
}
void RiempPiastre(double *f,int N,double ro0,double ro1,int colc1,int colc2,int l,double mll){
    
    for (int i=(N-l)/2; i<=N-(N-l)/2; i++) {
        f[i*N+colc1]=ro0/mll;
        f[i*N+colc2]=ro1/mll;
    }
}
void RiempSfere(double *f,int N,double ro0,double ro1,int col,int col1,int rig,int rig1,double r2_1,double r2_2,double mll,char intsf1,char intsf2){
    
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            if ( floor((i-rig)*(i-rig)+(j-col)*(j-col))<floor(r2_1) ) {
                f[i*N+j]=ro0/mll;
            }
        }
    }
       if (intsf1=='t') {
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            if ( floor((i-rig)*(i-rig)+(j-col)*(j-col))<floor(r2_1-36) ) {
                f[i*N+j]=0;
            }
        }
    }
       }
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            if ( floor((i-rig1)*(i-rig1)+(j-col1)*(j-col1))<floor(r2_2) ) {
                f[i*N+j]=ro1/mll;
            }
        }
    }
       if (intsf2=='t') {
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            if ( floor((i-rig1)*(i-rig1)+(j-col1)*(j-col1))<floor(r2_2-36) ) {
                f[i*N+j]=0;
            }
        }
    }
       }
}
void StampaSfere(int N,double *V,double *Ex,double *Ey){
    ofstream Pot;
    ofstream Cam;
    Pot.open("/Users/alessandrobreccia/Desktop/3_ANNO_FISICA/Metodi Computazioali/DatiSim/V.txt",ios::out);
        Pot.precision(5);
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
        Pot<< j << "    "<<i<<"   "<<-V[i*N+j]<<'\n';
        }
        Pot<<'\n';
    }
    // E
    Cam.open("/Users/alessandrobreccia/Desktop/3_ANNO_FISICA/Metodi Computazioali/DatiSim/E.txt",ios::out);
        Cam.precision(5);
      for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            if ( sqrt(pow(Ex[i*N+j],2)+pow(Ey[i*N+j],2))<6 ) {
        Cam<< j <<"    "<<i<<"   "<<Ex[i*N+j]<<"    "<<Ey[i*N+j]<<'\n';
            }
        }
        Cam<<'\n';
      }
}
void StampaPiastre(int N,double *V,double *Ex,double *Ey){
    ofstream Pot;
    ofstream Cam;
    Pot.open("/Users/alessandrobreccia/Desktop/3_ANNO_FISICA/Metodi Computazioali/DatiSim/VC.txt",ios::out);
        Pot.precision(5);
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
        Pot<< j << "    "<<i<<"   "<<-V[i*N+j]<<'\n';
        }
        Pot<<'\n';
    }
    // E
    Cam.open("/Users/alessandrobreccia/Desktop/3_ANNO_FISICA/Metodi Computazioali/DatiSim/EC.txt",ios::out);
        Cam.precision(5);
      for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
        Cam<< j <<"    "<<i<<"   "<<Ex[i*N+j]<<"    "<<Ey[i*N+j]<<'\n';
        }
        Cam<<'\n';
      }
}
void InizF (double *f,int N){
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            f[i*N+j]=0;
        }
    }
}
void InizM1 (double *V,double *f,int N){
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            V[i*N+j]=f[i*N+j];
        }
    }
}
