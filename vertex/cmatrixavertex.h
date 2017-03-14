/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE CMatrixA.h
   Class template for the nx*nx by nx*nx block tridiagonal matrix

                              | T -I          |
                              |-I  T -I       |
                         OP = |   -I  T       |
                              |        ...  -I|
                              |           -I T|

   derived from the standard central difference discretization
   of the 2 dimensional convection-diffusion operator
                       (Laplacian u) + rho*(du/dx) 
   on a unit square with zero boundary condition.
   T is a nx by nx tridiagonal matrix with DD on the diagonal,
   DL on the subdiagonal, and DU on the superdiagonal.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef CMATRIXA_H
#define CMATRIXA_H

#include "arcomp.h"
#include "blas1c.h"
#include "arcomp.h"
#include "matprod.h"
#include <unordered_map>
#include <complex>


 int name_hash(const int &a) ;


template<class T>
class CompMatrixA: public MatrixWithProduct<std::complex<T> > {

 private:

	int n;

  void Convert(int L, int j, int w[], int i);
                //


 public:


 
  void MultMv(std::complex<T>* v, std::complex<T>* w);

  CompMatrixA(int n);

}; // CompMatrixA.




template<class T>
void CompMatrixA<T>::Convert(int L, int j, int w[], int i)
{

	if (j!=0){

		if (j%2==0){
			w[2*L-i-1]=0;
		}
	   else {
	   w[2*L-i-1]=1;
	   }
		j=j/2;
        i=i+1;
		Convert(L, j, w, i);

	}
}







template<class T>
void CompMatrixA<T>::MultMv(std::complex<T>* v, std::complex<T>* w)
/*
  Matrix-vector subroutine. Computes w <- M*v.
*/

{

struct node                                   //
{                                             //
    std::complex<double> vectcomp;                 //new data structure that will be inserted into hash table
    int code;                                 //
    int spin [4];                             //
};                                            //





                    
	double x1=-0.158384440325;                  //
	double x2=x1;                               //

                                                //
	double a, b, c;                             //
	std::complex<double> q;                          //
	std::complex<double> qinverse;                   // coupling constants
	q={0.309016994375, 0.951056516295};         //
	a=0.309016994375;                           //
	b=0.951056516295;                           //
	qinverse={a/((a*a)+(b*b)), -b/((a*a)+(b*b))};   //

	int  i,r;
	
	int L=2;				//
	int powerofL;			//defining system size
	powerofL = pow(2,2*L);	//

	int z[2*L];             //spin vector

int p = 5;             //number of buckets in hash table, must change parameter in hash function to be this value too

std::unordered_map<int, node, decltype(&name_hash)> VTable (p, name_hash);           //
std::unordered_map<int, node, decltype(&name_hash)> WTable (p, name_hash);           //
std::unordered_map<int, node, decltype(&name_hash)> addWTable (p, name_hash);        // Setting up the required hash tables

typename std::unordered_map<int, node, decltype(&name_hash)>::iterator it;                    // Setting up the required iterators
typename std::unordered_map<int, node, decltype(&name_hash)>::iterator it2;                   //

int v1a[4]={1,1,0,0};                  //
int v1b[4]={1,0,1,0};                  //Initial spin vectors, for the S=0 sector only v1a is important
/*int v1c[4]={1,0,0,1};                  //
int v1d[4]={0,1,0,1};                  //
*/

int sizev=sizeof(v);
int sizew=sizeof(w);

/*std::cout<<sizew<<","<<sizev<<std::endl;
for(int i=0;i<=20;i++)
{
w[i]=i;
}

for(int i=0;i<=20;i++)
{
std::cout<<w[i]<<std::endl;
}
*/

/*std::complex<double> neww[sizew];


for(i=0;i<=7;i++){
neww[i]=v[i];
}

if(count==1){
for(i=0;i<=7;i++){
neww[i]=0;
}
neww[0]={1,0};
neww[1]={0,0};
neww[2]={0,0};
}
*/


						//
int total;					//
int r1=0;					//
for(int i=0; i<=15;i++){			// Change max i to highest binary number produced (i.e. all 1's on left, all 0's on right)
r=0;	
						//
for(int r3=0; r3<2*L; r3++){			//
z[r3]=0;					//
}						// This inserts the input vectors into VTable

						//
total=0;					//
Convert(L, i, z, r);				//

					//
for(int r2=0; r2<2*L;r2++){			//
total=total+z[r2];				//
}

						//
if (total==L){					//
VTable[i].vectcomp=v[r1];			//
r1++;						//
for(int j=0; j<2*L; j++){			//
VTable[i].spin[j]=z[j];				//
}						//
VTable[i].code=i;				//
}
						//
}						//




/*for(it=VTable.begin();it!=VTable.end();it++)
{
std::cout<<"VTABLE"<<"code= "<< VTable[it->second.code].code<<" "<<VTable[it->second.code].vectcomp<<std::endl;
}*/

/*std::cout << sizeof(w)<<std::endl;

for(i=0;i<sizeof(w);i++){
std::cout << w[i] << std::endl;
}
 


std::cout << VTable[3].vectcomp << std::endl;
std::cout << VTable[5].vectcomp << std::endl;
std::cout << VTable[6].vectcomp<< std::endl;
std::cout << VTable[9].vectcomp << std::endl;
std::cout << VTable[10].vectcomp << std::endl;
std::cout << VTable[12].vectcomp << std::endl;




std::cin >> i;

*/



node w1a;
w1a.code=0;
for(int i=0; i<2*L;i++){

    w1a.code=2*w1a.code+v1a[i];

}
 w1a.vectcomp={0,0};
for(i=0; i<2*L; i++){
    w1a.spin[i]=v1a[i];
}

node w2a;
w2a.code=0;
for(int i=0; i<2*L;i++){

    w2a.code=2*w2a.code+v1b[i];

}
 w2a.vectcomp={0,0};
for(i=0; i<2*L; i++){
    w2a.spin[i]=v1b[i];
}

WTable.insert({w1a.code, w1a});





															 
for(i=(2*L-2); i>0; i=i-2){											//loop over even Xi's
for (it=VTable.begin(); it!=VTable.end(); it++){								//loop over components of input vector


		if (it->second.spin[i-1]==1 && it->second.spin[i]==0)  {				//first spin up, second spin down

        WTable[it->first].vectcomp=VTable[it->first].vectcomp+x2*qinverse*VTable[it->first].vectcomp; //Change in component of vector due to this component
                                                                                                                        //

		for(int m=0;m<=(2*L-1);m++){											//initially sets z to spin vector in the relevant node
			z[m]=VTable[it->first].spin[m];								//
		}																		//
		z[i-1]=0;																//changes z to be the vector that the current component of the input vector "also affects"
		z[i]=1;																	//

	int zcode=0;                                                //
	for(int i=0; i<2*L;i++){                                    //finds the key of the new vector
    zcode=2*zcode+z[i];                                         //
}


   addWTable[zcode].code=zcode;                                             //Sets the key of the relevant node in the intermediate hash table
   addWTable[zcode].vectcomp=-(x2)*(VTable[it->first].vectcomp);      //"    "   vectcomp "           "           "           "


    }
	   if (it->second.spin[i-1]==0 && it->second.spin[i]==1){						//first spin down, second spin up

        WTable[it->first].vectcomp=VTable[it->first].vectcomp+x2*q*VTable[it->first].vectcomp; //Change in component of vector due to this component

		for(int m=0;m<=(2*L-1);m++){											//initially sets z to spin vector in the relevant node
			z[m]=VTable[it->first].spin[m];							    //
		}																		//
		z[i-1]=1;																//changes z to be the vector that the current component of the input vector "also affects"
		z[i]=0;																	//

	int zcode=0;                                            //
	for(int i=0; i<2*L;i++){                                //finds the key of the new vector
    zcode=2*zcode+z[i];                                     //
}

   addWTable[zcode].code=zcode;                                             //Sets the key of the relevant node in the intermediate hash table
   addWTable[zcode].vectcomp=-(x2)*(VTable[it->first].vectcomp);      //"         vectcomp "      "           "           "


	   }

	if (it->second.spin[i-1]==it->second.spin[i]){													//Both spins the same
	WTable[it->first].vectcomp=VTable[it->first].vectcomp;                              //Keeps current component the same


}


}																	//end of loop over V hash table



for (it=addWTable.begin(); it!=addWTable.end(); it++){                                                  //
    r=0;                                                                                                //
     it2=WTable.find(it->first) ;                                                                 //
     if(it2==WTable.end()) {                                                                            //
                                                                                                        //
                                                                                                        //
         WTable[it->first].code=addWTable[it->first].code;                                  //  Changes components of W affected by "flipped components"
         WTable[it->first].vectcomp=0;                                                            //
                                                                                                        //
         Convert(L, it->first, WTable[it->first].spin,r);                                   //
     }                                                                                                  //
                                                                                                        //
                                                                                                        //
    WTable[it->first].vectcomp=WTable[it->first].vectcomp+addWTable[it->first].vectcomp;//
                                                                                                        //
}                                                                                                       //




for (it=WTable.begin(); it!=WTable.end(); it++){                                //
        r=0;                                                                    //
     it2=VTable.find(it->first) ;                                         //
     if(it2==VTable.end())  {                                                   //
         VTable[it->first].code=WTable[it->first].code;             // Sets input for next round = output of last round
         Convert(L, it->first, VTable[it->first].spin,r);           //
     }                                                                          //
                                                                                //
    VTable[it->first].vectcomp=WTable[it->first].vectcomp;          //
}                                                                               //



for (it=addWTable.begin(); it!=addWTable.end(); it++){              //
    addWTable[it->first].vectcomp=0;                          // Initialises addw to zero
}                                                                   //


}															        //end of even i loop

/*std::cout<<"after even loop"<<std::endl;
for(it=WTable.begin();it!=WTable.end();it++)
{
std::cout<<"WTABLE"<<WTable[it->first].code<<" " <<it->first<<" "<<WTable[it->first].vectcomp<<std::endl;
}
*/
for(i=(2*L-1); i>0; i=i-2){														//loop over odd Xi's
for (it=VTable.begin(); it!=VTable.end(); it++){												//loop over components of input vector
                                                                 //k keeps track of how many components of the vector we have got through in each mult
                                                                                //ensures that the calculation doesn't bother running for components that are zero
		if (it->second.spin[i-1]==1 && it->second.spin[i]==0)  {									//first spin up, second spin down

        WTable[it->first].vectcomp=(x2)*(VTable[it->first].vectcomp)+qinverse*VTable[it->first].vectcomp;

		for(int m=0;m<=(2*L-1);m++){											//initially sets z to kth column of spin matrix
			z[m]=VTable[it->first].spin[m];														//
		}																		//
		z[i-1]=0;																//changes z to be the vector that the kth component of the input vector "also affects"
		z[i]=1;																	//

	int zcode=0;
	for(int i=0; i<2*L;i++){
    zcode=2*zcode+z[i];
}

   addWTable[zcode].code=zcode;
   addWTable[zcode].vectcomp=-VTable[it->first].vectcomp;

		}
	   if (it->second.spin[i-1]==0 && it->second.spin[i]==1){										//first spin down, second spin up

        WTable[it->first].vectcomp=(VTable[it->first].vectcomp)*(x1)+q*VTable[it->first].vectcomp;

		for(int m=0;m<=(2*L-1);m++){											//initially sets z to kth column of spin matrix
			z[m]=VTable[it->first].spin[m];														//
		}																		//
		z[i-1]=1;																//changes z to be the vector that the kth component of the input vector "also affects"
		z[i]=0;																	//

	int zcode=0;
	for(int i=0; i<2*L;i++){
    zcode=2*zcode+z[i];
}

   addWTable[zcode].code=zcode;
   addWTable[zcode].vectcomp=-VTable[it->first].vectcomp;


	   }

	if (it->second.spin[i-1]==it->second.spin[i]){													//Both spins the same
	WTable[it->first].vectcomp=(x1)*(VTable[it->first].vectcomp);															//	"	"
}


}																	//end of k loop


for (it=addWTable.begin(); it!=addWTable.end(); it++){
    r=0;
     it2=WTable.find(it->first) ;
     if(it2==WTable.end())  {
         WTable[it->first].code=addWTable[it->first].code;
         WTable[it->first].vectcomp=0;

         Convert(L, it->first, WTable[it->first].spin,r);
     }

    WTable[it->first].vectcomp=WTable[it->first].vectcomp+addWTable[it->first].vectcomp;

}



for (it=WTable.begin(); it!=WTable.end(); it++){
        r=0;
     it2=VTable.find(it->first) ;
     if(it2==VTable.end()) {

         VTable[it->first].code=WTable[it->first].code;
         Convert(L, it->first, VTable[it->first].spin,r);
     }

    VTable[it->first].vectcomp=WTable[it->first].vectcomp;

}


for (it=addWTable.begin(); it!=addWTable.end(); it++){
    addWTable[it->first].vectcomp=0;
}

/*cout << WTable[12].vectcomp<<"\n"<<endl;
cin>>r;*/


}															//end of odd i loop




/*for(it=WTable.begin();it!=WTable.end();it++)
{
std::cout<<"WTABLE"<<WTable[it->first].vectcomp<<std::endl;
}*/


int r4=0;
int r3=0;
for (int i=0; i<=15;i++){
r=0;
for(int r4=0; r4<2*L; r4++){			
z[r4]=0;					
}

total=0;
Convert(L,i,z,r);

for(int r2=0; r2<2*L;r2++){			//
total=total+z[r2];				//
}
if(total==(L)){
w[r3]=WTable[i].vectcomp;
r3++;
}

}




/*for(int r3=0; r3<sizeof(w);r3++)
{
std::cout<<w[r3]<<std::endl;

}

*/





  return;

} //  MultMv


template<class T>
CompMatrixA<T>::CompMatrixA(int nval):
  MatrixWithProduct<std::complex<T> >(nval)
/*
  Constructor.
*/

{
n=nval;

} // constructor.



int name_hash(const int &a)
                 //
{                                           //
    int temp = a;                           //Hash function: 5 here is the number of buckets, change if p changes
    int temp2 = temp%5;                     //
                                            //
    return temp2;                           //
}



#endif // CMATRIXA_H
