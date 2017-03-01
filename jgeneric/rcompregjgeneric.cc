/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE RCompReg.cc.
   Example program that illustrates how to solve a complex
   standard eigenvalue problem in regular mode using the 
   ARrcCompStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in regular mode,
      where A is obtained from the standard central difference
      discretization of the convection-diffusion operator 
                     (Laplacian u) + rho*(du / dx)
      on the unit squre [0,1]x[0,1] with zero Dirichlet boundary
      conditions.

   2) Data structure used to represent matrix A:

      ARrcCompStdEig is a class thar requires the user to provide
      a way to perform the matrix-vector product w = Av. In this
      example a class called CompMatrixA was created with this
      purpose. CompMatrixA contains a member function, MultMv(v,w),
      that takes a vector v and returns the product Av in w.

   3) The reverse communication interface:

      This example uses the reverse communication interface, which
      means that the desired eigenvalues cannot be obtained directly
      from an ARPACK++ class.
      Here, the overall process of finding eigenvalues by using the
      Arnoldi method is splitted into two parts. In the first, a
      sequence of calls to a function called TakeStep is combined
      with matrix-vector products in order to find an Arnoldi basis.
      In the second part, an ARPACK++ function like FindEigenvectors
      (or EigenValVectors) is used to extract eigenvalues and
      eigenvectors.

   4) Included header files:

      File             Contents
      -----------      -------------------------------------------
      cmatrixa.h       The CompMatrixA class definition.
      arrscomp.h       The ARrcCompStdEig class definition.
      rcompsol.h       The Solution function.
      arcomp.h         The "arcomplex" (complex) type definition.

   5) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include <fstream>
#include <forward_list>
#include <cmath>
#include "arcomp.h"
#include "arrscomp.h"
#include "cmatrixajgeneric.h"
#include "rcompsol.h"

using namespace std;

int name_hash(const vector<int> &a);
void MultMat(vector<vector<int>>& B,vector<vector<int>>& C, int L, double loopweight, double loopweightblob, int defectno);

template<class T>
void Test(T type)
{

ofstream myfile;
myfile.precision(16);
myfile.open("FvsN67r1j4AF.dat");
double gamma;
int r1=1;
double loopweightblob;
int L;
int defectno=4;


for(gamma=0.5;gamma<=1.5;gamma=gamma+0.1)
{
loopweightblob=(sin((r1+1)*gamma))/(sin(r1*gamma));
for(L=6;L<=7;L++)
{
cout<<gamma<<endl;

double loopweight=2*cos(gamma);

int r=0;
vector < vector<int> > B;
vector < vector<int> > C;
vector < vector<int> > D;
vector < vector<int> > E;
vector < vector<int> > G;

vector<int> y (2*L,0);
vector<int> z (2*L,0);
vector<int> w (2*L,0);
vector<int> y2 (2*L,0);
vector<int> z2 (2*L,0);
vector<int> w2 (2*L,0);
vector<int> y3 (2*L,0);
vector<int> z3 (2*L,0);
vector<int> w3 (2*L,0);

B.resize(0);
if(L==3 && defectno==2)
{
y={-1,0,-3,-2,100,100};
z={-1,0,100,100,5,4};
B.push_back(y);
B.push_back(z);
}

else if(L==4 && defectno==2)
{
y={-1,0,-3,-2,-5,-4,100,100};
z={-1,0,100,100,5,4,7,6};
w={-1,0,-3,-2,5,4,100,100};
B.push_back(y);
B.push_back(z);
B.push_back(w);
}

else if(L==5 && defectno==2)
{
y={-1,0,-3,-2,-5,-4,-7,-6,100,100};
z={-1,0,100,100,5,4,7,6,9,8};
w={-1,0,-3,-2,5,4,100,100,9,8};
y2={-1,0,-3,-2,-5,-4,100,100,9,8};
B.push_back(y);
B.push_back(z);
B.push_back(w);
B.push_back(y2);
}

else if(L==6 && defectno==2)
{
y={-1,0,-3,-2,-5,-4,-7,-6,-9,-8,100,100};
z={-1,0,100,100,5,4,7,6,9,8,11,10};
w={-1,0,-3,-2,5,4,100,100,9,8,11,10};
y2={-1,0,-3,-2,-5,-4,100,100,9,8,11,10};
z2={-1,0,-3,-2,-5,-4,-7,-6,100,100,11,10};
B.push_back(y);
B.push_back(z);
B.push_back(w);
B.push_back(y2);
B.push_back(z2);
}

else if(L==7 && defectno==2)
{
y={-1,0,-3,-2,-5,-4,-7,-6,-9,-8,-11,-10,100,100};
z={-1,0,100,100,5,4,7,6,9,8,11,10,13,12};
w={-1,0,-3,-2,5,4,100,100,9,8,11,10,13,12};
y2={-1,0,-3,-2,-5,-4,100,100,9,8,11,10,13,12};
z2={-1,0,-3,-2,-5,-4,-7,-6,100,100,11,10,13,12};
w2={-1,0,-3,-2,-5,-4,-7,-6,-9,-8,11,10,100,100};
B.push_back(y);
B.push_back(z);
B.push_back(w);
B.push_back(y2);
B.push_back(z2);
B.push_back(w2);
}

else if(L==8 && defectno==2)
{
y={-1,0,-3,-2,-5,-4,-7,-6,-9,-8,-11,-10,-13,-12,100,100};
z={-1,0,100,100,5,4,7,6,9,8,11,10,13,12,15,14};
w={-1,0,-3,-2,5,4,100,100,9,8,11,10,13,12,15,14};
y2={-1,0,-3,-2,-5,-4,100,100,9,8,11,10,13,12,15,14};
z2={-1,0,-3,-2,-5,-4,-7,-6,100,100,11,10,13,12,15,14};
w2={-1,0,-3,-2,-5,-4,-7,-6,-9,-8,-11,-10,100,100,15,14};
y3={-1,0,-3,-2,-5,-4,-7,-6,-9,-8,100,100,13,12,15,14};
B.push_back(y);
B.push_back(z);
B.push_back(w);
B.push_back(y2);
B.push_back(z2);
B.push_back(w2);
B.push_back(y3);
}

if(L==3 && defectno==4)
{
y={-1,0,100,100,100,100};
z={-1,0,-100,100,100,100};
B.push_back(y);
B.push_back(z);
}

else if(L==4 && defectno==4)
{
y={-1,0,3,2,100,100,100,100};
z={-1,0,-3,-2,100,100,100,100};
w={-1,0,-3,-2,-100,100,100,100};
B.push_back(y);
B.push_back(z);
B.push_back(w);
}

else if(L==5 && defectno==4)
{
y={-1,0,-3,-2,-5,-4,100,100,100,100};
z={-1,0,100,100,100,100,7,6,9,8};
w={-1,0,-3,-2,5,4,100,100,100,100};
y2={-1,0,-3,-2,-5,-4,-100,100,100,100};
B.push_back(y);
B.push_back(z);
B.push_back(w);
B.push_back(y2);
}

else if(L==6 && defectno==4)
{
y={-1,0,-3,-2,-5,-4,-7,-6,-100,100,100,100};
z={-1,0,100,100,100,100,7,6,9,8,11,10};
w={-1,0,-3,-2,5,4,100,100,100,100,11,10};
y2={-1,0,-3,-2,-5,-4,100,100,100,100,11,10};
z2={-1,0,-3,-2,-5,-4,-7,-6,100,100,100,100};
B.push_back(y);
B.push_back(z);
B.push_back(w);
B.push_back(y2);
B.push_back(z2);
}

else if(L==7 && defectno==4)
{
y={-1,0,-3,-2,-5,-4,-7,-6,-9,-8,-100,100,100,100};
z={-1,0,100,100,100,100,7,6,9,8,11,10,13,12};
w={-1,0,-3,-2,5,4,100,100,100,100,11,10,13,12};
y2={-1,0,-3,-2,-5,-4,100,100,100,100,11,10,13,12};
z2={-1,0,-3,-2,-5,-4,-7,-6,100,100,100,100,13,12};
w2={-1,0,-3,-2,-5,-4,-7,-6,-9,-8,100,100,100,100};
B.push_back(y);
B.push_back(z);
B.push_back(w);
B.push_back(y2);
B.push_back(z2);
B.push_back(w2);
}

else if(L==8 && defectno==2)
{
y={-1,0,-3,-2,-5,-4,-7,-6,-9,-8,-11,-10,-100,100,100,100};
z={-1,0,100,100,100,100,7,6,9,8,11,10,13,12,15,14};
w={-1,0,-3,-2,5,4,100,100,100,100,11,10,13,12,15,14};
y2={-1,0,-3,-2,-5,-4,100,100,100,100,11,10,13,12,15,14};
z2={-1,0,-3,-2,-5,-4,-7,-6,100,100,100,100,13,12,15,14};
w2={-1,0,-3,-2,-5,-4,-7,-6,-9,-8,-11,-10,100,100,100,100};
y3={-1,0,-3,-2,-5,-4,-7,-6,-9,-8,100,100,100,100,15,14};
B.push_back(y);
B.push_back(z);
B.push_back(w);
B.push_back(y2);
B.push_back(z2);
B.push_back(w2);
B.push_back(y3);
}



C.resize(0);
int csize=0;

MultMat(B,C,L,loopweight,loopweightblob, defectno);


if(C.size()!=csize)
{
while(C.size()!=csize)
{
csize=C.size();

//cout << "csize = " << csize <<endl;
//cout << "C.size() = "<< C.size() <<endl;

C.resize(0);

MultMat(B,C,L,loopweight,loopweightblob, defectno);
}

}

/*cout << "csize = " << C.size() << endl;
for(int i=0;i<C.size();i++)
{
	for(int j=0; j<(2*L); j++)
	{
	cout<<C[i][j]<<" ";
	}
cout << endl;
}

cin >> r;
*/
// Defining a complex matrix.


	

  CompMatrixA<T> A(csize); // n = 10*10.

  // Creating a complex eigenvalue problem and defining what we need:
  // the four eigenvectors of A with largest magnitude.
  ARrcCompStdEig<T> prob(csize, 2L);

  // Finding an Arnoldi basis.


	int i=0;
  while (!prob.ArnoldiBasisFound()) {

    // Calling ARPACK FORTRAN code. Almost all work needed to
    // find an Arnoldi basis is performed by TakeStep.

    prob.TakeStep();

    if ((prob.GetIdo() == 1)||(prob.GetIdo() == -1)) {

      // Performing matrix-vector multiplication.
      // In regular mode, w = Av must be performed whenever
      // GetIdo is equal to 1 or -1. GetVector supplies a pointer
      // to the input vector, v, and PutVector a pointer to the
      // output vector, w.

	/*int s;
	s=sizeof(prob.GetVector());
	complex<double>* v[s];	
	v[0]=prob.GetVector();*/

      //A.MultMv(prob.GetVector(), prob.PutVector()); This line was here originally, revert to this if necessary
	
	A.MultMv(prob.GetVector(), prob.PutVector(), C, L, loopweight, loopweightblob);
	//cout << " y= " << loopweightblob <<endl;
	//i++;
	
    }

  }

  // Finding eigenvalues and eigenvectors.

  prob.FindEigenvalues();

  // Printing solution.

  Solution(prob);

double F;
F=-(log(abs(prob.Eigenvalue(1).real())))/(2*L);
if(L==3)
{
myfile << gamma << "  ";
myfile << loopweightblob << "	";
}
myfile <<F<< "   ";

}				//End of for loop over L
myfile<< "\n";
}				//End of loop over loopweightblob

myfile.close(); 


} // Test.


int main()
{

  // Solving a single precision problem with n = 100.

/*#ifndef __SUNPRO_CC

  Test((float)0.0);

#endif*/






  // Solving a double precision problem with n = 100.
cout.precision(16);
  Test((double)0.0);

} // main


void MultMat(vector<vector<int>>& B,vector<vector<int>>& C, int L, double loopweight, double loopweightblob, int defectno)
/*
  Matrix-vector subroutine. Computes w <- M*v.
*/

{



int r;
//double x1=-0.414213562373095;
double x1=1;
double x2=x1;
vector<int> z (2*L,0);
vector<int> y (2*L,0);



int p=1709;
struct node                                   //
{                                             //
    std::complex<double> vectcomp;                 //new data structure that will be inserted into hash table                                //
    vector<int> loop;                             //
}; 

std::unordered_map<vector<int>, node, decltype(&name_hash)> VTable (p, name_hash);           //
std::unordered_map<vector<int>, node, decltype(&name_hash)> WTable (p, name_hash);           //
std::unordered_map<vector<int>, node, decltype(&name_hash)> addWTable (p, name_hash);        // Setting up the required hash tables

typename std::unordered_map<vector<int>, node, decltype(&name_hash)>::iterator it;                    // Setting up the required iterators
typename std::unordered_map<vector<int>, node, decltype(&name_hash)>::iterator it2;                   //




for(int i=0;i<(B.size());i++)
{
VTable[B[i]].vectcomp=1;
VTable[B[i]].loop=B[i];
}



/*
for(it=VTable.begin();it!=VTable.end();it++)
{
for(int i=0; i<2*L;i++)
{
cout<<VTable[it->first].loop[i];
}
cout<<VTable[it->first].vectcomp<<endl;

}
 cin>>r;
*/

/*for(int i=0;i<(C.size());i++)
{
cout << "v" << i << "=" << v[i] << endl;
}

cin >> r;
*/




/*node w1;							//
w1.vectcomp={1,0}; 						//
vector<int> looptemp(w1a,w1a+sizeof(w1a)/sizeof(int));		//insert initial node into WTable
w1.loop = looptemp;						//
WTable.insert({w1.loop, w1});					//
*/

						
for(int i=(2*L-3); i>0; i=i-2){					//loop over even Xi (even when considering first position as position 1, not 0)
for (it=VTable.begin(); it!=VTable.end(); it++){		//loop over components of input vector


if (VTable[it->first].loop[i]==i+1 || VTable[it->first].loop[i]==-(i+1)){
	if(VTable[it->first].loop[i]<=0 || VTable[it->first].loop[i+1]<=0)
		{
		WTable[it->first].vectcomp=VTable[it->first].vectcomp;
		WTable[it->first].loop=VTable[it->first].loop;

		z=WTable[it->first].loop;
		z[i]=abs(z[i]);
		z[i+1]=abs(z[i+1]);

		it2=addWTable.find(z);
		if(it2==addWTable.end())
			{
				addWTable[z].vectcomp=0;
			}

		addWTable[z].vectcomp=addWTable[z].vectcomp+x2*loopweightblob*VTable[it->first].vectcomp;
		addWTable[z].loop=z;
		}		
	else
		{
		WTable[it->first].vectcomp=VTable[it->first].vectcomp+x1*loopweight*VTable[it->first].vectcomp;
		WTable[it->first].loop=VTable[it->first].loop;
		}

}									//end of loop[i]=i+1 if statement


else
{
WTable[it->first].vectcomp=VTable[it->first].vectcomp;
WTable[it->first].loop=VTable[it->first].loop;

z=WTable[it->first].loop;

if (abs(z[i])==100 && abs(z[i+1])==100)
{
z[i]=i+1;
z[i+1]=i;
}

else			
{
if(abs(z[i+1])==100)
{
if(z[i]<=0 || z[i+1]<=0)
{
z[abs(z[i])]=-abs(z[i+1]);
}

else
z[z[i]]=z[i+1];
}

else if (abs(z[i])==100)
{

if(z[i]<=0 || z[i+1]<=0)
{
z[abs(z[i+1])]=-abs(z[i]);
}

else
z[z[i+1]]=z[i];

}


else
{
if(z[i]<=0 || z[i+1]<=0)
{
z[abs(z[i+1])]=-abs(z[i]);
z[abs(z[i])]=-abs(z[i+1]);
}

else
{
z[z[i+1]]=z[i];
z[z[i]]=z[i+1];
}

}

z[i]=i+1;
z[i+1]=i;

}			//end of else if not both 100




it2=addWTable.find(z);

if(it2==addWTable.end()){
addWTable[z].vectcomp=0;
}

addWTable[z].vectcomp=addWTable[z].vectcomp+x1*VTable[it->first].vectcomp;
addWTable[z].loop=z;						
}								//end of else statement

}								//end of input components loop

for(it=addWTable.begin();it!=addWTable.end();it++)					//
{											//
it2=WTable.find(it->first);								//Adding addW components
if(it2==WTable.end())									//
{											//
WTable[it->first].loop=addWTable[it->first].loop;					//
}											//
WTable[it->first].vectcomp=WTable[it->first].vectcomp+addWTable[it->first].vectcomp;	//
}											//

for(it=WTable.begin();it!=WTable.end();it++)				//
{									//
it2=VTable.find(it->first);						//
if(it2==VTable.end())							//Sets V=W
{									//
VTable[it->first].loop=WTable[it->first].loop;				//
}									//
VTable[it->first].vectcomp=WTable[it->first].vectcomp;			//
									//
}									//


for (it=addWTable.begin(); it!=addWTable.end(); it++){              //
    addWTable[it->first].vectcomp=0;                          	    // initialises addW to zero
} 								    //


/*for(it=WTable.begin();it!=WTable.end();it++)			//TEST
{
cout << "i = " << i << endl;
cout << WTable[it->first].vectcomp << endl;
for(int i=0;i<2*L+2;i++)
{
cout << WTable[it->first].loop[i]<<" ";
}
cout << "\n" << endl;
}
cin >> r;
*/

}								//end of even i loop

for(int i=(2*L-2); i>=0; i=i-2){					//loop over odd Xi (odd when considering first position as position 1, not 0)
for (it=VTable.begin(); it!=VTable.end(); it++){		//loop over components of input vector


if (VTable[it->first].loop[i]==i+1 || VTable[it->first].loop[i]==-(i+1)){
	if(VTable[it->first].loop[i]<=0 || VTable[it->first].loop[i+1]<=0)
		{
		WTable[it->first].vectcomp=x1*VTable[it->first].vectcomp;
		WTable[it->first].loop=VTable[it->first].loop;
		
		z=WTable[it->first].loop;
		if(i!=0){
		z[i]=abs(z[i]);
		z[i+1]=abs(z[i+1]);
		}

		it2=addWTable.find(z);
		if(it2==addWTable.end())
			{
				addWTable[z].vectcomp=0;
			}
		addWTable[z].vectcomp=addWTable[z].vectcomp+loopweightblob*VTable[it->first].vectcomp;
		addWTable[z].loop=z;
		}		
	else
		{
		WTable[it->first].vectcomp=x1*VTable[it->first].vectcomp+loopweight*VTable[it->first].vectcomp;
		WTable[it->first].loop=VTable[it->first].loop;
		}

}									//end of loop[i]=i+1 if statement

else
{
WTable[it->first].vectcomp=x2*VTable[it->first].vectcomp;
WTable[it->first].loop=VTable[it->first].loop;

z=WTable[it->first].loop;

if (abs(z[i])==100 && abs(z[i+1])==100)
{
z[i]=i+1;
z[i+1]=i;
}

else			
{
if(abs(z[i+1])==100)
{
if(z[i]<=0 || z[i+1]<=0)
{
z[abs(z[i])]=-abs(z[i+1]);
}

else
z[z[i]]=z[i+1];
}

else if (abs(z[i])==100)
{

if(z[i]<=0 || z[i+1]<=0)
{
z[abs(z[i+1])]=-abs(z[i]);
}

else
z[z[i+1]]=z[i];

}


else
{
if(z[i]<=0 || z[i+1]<=0)
{
z[abs(z[i+1])]=-abs(z[i]);
z[abs(z[i])]=-abs(z[i+1]);
}

else
{
z[z[i+1]]=z[i];
z[z[i]]=z[i+1];
}

}

z[i]=i+1;
z[i+1]=i;

}			//end of else if not both 100

if(i==0)
z[0]=-z[0];


it2=addWTable.find(z);

if(it2==addWTable.end()){
addWTable[z].vectcomp=0;
}
addWTable[z].vectcomp=addWTable[z].vectcomp+VTable[it->first].vectcomp;
addWTable[z].loop=z;
}
						

}								//end of input components loop

for(it=addWTable.begin();it!=addWTable.end();it++)					//
{											//
it2=WTable.find(it->first);								//Adding addW components
if(it2==WTable.end())									//
{											//
WTable[it->first].loop=addWTable[it->first].loop;					//
}											//
WTable[it->first].vectcomp=WTable[it->first].vectcomp+addWTable[it->first].vectcomp;	//
}											//

for(it=WTable.begin();it!=WTable.end();it++)				//
{									//
it2=VTable.find(it->first);						//
if(it2==VTable.end())							//Sets V=W
{									//
VTable[it->first].loop=WTable[it->first].loop;				//
}									//
VTable[it->first].vectcomp=WTable[it->first].vectcomp;			//
									//
}									//

for (it=addWTable.begin(); it!=addWTable.end(); it++){              //
    addWTable[it->first].vectcomp=0;                          	    // initialises addW to zero
} 								    //

/*for(it=WTable.begin();it!=WTable.end();it++)			//TEST
{
cout << "i = " << i << endl;
cout << WTable[it->first].vectcomp << endl;
for(int i=0;i<2*L;i++)
{
cout << WTable[it->first].loop[i];
}
cout << "\n" << endl;
}
cin >> loopweight;
*/

/*for(it=WTable.begin();it!=WTable.end();it++)			//TEST
{
cout << "i = " << i << endl;
cout << WTable[it->first].vectcomp << endl;
for(int i=0;i<2*L+2;i++)
{
cout << WTable[it->first].loop[i]<<" ";
}
cout << "\n" << endl;
}
cin >> r;
*/

}								//end of odd i loop


//cout << WTable.size()<<endl;
//cin>>r;



int defects=0;
int i=0;
for(it=WTable.begin();it!=WTable.end();it++)
{
defects=0;

for(int j=0;j<(2*L);j++)
{
if (abs(WTable[it->first].loop[j])==100)
	{
defects++;
//cout << defects <<endl;
	}

}
if (defects==defectno)
C.push_back(WTable[it->first].loop);
}



B.resize(C.size());
for(int i=0;i<C.size();i++)
{
B[i]=C[i];
}




/*
for(it=WTable.begin();it!=WTable.end();it++)
{
for(int i=0; i<2*L;i++)
{
cout<<WTable[it->first].loop[i];
}

cout<< WTable[it->first].vectcomp<<endl;


}

for(int i=0;i<B.size();i++)
{
cout << "w" << i <<"=" << w[i]<<endl;
}

cin >> r;
*/
  return;

} //  MultMat



int name_hash(const vector<int> &a)                 
{                                           
    int total = 0;
     for(int i=0; i<(sizeof(a)/sizeof(int));i++)
{

	if (a[i]!=-1)
	total=((sizeof(a)/sizeof(int))*total+a[i])%1709;

}
                          
int totalmod;
totalmod=total%1709;   						//Change this depending on number of expected states                 
                                            
 return totalmod;                           
}






