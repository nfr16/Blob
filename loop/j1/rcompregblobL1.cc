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
#include "cmatrixa.h"
#include "rcompsol.h"

using namespace std;

void Permutations(int pos, int L, int open, int close, vector<vector<int>>& B, vector<int>& y, vector<int>& z);
void Convert(int n, int j, vector<int>& w, int i);
void Outloop(int pos, int& max, vector<vector<int>>& B, int i, int L);

template<class T>
void Test(T type)
{

ofstream myfile;
myfile.precision(16);
myfile.open("FvsLAF.dat");
double loopweightblob;
int L;



for(loopweightblob=0;loopweightblob<=2;loopweightblob=loopweightblob+0.1)
{
for(L=3;L<=6;L++)
{


double loopweight=1.414213562373095;

int r=0;
vector < vector<int> > B;
vector < vector<int> > C;
vector < vector<int> > D;
vector < vector<int> > E;

vector<int> y (2*L,0);
vector<int> ynew (2*L+1,0);
vector<int> z (2*L,0);
vector<int> w (2*L,0);

Permutations(0,L,0,0,B,y,z);

/*for(int i=0;i<B.size();i++){
for(int j=0;j<(2*L);j++){
cout << B[i][j]<<" ";
}
cout << endl;
}
cin>>r;
*/


for(int i=0;i<B.size();i++)
{
for(int j=0; j<=2*L; j++)
	{	
		int totalopen=0;
		int totalclose=0;
		for(int k=0;k<j;k++)
		{
			if (B[i][k]==0)
			totalopen++;
			if (B[i][k]==1)
			totalclose++;	
		}
		
		if(j==0||totalopen==totalclose)		
		{
		for(int k=0;k<j;k++)
		{
			ynew[k]=B[i][k];
		}
		
		ynew[j]=100;

		for(int k=j+1;k<=2*L;k++)
		{
			ynew[k]=B[i][k-1];
		}
		E.push_back(ynew);
		}
	}
}


/*for(int i=0;i<E.size();i++){
for(int j=0;j<(2*L+1);j++){
cout << E[i][j]<<" ";
}
cout << endl;
}
cin>>r;
*/


int tzero=0;
int tone=0;
for(int i=0; i<E.size();i++)
{
tzero=0;
tone=0;
for(int j=0; j<(2*L+1);j++)
	{
		if(E[i][j]==100)			
		E[i][j]=-100;
		if(E[i][j]==-100)			
		break;
		if (E[i][j]==0)
		tzero++;
		if (E[i][j]==1)
		tone++;
		if(tone==tzero && E[i][j]!=100)
		E[i][j]=-1;
		if(tone==tzero)
		break;
	}

}


/*for(int i=0;i<E.size();i++){
for(int j=0;j<(2*L+1);j++){
cout << E[i][j]<<" ";
}
cout << endl;
}
cin>>r;
*/

int esize=E.size();
int combnum;
int pos;
int outloopnum=0;
int max;
int posinit=0;
int cond1=0;
int cond2=0;
for(int i=0; i<esize;i++)
{
for(int k=0; k<(2*L+1);k++)
{
  ynew[k]=E[i][k];
}
outloopnum=0;
if(ynew[posinit]!=-100)
{
Outloop(posinit, max, E, i, L);
if(abs(E[i][max+1])!=100)
pos=max+1;
else
pos=max+2;
}
else
{
pos=1;
max=0;
}

while(max!=2*L)
{
Outloop(pos, max, E, i, L);
outloopnum++;

if(abs(E[i][max+1])!=100)
pos=max+1;
else 
{
pos=max+2;
max=max+1;
}

}

combnum=pow(2, outloopnum);
w.resize(outloopnum);

if(combnum!=1)
{
for(int k=1; k<combnum; k++)
{

for(int j=0; j<outloopnum; j++)
{
w[j]=0;
}
r=0;
Convert(outloopnum, k, w, r);


if(ynew[0]!=-100)
{
for(int j=0; j<outloopnum; j++)
	{
if (w[j]==1)
		{	r=0;
			Outloop(r, max, E, i, L);
			if(ynew[max+1]!=100){
			for(int k1=0; k1<=j; k1++)
				{	
					Outloop(max+1, max, E, i, L);
					if (ynew[max+1]==100 && k1!=j)
					cond1=1;
					if (ynew[max+1]==100 && k1!=j)
					break;
				}
			if(cond1!=1)
			ynew[max]=-1;
			else
			{			
			ynew[max+1]=-100;
			cond1=0;
			}
					    }

		
			else 
			{
				ynew[max+1]=-100;
			}
		}
	}
for(int j=0; j<E.size();j++)
			{
				if (ynew==E[j])
				cond1=2;
			}
if (cond1!=2)
{
E.push_back(ynew);
cond1=3;
}



for(int k1=0; k1<(2*L+1);k1++)
{
  ynew[k1]=E[i][k1];
}

}				//end of if ynew not equal to -100

}				//end of k

}				//end of if combnum statement




}				//end of i loop

for(int i=0; i<E.size();i++)
{
if(abs(E[i][2*L])==100)
{
/*for(int k2=0;k2<(2*L+1);k2++)
{
cout<<E[i][k2]<<" ";
}
cout << "\n";
*/
ynew=E[i];
ynew[2*L]=-ynew[2*L];

/*for(int k2=0;k2<(2*L+1);k2++)
{
cout<<ynew[k2]<<" ";
}
cout << "\n";
*/
for(int j=0; j<E.size();j++)
			{
				if (ynew==E[j])
				{
					cond2=1;
			
				}
			}
if (cond2!=1)
{
E.push_back(ynew);
}
cond2=0;

}
}




/*for(int i=0;i<E.size();i++){
for(int j=0;j<(2*L+1);j++){
cout << E[i][j]<<" ";
}
cout << endl;
}

cout << E.size();
cin>>r;
*/



int zerocount=0;
int onecount=0;
for(int k=0;k<E.size(); k++)
{
C.push_back(E[k]);

	for(int i=0; i<=2*L;i++)
	{
		if(abs(E[k][i])==100)
		C[k][i]=E[k][i];		
				
		if(E[k][i]==0)
		{	zerocount=0;
			onecount=0;
			for(int j=1; j<(2*L+1-i);j++)
			{
				if(E[k][i+j]==0)
				zerocount++;

				if(E[k][i+j]==1 || E[k][i+j]==-1)
				onecount++;

				if(onecount>zerocount)
				{
				  if(E[k][i+j]==1)
					{
						C[k][i+j]=i;
						C[k][i]=i+j;
					}
				  if(E[k][i+j]==-1)
					{
						C[k][i+j]=-i;
						C[k][i]=-(i+j);	
					}
				}	
				if(onecount>zerocount)
				break;				
			}
		}
	}
}



/*for(int i=0;i<C.size();i++){
for(int j=0;j<=(2*L);j++){
cout << C[i][j]<<" ";
}
cout << " "<<i <<endl;
}
cin >> r;

cout<<C.size()<<endl;
cout << "input r"<<endl;
cin >> r;
*/




// Defining a complex matrix.


	int csize=C.size();

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
	//cout << i <<endl;
	//i++;
	
    }

  }

  // Finding eigenvalues and eigenvectors.

  prob.FindEigenvalues();

  // Printing solution.

  Solution(prob);

double F;
F=-(log(abs(prob.Eigenvalue(1).real())))/(2*L);
if(L==4)
myfile << loopweightblob << "	";

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


void Permutations(int pos, int L, int open, int close, vector<vector<int>>& B, vector<int>& y, vector<int>& z)
{

if(close==L){
B.push_back(y);
}
else
{
	if(open>close)
	{
		y[pos]=1;
		Permutations(pos+1, L, open, close+1, B,y,z);
	}
	if(open<L)
	{
		y[pos]=0;
		Permutations(pos+1, L, open+1, close, B,y,z);
	}
}

}

void Outloop(int pos, int& max, vector<vector<int>>& B, int i, int L)
{
int tzero=0;
int tone=0;
for(int j=0; j<(2*L+1-pos);j++)
{if (B[i][pos+j]==0)
tzero++;

if (B[i][pos+j]==1 || B[i][pos+j]==-1 )
tone++;

if (tone==tzero)
max=pos+j;

if (tone==tzero)
break;
}
}

void Convert(int n, int j, vector<int>& w, int i)
{

	if (j!=0){

		if (j%2==0){
			w[n-i-1]=0;
		}
	   else {
	   w[n-i-1]=1;
	   }
		j=j/2;
        i=i+1;
		Convert(n, j, w, i);

	}
}









