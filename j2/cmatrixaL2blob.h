/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

  

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
#include <vector>
#include <algorithm>

using namespace std;

 int name_hash(const vector<int> &a) ;


template<class T>
class CompMatrixA: public MatrixWithProduct<std::complex<T> > {

 private:

	int n;
	void Convert(int L, int j, vector<int>& w, int i);
	void Convert2(int L, int i, int j, vector<int>& w, int r, int limit);

 
                //


 public:


 
  void MultMv(std::complex<T>* v, std::complex<T>* w, vector<vector<int>> C, int L, double loopweight, double loopweightblob);

  CompMatrixA(int n);

}; // CompMatrixA.






template<class T>
void CompMatrixA<T>::Convert(int L, int j, vector<int>& w, int i)
{

	if (j!=0){
	
		int k;
		k=j%(2*L);
		w[2*L-i-1]=k;
		
		j=j/(2*L);
        	i=i+1;
		Convert(L, j, w, i);

	}
}




template<class T>
void CompMatrixA<T>::MultMv(std::complex<T>* v, std::complex<T>* w, vector<vector<int>> C, int L, double loopweight, double loopweightblob)
/*
  Matrix-vector subroutine. Computes w <- M*v.
*/

{



int r;
//double x1=-0.414213562373095;
//double x1=1;
double x1;
x1=(sqrt(4-(loopweight)*(loopweight))-2)/(loopweight);
double x2=x1;
vector<int> z (2*L+2,0);
vector<int> y (2*L+2,0);



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




for(int i=0;i<(C.size());i++)
{
VTable[C[i]].vectcomp=v[i];
VTable[C[i]].loop=C[i];
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

						
for(int i=(2*L-1); i>0; i=i-2){					//loop over even Xi (even when considering first position as position 1, not 0)
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

for(int i=(2*L); i>=0; i=i-2){					//loop over odd Xi (odd when considering first position as position 1, not 0)
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


/*cout << WTable.size()<<endl;
cin>>r;
*/
for(int i=0;i<(C.size());i++)
{
w[i]=WTable[C[i]].vectcomp;
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


#endif // CMATRIXA_H
