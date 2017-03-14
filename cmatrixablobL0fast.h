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




struct node                                   //
{                                             //
    std::complex<double> vectcomp;                 //new data structure that will be inserted into hash table                                //
    vector<int> loop;                             //
}; 


template<class T>
void CompMatrixA<T>::MultMv(std::complex<T>* v, std::complex<T>* w, vector<vector<int>> C, int L, double loopweight, double loopweightblob)
/*
  Matrix-vector subroutine. Computes w <- M*v.
*/

{



int r;
//double x1=1;
//double x1=-0.414213562373095;
double x1;
x1=(sqrt(4-(loopweight)*(loopweight))-2)/(loopweight);
double x2=x1;
vector<int> z (2*L,0);
vector<int> y (2*L,0);
vector<int> s  (2*L,0);
int p;




if(L==7)
{
p=1709;
}   



if(L==8)
{
p=6449;
}

                  //


if(L==9)
{
p=24317;
}


                  //


if(L==10)
{
p=92377;
}




if(L==11)
{                   //
p=352711;
}




std::unordered_map<vector<int>, node, decltype(&name_hash)> VTable (p, name_hash);           //
std::unordered_map<vector<int>, node, decltype(&name_hash)> WTable (p, name_hash);           // Setting up the required hash tables
typename std::unordered_map<vector<int>, node, decltype(&name_hash)>::iterator it;                    // Setting up the required iterators
typename std::unordered_map<vector<int>, node, decltype(&name_hash)>::iterator it2;                   //   



for(int i=0;i<(C.size());i++)
{
VTable[C[i]].vectcomp=v[i];
VTable[C[i]].loop=C[i];
}

for(int i=0;i<(C.size());i++)
{
WTable[C[i]].vectcomp=0;
WTable[C[i]].loop=C[i];
}




/*
for(it=VTable.begin();it!=VTable.end();it++)
{
for(int i=0; i<2*L;i++)
{
cout<<VTable[it->first].loop[i];
}
cout<<"	  "<<VTable.bucket(it->first)<<endl;

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

s=it->first;
if (VTable[s].loop[i]==i+1 || VTable[s].loop[i]==-(i+1)){
	if(VTable[s].loop[i]<=0 || VTable[s].loop[i+1]<=0)
		{
		WTable[s].vectcomp=WTable[s].vectcomp+VTable[s].vectcomp;
		WTable[s].loop=VTable[s].loop;

		z=WTable[s].loop;
		z[i]=abs(z[i]);
		z[i+1]=abs(z[i+1]);

	
		WTable[z].vectcomp=WTable[z].vectcomp+x2*loopweightblob*VTable[s].vectcomp;
		//WTable[z].loop=z;
		}		
	else
		{
		WTable[s].vectcomp=WTable[s].vectcomp+VTable[s].vectcomp+x1*loopweight*VTable[s].vectcomp;
		//WTable[it->first].loop=VTable[it->first].loop;
		}

}									//end of loop[i]=i+1 if statement


else
{
WTable[s].vectcomp=WTable[s].vectcomp+VTable[s].vectcomp;
//WTable[it->first].loop=VTable[it->first].loop;

z=WTable[s].loop;

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

z[i]=i+1;
z[i+1]=i;


WTable[z].vectcomp=WTable[z].vectcomp+x1*VTable[s].vectcomp;
//WTable[z].loop=z;						
}								//end of else statement


}								//end of input components loop


for(it=WTable.begin();it!=WTable.end();it++)				//
{									//Sets V=W
s=it->first;								//
									//									//
VTable[s].vectcomp=WTable[s].vectcomp;					//
}									//

for(it=WTable.begin();it!=WTable.end();it++)				//
{									//Sets W=0
s=it->first;								//						
WTable[s].vectcomp=0;							//
}
 								    


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

/*cout << i << endl;
for(it2=VTable.begin();it2!=VTable.end();it2++)
{
for(int i1=0; i1<2*L;i1++)
{
cout<<VTable[it2->first].loop[i1];
}
cout<<VTable[it2->first].vectcomp<<endl;

}
 cin>>r;*/




}								//end of even i loop

for(int i=(2*L-2); i>=0; i=i-2){					//loop over odd Xi (odd when considering first position as position 1, not 0)
for (it=VTable.begin(); it!=VTable.end(); it++){		//loop over components of input vector

s=it->first;

if (VTable[s].loop[i]==i+1 || VTable[s].loop[i]==-(i+1)){
	if(VTable[s].loop[i]<=0 || VTable[s].loop[i+1]<=0)
		{
		WTable[s].vectcomp=WTable[s].vectcomp+x1*VTable[s].vectcomp;
		//WTable[it->first].loop=VTable[it->first].loop;
		
		z=WTable[s].loop;
		if(i!=0){
		z[i]=abs(z[i]);
		z[i+1]=abs(z[i+1]);
		}

		
		WTable[z].vectcomp=WTable[z].vectcomp+loopweightblob*VTable[s].vectcomp;
		//WTable[z].loop=z;
		}		
	else
		{
		WTable[s].vectcomp=WTable[s].vectcomp+x1*VTable[s].vectcomp+loopweight*VTable[s].vectcomp;
		//WTable[it->first].loop=VTable[it->first].loop;
		}

}									//end of loop[i]=i+1 if statement

else
{
WTable[s].vectcomp=WTable[s].vectcomp+x2*VTable[s].vectcomp;
//WTable[it->first].loop=VTable[it->first].loop;

z=WTable[s].loop;

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

z[i]=i+1;
z[i+1]=i;

if(i==0)
z[0]=-z[0];



WTable[z].vectcomp=WTable[z].vectcomp+VTable[s].vectcomp;
//WTable[z].loop=z;
}

						

}								//end of input components loop


/*cout << i << endl;
for(it2=WTable.begin();it2!=WTable.end();it2++)
{
for(int i1=0; i1<2*L;i1++)
{
cout<<WTable[it2->first].loop[i1];
}
cout<<WTable[it2->first].vectcomp<<endl;

}
 cin>>r;
*/


for(it=WTable.begin();it!=WTable.end();it++)				//
{									//
s=it->first;								//
VTable[s].vectcomp=WTable[s].vectcomp;					//Sets V=W
									//
}									//

for(it=WTable.begin();it!=WTable.end();it++)				//
{									//Sets W=0
s=it->first;								//
WTable[s].vectcomp=0;							//
}



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

}								//end of odd i loop


/*cout << WTable.size()<<endl;
cin>>r;
*/
for(int i=0;i<(C.size());i++)
{
w[i]=VTable[C[i]].vectcomp;
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

if(a.size()==22)
{
     for(int i=0; i<21;i++)
{

	if (a[i]>=0)
	total=(22*total+abs(a[i]))%352711;
	else
	total=(22*total+abs(a[i])+21)%352711;

}
                          
int totalmod;
totalmod=total%352711;   						//Change this depending on number of expected states                 
                                            
 return totalmod; 
}				//end of if 22

if(a.size()==20)
{
     for(int i=0; i<19;i++)
{

	if (a[i]>=0)
	total=(20*total+abs(a[i]))%92377;
	else
	total=(20*total+abs(a[i])+19)%92377;

}
                          
int totalmod;
totalmod=total%92377;   						//Change this depending on number of expected states                 
                                            
 return totalmod; 
}				//end of if 20

if(a.size()==18)
{
     for(int i=0; i<17;i++)
{

	if (a[i]>=0)
	total=(18*total+abs(a[i]))%24317;
	else
	total=(18*total+abs(a[i])+17)%24317;

}
                          
int totalmod;
totalmod=total%24317;   						//Change this depending on number of expected states                 
                                            
 return totalmod;  
}				//end of if 18

if(a.size()==16)
{
    for(int i=0; i<15;i++)
{

	if (a[i]>=0)
	total=(16*total+abs(a[i]))%6449;
	else
	total=(16*total+abs(a[i])+15)%6449;

}
                          
int totalmod;
totalmod=total%6449;   						//Change this depending on number of expected states                 
                                            
 return totalmod;  
}				//end of if 16

if(a.size()==14)
{

     for(int i=0; i<13;i++)
{
	
	if (a[i]>=0)
	total=(14*total+abs(a[i]))%1709;
	else
	total=(14*total+abs(a[i])+13)%1709;

}
                          
int totalmod;
totalmod=total%1709;   						//Change this depending on number of expected states                 
                                            
 return totalmod;

}			//end of if 14
                          
}


#endif // CMATRIXA_H
