for(int i=(2*L-3); i>0; i=i-2){					//loop over even Xi (even when considering first position as position 1, not 0)
for (it=VTable.begin(); it!=VTable.end(); it++){		//loop over components of input vector


if (VTable[it->first].loop[i]==i+1 || VTable[it->first].loop[i]==-(i+1)){
	if(VTable[it->first].loop[i]<=0 || VTable[it->first].loop[i+1]<=0)
		{
		WTable[it->first].vectcomp=WTable[it->first].vectcomp+VTable[it->first].vectcomp;
		WTable[it->first].loop=VTable[it->first].loop;

		z=WTable[it->first].loop;
		z[i]=abs(z[i]);
		z[i+1]=abs(z[i+1]);

	
		WTable[z].vectcomp=WTable[z].vectcomp+x2*loopweightblob*VTable[it->first].vectcomp;
		//WTable[z].loop=z;
		}		
	else
		{
		WTable[it->first].vectcomp=WTable[it->first].vectcomp+VTable[it->first].vectcomp+x1*loopweight*VTable[it->first].vectcomp;
		//WTable[it->first].loop=VTable[it->first].loop;
		}

}									//end of loop[i]=i+1 if statement


else
{
WTable[it->first].vectcomp=WTable[it->first].vectcomp+VTable[it->first].vectcomp;
//WTable[it->first].loop=VTable[it->first].loop;

z=WTable[it->first].loop;

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


WTable[z].vectcomp=WTable[z].vectcomp+x1*VTable[it->first].vectcomp;
//WTable[z].loop=z;						
}								//end of else statement


}								//end of input components loop


for(it=WTable.begin();it!=WTable.end();it++)				//
{									//Sets V=W
									//									//
VTable[it->first].vectcomp=WTable[it->first].vectcomp;			//
}									//


 								    


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

}								//end of even i loop
