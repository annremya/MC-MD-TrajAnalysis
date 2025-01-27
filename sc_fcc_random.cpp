#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include "para.h"
#include <math.h>
int main (int argc, char *argv[])
{
    srand(time(0));
    if(argc > 1) {
        ind = atoi(argv[1]);
	start=atoi(argv[2]);
	N = atoi(argv[3]);
	rho = atof(argv[4]);
    }

box = pow ((double(N)/rho), 1./2.);
cout<<"ind  "<<ind<<"start   "<<start<<"N   "<<N<<"rho   "<<rho<<"box   "<<box<<"\n";
initialcondition();
return 0;
}

void initialcondition()
{
    if (start ==0){   //simple cubic
	int index=0;
        int baseunit = int(pow (N, 1./2.))+1;
	double lengthcell = box/baseunit;
	//for (int iz=0; iz<baseunit; iz++) {
	  for(int iy=0; iy<baseunit; iy++) {
	     for(int ix=0; ix<baseunit; ix++){
		colloid[index].x=double(ix)*lengthcell;
		colloid[index].y=double(iy)*lengthcell;
		//colloid[index].z=double(iz)*lengthcell;
		index+=1;
            // }
	}
	}	 
}
	else if(start == 1){                          //fcc crystal
        int index = 0;
        int baseunit;
	baseunit = ceil(pow(N/4, 1.0/3.0));
        double lengthcell = double(box)/double(baseunit);
	cout<<"length"<<lengthcell<<"\n";
//	double facto = 1.0/sqrt(2.0);
        for (int iz=0; iz < baseunit ; iz++) {
            for (int iy=0; iy < baseunit ; iy++) {
                for (int ix = 0; ix < baseunit; ix++) {
                    colloid[index].x = double(ix)*lengthcell;
                    colloid[index].y = double(iy)*lengthcell;
                    colloid[index].z = double(iz)*lengthcell;
		    colloid[index+1].x = colloid[index].x+lengthcell/2.0;
                    colloid[index+1].y = colloid[index].y+lengthcell/2.0;
                    colloid[index+1].z = colloid[index].z;
                    colloid[index+2].x = colloid[index].x+lengthcell/2.0;
                    colloid[index+2].y = colloid[index].y;
                    colloid[index+2].z = colloid[index].z+lengthcell/2.0;
	            colloid[index+3].x = colloid[index].x;
                    colloid[index+3].y = colloid[index].y+lengthcell/2.0;
                    colloid[index+3].z = colloid[index].z+lengthcell/2.0;
                    index+=4;
                }
            }
        }
    }               
        
    else if (start == 2) {                 //create random initial configuration// begin //
    
    colloid[0].x = (randnum()-0.5)*box;
    colloid[0].y = (randnum()-0.5)*box;            
    colloid[0].z = (randnum()-0.5)*box;
    double r1, r2;
    Cordinates rij, t; 
    double hbox = box/2.0;
    for (int i=1;i<N; i++) {
        do {
            t.x = (randnum()-0.5)*box;
            t.y = (randnum()-0.5)*box;
            t.z = (randnum()-0.5)*box;
            for (int j=0; j<i; j++) {
		rij.x = t.x - colloid[j].x;
		rij.y = t.y - colloid[j].y;
		rij.z = t.z - colloid[j].z;
                if (rij.x >  hbox) rij.x -= box;
        	if (rij.x < -hbox) rij.x += box;
        	if (rij.y >  hbox) rij.y -= box;
        	if (rij.y < -hbox) rij.y += box;
        	if (rij.z >  hbox) rij.z -= box;
        	if (rij.z < -hbox) rij.z += box;
                r2 = rij.x*rij.x+rij.y*rij.y+rij.z*rij.z;
                r1 = sqrt(r2);
                if(r1<1.0) break;
            }
        }while (r1<1.0);
        
        colloid[i].x = t.x;
        colloid[i].y = t.y;
        colloid[i].z = t.z;
    }
    }

   	else if (start == 3) {               //read from file// begin //

        char dum[100];
        FILE * fin;
        char IntStr[80];
        sprintf( IntStr, "configin.%d.xyz", ind);
        fin = fopen (IntStr, "r");
        fscanf(fin, "%d", &N);
        fscanf(fin, " ");
        for (int i = 0; i<N; i++) {
            fscanf(fin, "%lf%lf%lf", &colloid[i].x, &colloid[i].y, &colloid[i].z);
        }
}
       
    writeconf();    
return;
}

double randnum()
{
    return double(rand())/double(RAND_MAX);
}

void writeconf()
{
    char IntStr[80];
    ofstream of;
    sprintf( IntStr, "config%d.xyz", ind);
    of.open (IntStr, ofstream::out | ofstream::trunc);
    if (of.is_open())
    {
        of << N<< "\n";
        of << "\n";
	if ((start==0) || (start ==2)){
        for (int i=0; i<N; i++)
        {
            of << 'h'<<"  "<<colloid[i].x << "  " << colloid[i].y << "  "<< colloid[i].z << "\n";
        }    
    } 
	else if (start ==1) {
	 for (int i=0; i<N; i++)
        {
	   if (i%4 ==0) {
            of <<"he"<<"  "<<colloid[i].x << "  " << colloid[i].y << "  "<< colloid[i].z << "\n";
	   }else {of << 'h'<<"  "<<colloid[i].x << "  " << colloid[i].y << "  "<< colloid[i].z << "\n";}
	}
        }

        }
	else {cerr << "unable to open file for config output \n";}
    of.close();
    
    return;
}
