#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string>
#include<fstream>
#include<iostream>

#define frames 100
#define vectsub(a,b,c) c.x = b.x-a.x; c.y = b.y-a.y; c.z = b.z-a.z;
#define vectinnerpdct(a,b) (a.x*b.x)+(a.y*b.y)+(a.z*b.z);

using namespace std;

//------------------
string GetStdoutFromCommand(string cmd);

string GetStdoutFromCommand(string cmd) 
{

	string data;
	FILE * stream;
	const int max_buffer = 256;
	char buffer[max_buffer];
	cmd.append(" 2>&1");

	stream = popen(cmd.c_str(), "r");
	if (stream) {
		while (!feof(stream))
		if (fgets(buffer, max_buffer, stream) != NULL) data.append(buffer);
		pclose(stream);
	}
	return data;
}
//------------------

const int N=2000;
const int Np=1;

class coordinates
{
	public: double x,y,z;
};

int main(int argc, char *argv[])
{	
	int	a, binindex;
	int     nA, nB, nC, npart, dum, ngr, proconfig;
	int	moltp[N];
	char	filename[100];
	double box, hbox, delg, distance;
	coordinates r[N], e[N], dr;

	int	maxbin=1000;
	double hist[maxbin];
	double grpp[maxbin], grww[maxbin], groo[maxbin];
	double grpw[maxbin], grpo[maxbin], grwo[maxbin];


	//proconfig  = atoi(argv[1]);
	//nA	   = atoi(argv[2]);
	//nB	   = atoi(argv[3]);
	
	nA = 50; nB= 500; nC = 1000;
	//nA = 150; nB= 500; nC = 1000;
	//nA = 250; nB= 500; nC = 1000;
	//nA = 500; nB= 500; nC = 1000;
	npart = nA+nB+nC;

	//string ls = GetStdoutFromCommand("find ../eq/ -name '*.dat' -type f | wc -l");
	//cout << ls << endl; 
	//std::string::size_type sz;
	//int max_steps  = std::stoi (ls, nullptr, 0);
	//max_steps -= 1;
	//cout << "max_steps = " << max_steps << endl;
		
	//cout << "Enter the value of first frame and last frame: \n";
	int first_f, last_f, step_f;
	//cin >> first_f;
	//cin >> step_f;
	//cin >> last_f;
	proconfig = 100;
	step_f = 1;
	first_f = 0; last_f = proconfig;
	
	//initialise
	ngr = 0;
	for (int i=0; i<maxbin; i++) 
	{
	    hist[i] = 0;
	    grpp[i] = 0; grww[i] = 0; groo[i] = 0;
	    grpw[i] = 0; grpo[i] = 0; grwo[i] = 0;	    
	    
	}   
	//sampling
	
//	for (int k=0; k<proconfig; k++)
	for (int k=first_f; k<last_f; k=k+step_f)
	{
		if(k%10 == 0) cout << k << endl;
		ifstream config;
		sprintf(filename,"../pro/config_%d.dat",k);
		//sprintf(filename,"../eq/config_%d.dat",k);		
		config.open(filename, ifstream::in);
		//cout << k << " " << filename << endl;
		
		
//		char command[100];
//		sprintf(command, "wc -l ../eq/config_%d.dat", k);
//		ls = GetStdoutFromCommand(command);

//		int max_lines  = std::stoi (ls, nullptr, 0);
//		int frames = max_lines/(npart+1);
//		cout << "max_lines = " << max_lines << endl;
//		cout << "frames    = " << frames << endl;
		
		
		//if(last_f == 0) {cout << "Enter number of frames to be processed \n";}
				
		if(config.is_open())
		{
			a=0;
			while(a<frames)
			{
				config >> dum >> dum >> dum >> box;
				for(int j=0;j<npart;j++)
				{
					config  >> r[j].x >> r[j].y >> r[j].z
						>> e[j].x >> e[j].y >> e[j].z
						>> moltp[j];
				}
				if ((k==first_f) && (a==0)) 
				{
					delg = box/(2*maxbin);
					cout << "delg is " << delg << endl;
					cout << "box is  " << box  << endl;
				}
				hbox = box/2.;
				for (int i=0; i<npart-1; i++) 
				{
				    for (int j=i+1; j<npart; j++) 
				    {
					vectsub(r[i],r[j],dr);
					if (dr.x >  hbox) dr.x -= box;
					if (dr.x < -hbox) dr.x += box;
					if (dr.y >  hbox) dr.y -= box;
					if (dr.y < -hbox) dr.y += box;
					if (dr.z >  hbox) dr.z -= box;
					if (dr.z < -hbox) dr.z += box;
					distance = vectinnerpdct(dr,dr);
					distance = sqrt(distance);
					if (distance < hbox) 
					{
						binindex = int(distance/delg)+1;
						if(binindex <= maxbin) 
						{
							hist[binindex] += 2;
							if 	((moltp[i] == 0) && (moltp[j] == 0)) grpp[binindex] +=2;
							else if ((moltp[i] == 1) && (moltp[j] == 1)) grww[binindex] +=2;
							else if ((moltp[i] == 2) && (moltp[j] == 2)) groo[binindex] +=2;
							else if ((moltp[i] == 0) && (moltp[j] == 1)) grpw[binindex] +=1;
							else if ((moltp[i] == 0) && (moltp[j] == 2)) grpo[binindex] +=1;
							else if ((moltp[i] == 1) && (moltp[j] == 2)) grwo[binindex] +=1;
							
						}
					}
				    }
				}
				a++;
				ngr++;
				continue;
			}
		}
		else {cerr << "Cannot open config_" << k << ".dat file for input.\n";}
	}
	cout << "ngr is " << ngr << endl;
	
	//calculate rdf
	double rlow, rup, vb, nid, nidA, nidB, nidC, dist;
	double V     = box*box*box;
	double constant = 3.141519;
	ofstream g_out;
	sprintf(filename,"rdf.dat");
	g_out.open(filename, ofstream::out | ofstream::trunc); 
	if(g_out.is_open())
	{
		for(int h=1;h<maxbin;h++)
		{
			//cout << h << "  " << maxbin << endl;
			dist   = delg*(double(h)+0.5);
			rlow   = double(h)*delg;
			rup    = rlow + delg;
			vb     = pow(rup,3.) - pow(rlow,3.); //volume of bin between rlow and rup
			
			nid    = (4./3.)*constant*vb*double(npart)/V;
			nidA   = (4./3.)*constant*vb*double(nA)/V;
			nidB   = (4./3.)*constant*vb*double(nB)/V;
			nidC   = (4./3.)*constant*vb*double(nC)/V;
			
			hist[h]= double(hist[h])/double(ngr)/double(npart)/nid;
			
			grpp[h]= double(grpp[h])/double(ngr)/double(nA)/nidA;
			grww[h]= double(grww[h])/double(ngr)/double(nB)/nidB;
			groo[h]= double(groo[h])/double(ngr)/double(nC)/nidC;
			
			grpw[h]= double(grpw[h])/double(ngr)/double(nA)/nidB;
			grpo[h]= double(grpo[h])/double(ngr)/double(nA)/nidC;
			grwo[h]= double(grwo[h])/double(ngr)/double(nB)/nidC;
			
			g_out << dist << "\t"   << hist[h] << "  "
						<< grpp[h] << "  " << grww[h] << "  " << groo[h] << "  "
						<< grpw[h] << "  " << grpo[h] << "  " << grwo[h] << "\n";
			
		}
	} else { cerr << "Cannot open rdf file\n";}
	g_out.close();
	
	sprintf(filename, "gnuplot rdfplot");
	system(filename);

return 0;
}





