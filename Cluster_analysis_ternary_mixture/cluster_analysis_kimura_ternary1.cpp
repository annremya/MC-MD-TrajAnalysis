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


const int N=2010;
const int Np=1;

class coordinates
{
	public: double x,y,z;
};

int	a, binindex;
int     nA, nB, nC, npart, dum, frame_count, proconfig;
int	moltp[N];
char	filename[100];
double box, hbox, delg, distance;
coordinates r[N], e[N], dr;
//--------------------for cluster analysis-----------------------------------------------
float dcrit = 2.1;
int cluster_nr; // number of clusters (by default it is 1)
int biggest_cluster, biggest_cluster_index;
int cluster_status[N], cluster_index[N];	//array size shud be => N to avoid memory issues
int particle_in_cluster[N];
int number_cluster[N]; //number of clusters with a particular number of particles

void cluster_analysis();
void cluster_search(int j);
void writexyzfile_forcluster(); //write xyz file for individual clusters
void clust_size_dist();	//sample the cluster size distribution

void cluster_analysis()
{
	//cluster status if 0, i is not asssigned to a cluster, if 1, i is assigned to a cluster
	for(int i=0; i<nA+nB; i++){
		cluster_status[i]=0;
	}
	cluster_nr = 1;
	//cluster_search
	for(int i=0; i<(nA+nB); i++){
		if(cluster_status[i]==0){
			cluster_status[i] = 1;
			cluster_search(i);
			cluster_index[i] = cluster_nr;	//assign i to current cluster index
			//-----------
			int check=0;
			for(int a=0; a<(nA+nB); a++){
				if((i!=a)&&(cluster_index[a]==cluster_nr) && (moltp[i]!=moltp[a])){
					check = 1;
					break;
				} 
			}
			if(check == 0){
				for(int a=0; a<(nA+nB); a++){	
					if(cluster_index[a]==cluster_nr){
						cluster_index[a] = 0;
						
					}
				}
				cluster_nr--;
			}
			//------------------
			cluster_nr++;
		}
	}
	//find number of paricles in each cluster
	//int particle_in_cluster[cluster_nr+1];  //made global
	for(int i=1; i<=cluster_nr; i++){
		particle_in_cluster[i]=0;
		
	}
	for(int i=0; i<(nA+nB); i++)
	{
		particle_in_cluster[cluster_index[i]]++;
		
	}
	//update the number of clusters of a particular size
	for(int i=1; i<=cluster_nr; i++){
		number_cluster[particle_in_cluster[i]]++;
	}
	//find the biggest_cluster
	biggest_cluster = 1;
	for(int i=1; i<=cluster_nr; i++){
		if(particle_in_cluster[i]>biggest_cluster){
			biggest_cluster       = particle_in_cluster[i];
			biggest_cluster_index = i;
		}
	}
	char filename[100];
	ofstream bigcl;
	sprintf(filename,  "biggest_cluster.dat");
	bigcl.open(filename, ofstream::out | ofstream::app);
	if(bigcl.is_open()){
		bigcl << frame_count << "  " << biggest_cluster_index << "  " << biggest_cluster << endl;
	} else {cerr << "unable to open bigcl file \n";}
	bigcl.close();
}

void writexyzfile_forcluster()
{	
	int clust_N=0;
	//write_xyz of clusters
	ofstream xyzclust;
	sprintf(filename, "mkdir cluster");
	system(filename);
	for(int i=1; i<=cluster_nr; i++){
		if(particle_in_cluster[i]>1) clust_N+= particle_in_cluster[i];
		
		if(particle_in_cluster[i]>5){
			sprintf(filename,"cluster/cluster%d.xyz",i);
			xyzclust.open(filename, ofstream::out | ofstream::trunc);
			if(xyzclust.is_open()){
				xyzclust << particle_in_cluster[i]*2 << "\n";
				xyzclust << "Number of paricles = " << particle_in_cluster[i] << " and box = " << box << "\n";
				for(int j=0; j<(nA+nB); j++){
					if(cluster_index[j]==i){
						//printcubexyz(j,xyzclust);
						if(moltp[j] == 0)
						{
							xyzclust  << "He1"			<< "  "
								 << r[j].x			<< "  "
								 << r[j].y			<< "  "
								 << r[j].z			<< "\n";
							xyzclust  << "He2"			<< "  "
								 << r[j].x + e[j].x/20.	<< "  "
								 << r[j].y + e[j].y/20.	<< "  "
								 << r[j].z + e[j].z/20.	<< "\n";
						 }
						 else if(moltp[j] == 1)
						 {
							xyzclust  << "H1"			<< "  "
								 << r[j].x			<< "  "
								 << r[j].y			<< "  "
								 << r[j].z			<< "\n";
							xyzclust  << "H2"			<< "  "
								 << r[j].x + e[j].x/20.	<< "  "
								 << r[j].y + e[j].y/20.	<< "  "
								 << r[j].z + e[j].z/20.	<< "\n";
						 }

					}
				}
			} else {cerr << "Unable to open xyzclust file number " << i << endl;}
			xyzclust.close();
		}
	}
	sprintf(filename,"config.xyz");
	xyzclust.open(filename, ofstream::out | ofstream::trunc);
	if(xyzclust.is_open()){
		xyzclust << clust_N*2 << "\n";
		xyzclust << "Box = " << box << "\n";
		for(int j=0; j<(nA+nB); j++){
			//if(cluster_index[j]==i){
			if((particle_in_cluster[cluster_index[j]]>1)&&(cluster_index[j]>0)){
				//printcubexyz(j,xyzclust);
				if(moltp[j] == 0)
				{
					xyzclust  << "He1"			<< "  "
						 << r[j].x			<< "  "
						 << r[j].y			<< "  "
						 << r[j].z			<< "\n";
					xyzclust  << "He2"			<< "  "
						 << r[j].x + e[j].x/20.	<< "  "
						 << r[j].y + e[j].y/20.	<< "  "
						 << r[j].z + e[j].z/20.	<< "\n";
				 }
				 else if(moltp[j] == 1)
				 {
					xyzclust  << "H1"			<< "  "
						 << r[j].x			<< "  "
						 << r[j].y			<< "  "
						 << r[j].z			<< "\n";
					xyzclust  << "H1"			<< "  "
						 << r[j].x + e[j].x/20.	<< "  "
						 << r[j].y + e[j].y/20.	<< "  "
						 << r[j].z + e[j].z/20.	<< "\n";
				 }

			}
		}
	} else {cerr << "Unable to open xyzclust file number " << endl;}
	xyzclust.close();
}

void cluster_search(int j)
{
	//cout << "cluster search for " << j << endl;
	double hbox, distance;
	coordinates dr;
	hbox = box/2.;
	for(int i=0; i<(nA+nB); i++){
		vectsub(r[i],r[j],dr);
		if (dr.x >  hbox) dr.x -= box;
		if (dr.x < -hbox) dr.x += box;
		if (dr.y >  hbox) dr.y -= box;
		if (dr.y < -hbox) dr.y += box;
		if (dr.z >  hbox) dr.z -= box;
		if (dr.z < -hbox) dr.z += box;
		distance = vectinnerpdct(dr,dr);
		distance = sqrt(distance);
		
		if((i!=j) && (distance<=dcrit) && (cluster_status[i]==0)){
			//cout << i << "  " << j << "  " << distance << endl;
			//cout << r1 << endl;
			cluster_status[i] = 1;
			cluster_index[i]  = cluster_nr;
			cluster_search(i);
		}
	}
}
//------------------
void clust_size_dist()
{
	double cluster_size, cluster_size_sum=0., cluster_sum=0., cluster_size_sum2=0.;
	ofstream clust_size_file;
	char file[100];
	sprintf(file,"size_dist_prod.dat");
	clust_size_file.open(file, ofstream::out | ofstream::app);
	if(clust_size_file.is_open()){
		for (int i=1; i<=(nA+nB); i++) {
			cluster_size = number_cluster[i]/(double(frame_count)*box*box*box);
			clust_size_file << i << "\t" << cluster_size << endl;	// n vs <n>
			
			cluster_size_sum += (i*cluster_size);
			cluster_sum      += cluster_size;
			cluster_size_sum2+= (i*i*cluster_size);
		}
		double averagesize = cluster_size_sum / cluster_sum;
		double meansquare  = cluster_size_sum2 / cluster_sum;
		double stdsize     = sqrt(meansquare - (averagesize*averagesize))/sqrt(frame_count);
		
		clust_size_file << "#Average Size"      << " = " << averagesize	<< "  "
				<< "Standard deviation" << " = " << stdsize 		<< endl;
	} else {cerr << "Can't open clust_size_dist file \n";}
	clust_size_file.close();
}
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

	
int main(int argc, char *argv[])
{	
	proconfig = 100;
	//proconfig  = atoi(argv[1]);
	//nA	   = atoi(argv[2]);
	//nB	   = atoi(argv[3]);
	
	nA = 50; nB= 500; nC = 1000;
	//nA = 150; nB= 500; nC = 1000;
	//nA = 250; nB= 500; nC = 1000;
	//nA = 500; nB= 500; nC = 1000;
	npart = nA+nB+nC;

//	string ls = GetStdoutFromCommand("find ../eq/ -name '*.dat' -type f | wc -l"); 
//	std::string::size_type sz;
//	int max_steps  = std::stoi (ls, nullptr, 0);
//	max_steps -= 1;
//	cout << "max_steps = " << max_steps << endl;
//		
//	cout << "Enter the value of first frame and last frame: \n";
	int first_f, last_f, step_f;
	//cin >> first_f;
	//cin >> step_f;
	//cin >> last_f;
	
	first_f = 0; last_f = proconfig;
	step_f = 1;
	
	//initialise
	frame_count = 0;  
	//sampling
	
//	for (int k=0; k<proconfig; k++)
	for (int k=first_f; k<last_f; k=k+step_f)
	{
		if (k%10 == 0) cout << k << endl;
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
				if ((k==first_f) && (a==0)) 	//sample only last 20 frames
				{
					for(int i=0; i<N; i++){
						number_cluster[i]=0;	//initialise no. of clusters of a particular size
					}
				}
				a++;
				//if(a >= (frames-20)){
				cluster_analysis();
				frame_count++;
				//}
				continue;
			}
		}
		else {cerr << "Cannot open config_" << k << ".dat file for input.\n";}
	}
	cout << "frame_count is " << frame_count << endl;
	clust_size_dist();
	writexyzfile_forcluster();
	
	sprintf(filename, "gnuplot clustplot");
	system(filename);

return 0;
}







