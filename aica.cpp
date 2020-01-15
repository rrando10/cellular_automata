#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <random>
#include <string>
#include <fstream>
using namespace std;

vector<vector<int> > grid;
double Correlation[15];
double MutualInfo[15];
double JointEntropy[15];
double Entropy;

int getDist(int x1,int y1,int x2,int y2)
{
	int xdiff = abs(x1 - x2);
	int ydiff = abs(y1 - y2);

	if(xdiff > 15)
		xdiff = (30 - xdiff);
	if(ydiff > 15)
		ydiff = (30 - ydiff);

	return (xdiff + ydiff);
}

int randToggle()
{
	random_device rd;
	mt19937 gen(rd());
	bernoulli_distribution d(0.3);
	
	if(d(gen))
		return 1;
	else
		return -1;
}

void getEntropy()
{
	double Bsum = 0;
	double posP, negP;
	double P1,P2;

	//add up converted cell values
	for(int i = 0; i < 30; i++)
		for(int j = 0; j < 30; j++)
			Bsum += ((1 + grid[i][j])/2);

	//calculate probabilites of states
	posP = (Bsum / 900);
	negP = (1 - posP);

	//********** avoid 0 ************
	if(posP == 0)
		P1 = 0;
	else
		P1 = (posP * log2(posP));

	if(negP == 0)
		P2 = 0;
	else
		P2 = (negP * log2(posP));
	//******************************
	
	//calculate and set overall entropy
	Entropy = (-1) * (P1 + P2);
}

void getJointE()
{
	int l;
	int i_x,i_y;
	int j_x,j_y;
	double pBSum = 0;
	double nBSum = 0;
	double p1,p2,p3;
	double P1,P2,P3;

	for(l = 0; l < 15; l++){
		
		pBSum = 0;
		nBSum = 0;

		for(i_x = 0; i_x < 30; i_x++)						//iterate through all cells - i 
			for(i_y = 0; i_y < 30; i_y++)
				for(j_x = 0; j_x < 30; j_x++)
					for(j_y = 0; j_y < 30; j_y++){			//iterate through all cells - j
						if((i_y == j_y) && (j_x <= i_x))
							continue;

						if(getDist(i_x,i_y,j_x,j_y) == l){
							pBSum += (((1.0 + grid[i_x][i_y])/2) * ((1 + grid[j_x][j_y])/2));
							nBSum += (((1.0 - grid[i_x][i_y])/2) * ((1 - grid[j_x][j_y])/2));
						}
					}

		if(l == 0)
			JointEntropy[0] = Entropy;
		
		else{
			p1 = ((2.0/(900*(4*l))) * pBSum); //P(+1,+1)
			p2 = ((2.0/(900*(4*l))) * nBSum); //P(-1,-1)
			p3 = (1 - p1 - p2);				//P(+1,-1) or P(-1,+1)

//			cout<<"Joint("<<l<<"): pBSum = "<<pBSum<<endl;
//			cout<<"Joint("<<l<<"): nBSum = "<<nBSum<<endl;
//			cout<<"Joint("<<l<<"): p1 = "<<p1<<endl;
//			cout<<"Joint("<<l<<"): p2 = "<<p2<<endl;
//			cout<<"Joint("<<l<<"): p3 = "<<p3<<endl;
//			//********** avoid 0 ************
			if(p1 == 0)
				P1 = 0;
			else
				P1 = (p1 * log2(p1));

			if(p2 == 0)
				P2 = 0;
			else
				P2 = (p2 * log2(p2));

			if(p3 <= 0)
				P3 = 0;
			else
				P3 = (p3 * log2(p3));
			//*****************************
//			cout<<"Joint("<<l<<"): P1 = "<<P1<<endl;
//			cout<<"Joint("<<l<<"): P2 = "<<P2<<endl;
//			cout<<"Joint("<<l<<"): P3 = "<<P3<<endl<<endl;
			JointEntropy[l] = (-1) * (P1 + P2 + P3);
		}
	}
}

void getMI()
{
	for(int l = 0; l < 15; l++)
		MutualInfo[l] = (2 * Entropy - JointEntropy[l]);
}

void getCorrelation()
{
	int l;
	int i_x,i_y;
	int j_x,j_y;
	double sum1,sum2;

	//loop1: distances 0-14
	for(l = 0; l < 15; l++){
		
		sum1 = 0;
		sum2 = 0;

		//loop2: x val for cell i
		for(i_x = 0; i_x < 30; i_x++){
			for(i_y = 0; i_y < 30; i_y++){
				
				//add value of cell i to sum
				sum2 += grid[i_x][i_y];
				
				//loop4: x val for cell j
				for(j_x = 0; j_x < 30; j_x++){
					for(j_y = 0; j_y < 30; j_y++){
						
						//duplicate condition
						if((i_y == j_y) && (j_x <= i_x))
							continue;
						
						//add product of cell i,j to total
						if(getDist(i_x,i_y,j_x,j_y) == l)
							sum1 += (grid[i_x][i_y] * grid[j_x][j_y]);
					}
				}
			}
		}

		//calculate correlation
		if(l == 0)
			Correlation[0] = abs(1 - pow(((1.0/900)*sum2),2));
		else{
			double nbh = (2.0/(900*(4*l)));
			double P1 = (nbh * sum1);

			double OoN = (1.0/900);
			double P2 = (OoN * sum2);

			Correlation[l] = abs(P1 - pow(P2,2));  
		}
	}
}

void initGrid()
{	
	//declare row vector for push_back
	vector<int> row(30);
	
	//push 30 rows (vectors) onto vector
	for(int i = 0; i < 30; i++){

		//randomly fill row with -1 and 1
		for(int j = 0; j < 30; j++)
			row[j] = randToggle();
	
		//add row to matrix
		grid.push_back(row);
	}
	cout<<endl;
}

void printGrid()
{
	cout<<endl;
	for(int i = 0; i < grid.size(); i++){
		for(int j = 0; j < grid[i].size(); j++){
			if(grid[i][j] == -1)
				cout<<" ";
			else
				cout<<"*";
		}
		cout<<endl;
	}
}

void updateGrid(int R1,int R2, double J1, double J2, int h)
{
	int updated[900];
	int rCell,r_x,r_y;
	int j_x,j_y;
	bool updating;
	bool unstable;
	int near = 0;
	int far = 0;
	int d;
	double p1,p2,state;

	srand(time(NULL));
	unstable = true;
	int n = 0;
	while(unstable){
		n = n+1;

		unstable = false;
		updating = true;
		for(int i = 0; i < 900; i++)
			updated[i] = 0;

		while(updating){

			rCell = rand() % 900;

			//if cell has not been updated
			if(updated[rCell] == 0){
				r_x = rCell / 30;
				r_y = rCell % 30;
				near = 0;
				far = 0;

				for(j_x = 0; j_x < 30; j_x++)
					for(j_y = 0; j_y < 30; j_y++)
						if((r_x != j_x)||(r_y != j_y)){
			
							d = getDist(r_x,r_y,j_x,j_y);

							if(d < R1)
								near += grid[j_x][j_y];
							else if((d >= R1)&&(d < R2))
								far += grid[j_x][j_y];
						}
				
				p1 = (J1 * near);
				p2 = (J2 * far);
				state = (h + p1 + p2);

				if((state < 0)&&(grid[r_x][r_y] != -1)){
					grid[r_x][r_y] = -1;
					unstable = true;
				}else if ((state >= 0)&&(grid[r_x][r_y] != 1)){
					grid[r_x][r_y] = 1;
					unstable = true;
				}

				updated[rCell] = 1;
				updating = false;

				for(int i = 0; i < 900; i++)
					if(updated[i] == 0)
						updating = true;
			}
		}
	}
	cout << n << " iterations completed\n";
}

void outputCSV(int R1, int R2, double J1, double J2, int h)
{
	ofstream ofs;
	string headers[3] = {"correlation-","jointentropy-","mutualinfo-"};
	double *data[3] = {Correlation,JointEntropy,MutualInfo};
	string file; 
	int n = 3;
	for(int i=0; i < 3; i++){ 
		file = headers[i] + to_string(n) + ".csv";
		ofs.open(file,ios::app);
		//ofs<<file<<endl;
		for(int l=0; l<15; l++)
			ofs << data[i][l] <<",";
		ofs<<endl;
		ofs.close();
	}

	ofs.open("manifest.csv",ios::app);
	ofs << J1 << "," << J2 << "," << R1 << "," << R2 << "," << h << "," << Entropy << endl;
	ofs.close();
}
void bumpCSVs()
{
	ofstream bump;
	string headers[3] = {"correlation-","jointentropy-","mutualinfo-"};
	
	string file; 
	int n = 3;
	
	for(int i=0; i < 3; i++){ 
		file = headers[i] + to_string(n) + ".csv";
		bump.open(file,ios::app);
		bump<<endl;
		bump.close();
	}

	bump.open("manifest.csv",ios::app);
	bump<<endl;
	bump.close();

}
void drawPic(int i)
{
	ofstream pic;
	string fname = "pic-" + to_string(i) +".pgm";
	pic.open(fname);

	pic << "P2" << endl;
	pic << "30 30" << endl;
	pic << "255" << endl;

	for(int x=0; x<30; x++){
		for(int y=0; y<30; y++){
			if(grid[x][y] == 1)
				pic << "0 ";
			else
				pic << "255 ";
		}
		pic << endl;
	}
	pic.close();
}
int main(int argc, char* argv[])
{
	int R1,R2,h;
	double J1,J2;

	J1 = 1;
	J2 = -0.1;
	R1 = 12;
	
	R2 = 14;
	h = 25;

	int n = 8;
	int record = 0;
	int draw = 0;

	for(int i=0; i<4; i++){
		initGrid();
		updateGrid(R1,R2,J1,J2,h);
		printGrid();

		getCorrelation();
		getEntropy();
		getJointE();
		getMI();
		if(record == 1)
			outputCSV(R1,R2,J1,J2,h);

		grid.clear();
	}

	if(record == 1)
		bumpCSVs();
	
	if(draw == 1)
		drawPic(n);
	
	return 1;
}
