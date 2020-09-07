/*
DIMACS format:
	http://www.satcompetition.org/2009/format-benchmarks2009.html
CC system:
	https://en.wikipedia.org/wiki/CC_system

At least 30 points are needed:
	there exists a set of 29 points in general position with no empty convex hexagon.
*/

#include<iostream>
#include<fstream>
#include<cstdio>
#include<cstring>
#include<vector>
#include<set>
#include<unordered_set>
#include<algorithm>

#include<assert.h> 

using namespace std;

struct VectorHash {
    size_t operator()(const vector<int>& v) const {
        hash<int> hasher;
        size_t seed = 0;
        for (int i : v) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
    }
};

const int MAXN=30;

ofstream fout;

int nbvar;

/* Variables */
int idx[MAXN+1][MAXN+1][MAXN+1]; // idx[p][q][r]: The triangle pqr is oriented
int idx3[MAXN+1][MAXN+1][MAXN+1];
int idx4[MAXN+1][MAXN+1][MAXN+1][MAXN+1]; 

unordered_set<vector<int>,VectorHash> clauses;

int new_var(){ // Set a new variable
	++nbvar;
	return nbvar;
}

void new_clauses(vector<int> d){ // Give a new clause
	sort(d.begin(),d.end());
	clauses.insert(d);
	
	//for(int x:d)
		//fout<<x<<" ";
	//fout<<"0\n";
}

void define_var_and(int nv,vector<int> cond){ // Define the new variable nv as (x1 & x2 & ... & xk) with cond={x1,...,xk}
	vector<int> tmp;
	for(int x:cond) tmp.push_back(-x);
	tmp.push_back(nv);
	new_clauses(tmp);
	
	for(int x:cond)	new_clauses({-nv,x});
}

void new_known(int x){ // Claim that the variable x is known where x may be negative
	new_clauses({x});
}

void cc_system(int n,int m,vector<int> partition){
	for(int i=1;i<=n;i++)
		for(int j=i+1;j<=n;j++)
			for(int k=j+1;k<=n;k++){
				int x=new_var();
				idx[i][j][k]=idx[j][k][i]=idx[k][i][j]=x;
				idx[i][k][j]=idx[j][i][k]=idx[k][j][i]=-x;
			}
/* known information */
	vector<int> layerpts;
	for(int k=1;k<=7;k++) layerpts.push_back(k);
	for(int u=0;u<7;u++)
		for(int v=(u+1)%7;(v+1)%7!=u;v=(v+1)%7)
			for(int w=(v+1)%7;w!=u;w=(w+1)%7){
				int uu=layerpts[u];
				int vv=layerpts[v];
				int ww=layerpts[w];
				new_known(idx[uu][vv][ww]);
			}
	for(int u=0;u<7;u++){
		int uu=layerpts[u];
		for(int k=1;k<=3;k++){
			int v=(u+k)%7;
			int vv=layerpts[v];
			new_known(idx[uu][vv][8]);
		}
		for(int k=4;k<=6;k++){
			int v=(u+k)%7;
			int vv=layerpts[v];
			new_known(-idx[uu][vv][8]);
		}
	}
	// n-m+1 ~ n are outside
	for(int p=9;p<=n;p++){
		vector<int> tmp;
		for(int i=1;i<=7;i++){
			int j=(i==7?1:i+1);
			tmp.push_back(idx[j][i][p]);
		}
		new_clauses(tmp);
	}
	for(int i=0,j,u=9,v;i<partition.size();i++,u=v){
		v=u+partition[i];
		j=(i+1)%7;
		for(int w=u;w<v;w++){
			new_known(idx[8][layerpts[i]][w]);
			new_known(-idx[8][layerpts[j]][w]);
		} 
	}
/* ----------------- */
	// Interiority: If tqr and ptr and pqt, then pqr.
	for(int p=1;p<=n;p++)
		for(int q=p+1;q<=n;q++)
			for(int r=p+1;r<=n;r++)
				for(int t=1;t<=n;t++){
					set<int> pts={p,q,r,t};
					if(pts.size()!=4) continue;
					new_clauses({-idx[t][q][r],-idx[p][t][r],-idx[p][q][t],idx[p][q][r]});
				}
	// Transitivity: If tsp and tsq and tsr, and tpq and tqr, then tpr.
	for(int p=1;p<=n;p++)
		for(int q=1;q<=n;q++)
			for(int r=1;r<=n;r++)
				for(int s=1;s<=n;s++)
					for(int t=1;t<=n;t++){
						set<int> pts={p,q,r,s,t};
						if(pts.size()!=5) continue;
						new_clauses({-idx[t][s][p],-idx[t][s][q],-idx[t][s][r],-idx[t][p][q],-idx[t][q][r],idx[t][p][r]});
					}
}

// idx3[p][q][r]   : the triangle pqr has no point inside
void var_pt_inside_triangle(int n){
	for(int p=1;p<=n;p++)
		for(int q=p+1;q<=n;q++)
			for(int r=p+1;r<=n;r++){
				vector<int> empty_cond;
				for(int s=1;s<=n;s++){
					set<int> pts={p,q,r,s};
					if(pts.size()!=4) continue;
					
					int x=new_var(); // x: s is inside triangle pqr
					define_var_and(x,{idx[p][q][s],idx[q][r][s],idx[r][p][s]});
					empty_cond.push_back(-x);
				}
				int x_empty=new_var();
				idx3[p][q][r]=idx3[q][r][p]=idx3[r][p][q]=x_empty;
				define_var_and(x_empty,empty_cond); 
			}
}

// idx4[p][q][r][s]: A special 4-point region is empty
/*   -------s------r
     xxxxxxx|
	 xxxxxxx|
	 xxxxxxx|
	 -------p------q 
*/ 
void var_4pt_region_empty(int n){
	for(int p=1;p<=n;p++)
		for(int q=1;q<=n;q++)
			for(int r=1;r<=n;r++)
				for(int s=1;s<=n;s++){
					set<int> pts={p,q,r,s};
					if(pts.size()!=4) continue;
					
					vector<int> empty_cond;
					for(int t=1;t<=n;t++){
						if(pts.find(t)!=pts.end()) continue;
						
						int x=new_var(); 
						define_var_and(x,{idx[p][q][t],idx[r][s][t],idx[p][s][t]});
						empty_cond.push_back(-x);
					}
					int y=new_var();
					idx4[p][q][r][s]=y;
					define_var_and(y,empty_cond);
				}
}

void mk_no6hole_710_outside_7parts(int m,vector<int> partition){
	int n=7+1+m;
	
	string s="no-6hole-710-outside-"+to_string(m)+"("+to_string(partition[0]);
	for(int i=1;i<partition.size();i++) s=s+","+to_string(partition[i]);
	s=s+").sat";
	fout.open(s);
	fout<<"c\n";
	//fout<<"p cnf 123 0\n"; 

	nbvar=0;
	clauses.clear();
	
	cc_system(n,m,partition);
	var_pt_inside_triangle(n);
	var_4pt_region_empty(n);
	
	// Restriction: No 6-hole.
	for(int p=1;p<=n;p++)
		for(int q=p+1;q<=n;q++)
			for(int r=p+1;r<=n;r++)
				for(int s=p+1;s<=n;s++){
					set<int> pts={p,q,r,s};
					if(pts.size()!=4) continue;
					
					int empty_r1=idx3[p][q][r],empty_r2=idx3[p][r][s];
					int left_empty=idx4[p][q][r][s],right_empty=idx4[r][s][p][q];
					new_clauses({-idx[p][q][r],-idx[q][r][s],-idx[r][s][p],-idx[s][p][q],-empty_r1,-empty_r2,left_empty,right_empty});
				}
	
	fout<<"p cnf "<<nbvar<<" "<<clauses.size()<<"\n"; 
	for(auto elm:clauses){
		for(int x:elm)
			fout<<x<<" ";
		fout<<"0\n";
	}
	fout.close();
	
	cerr<<n<<" ("<<nbvar<<")\n";
}

int main(){
	mk_no6hole_710_outside_7parts(14,{2,2,2,2,2,2,2});
	mk_no6hole_710_outside_7parts(14,{1,2,2,2,2,2,3});
	mk_no6hole_710_outside_7parts(14,{1,1,2,2,2,3,3});
} 



