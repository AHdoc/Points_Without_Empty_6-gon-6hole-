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

int nbvar,nbclauses,nbliterals;

/* Variables */
int idx[MAXN+1][MAXN+1][MAXN+1]; // idx[p][q][r]: The triangle pqr is oriented
int idx2[MAXN+1][MAXN+1][MAXN+1][MAXN+1]; 
int idx3[MAXN+1][MAXN+1][MAXN+1];
int idx4[MAXN+1][MAXN+1][MAXN+1][MAXN+1]; 

unordered_set<vector<int>,VectorHash> clauses;

int new_var(){ // Set a new variable
	++nbvar;
	return nbvar;
}

void new_clauses(vector<int> d){ // Give a new clause
	sort(d.begin(),d.end());
	if(clauses.find(d)==clauses.end()){
		++nbclauses;
		nbliterals+=d.size();
		clauses.insert(d);
	}
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

void new_known(int i,int j,int k){
	new_known(idx[i][j][k]);
}

void cc_system(int n,vector<int> hull,int m){
	for(int i=1;i<=n;i++)
		for(int j=i+1;j<=n;j++)
			for(int k=j+1;k<=n;k++){
				int x=new_var();
				idx[i][j][k]=idx[j][k][i]=idx[k][i][j]=x;
				idx[i][k][j]=idx[j][i][k]=idx[k][j][i]=-x;
			}
/* known information */
	for(int i=1,i2,j=0;j<hull.size();i=i2,j++){ 
		if(j>=1 && hull[j-1]==3 && hull[j]>=2){ // Symmetry Breaking for Triangles
			int p=i-3,q=i-2,x=i,y=i+1;
			new_known(p,y,x);
			new_known(q,x,y);
		} 
		if(j>=1 && hull[j-1]==4){ // Symmetry Breaking for Quadrilaterals
			int p=i-4,q=i-3,r=i-2,s=i-1;
			new_known(q,s,i);
			new_known(r,p,i);
		}
		if(j>=1 && hull[j-1]==5){ // Symmetry Breaking for Pentagons
			int p=i-5,q=i-4,r=i-3,s=i-2,t=i-1;
			int x1=new_var(),x6=new_var(),x11=new_var();
			new_known(q,s,i);
			new_known(s,p,i);
			new_known(r,t,i);
			new_clauses({idx[q][t][i],idx[p][r][i]}); 
		} 
		if(j>=1 && hull[j-1]==6){ // Symmetry Breaking for Hexagons
			int p=i-6,q=i-5,r=i-4,s=i-3,t=i-2,u=i-1;
			int x1=new_var(),x6=new_var(),x11=new_var();
			new_known(q,t,i);
			new_known(s,p,i);
		}
		
		i2=i+hull[j];
		vector<int> layerpts;
		for(int k=i;k<i2;k++)
			layerpts.push_back(k);
		int tot=layerpts.size();
		if(tot>=3){
			for(int u=0;u<tot;u++)
				for(int v=(u+1)%tot;(v+1)%tot!=u;v=(v+1)%tot)
					for(int w=(v+1)%tot;w!=u;w=(w+1)%tot){
						int uu=layerpts[u];
						int vv=layerpts[v];
						int ww=layerpts[w];
						new_known(uu,vv,ww);
					}
			for(int u=0;u<tot;u++){
				int v=(u+1)%tot;
				int p=layerpts[u],q=layerpts[v];
				for(int r=i2;r<=n-m;r++)
					new_known(p,q,r);
			}
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
/* ----------------- */
	// n-m+1 ~ n are outside
	for(int p=n-m+1;p<=n;p++){
		vector<int> tmp;
		for(int i=1;i<=hull[0];i++){
			int j=(i==hull[0]?1:i+1);
			tmp.push_back(idx[j][i][p]);
			new_clauses(tmp);
		}
	}
}

// idx2[p][q][r][s]: s is inside the triangle pqr
// idx3[p][q][r]   : the triangle pqr has no point inside
void var_pt_inside_triangle(int n,vector<int> hull){
	for(int p=1;p<=n;p++)
		for(int q=p+1;q<=n;q++)
			for(int r=p+1;r<=n;r++){
				vector<int> empty_cond;
				for(int s=1;s<=n;s++){
					set<int> pts={p,q,r,s};
					if(pts.size()!=4) continue;
					
					int x=new_var(); // x: s is inside triangle pqr
					
					idx2[p][q][r][s]=idx2[q][r][p][s]=idx2[r][p][q][s]=x;
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
void var_4pt_region_empty(int n,vector<int> hull){
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

void mk_no6hole_givenhullstructure_outside(vector<int> hull,int m){
	int n=m;
	for(int x:hull) n+=x;
	
	nbvar=0;
	nbclauses=0;
	nbliterals=0;
	clauses.clear();
	
	cc_system(n,hull,m);
	var_pt_inside_triangle(n,hull);
	var_4pt_region_empty(n,hull);
	
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
	
	string s="no-6hole-";
	for(int x:hull)
		s.push_back(char('0'+x));
	s+="-outside-"+to_string(m)+".sat";
	ofstream fout(s);
	
	fout<<"c\n";
	fout<<"p cnf "<<nbvar<<" "<<clauses.size()<<"\n"; 
	for(auto elm:clauses){
		for(int x:elm)
			fout<<x<<" ";
		fout<<"0\n";
	}
	fout.close();
	
	cerr<<n<<" ("<<nbvar<<","<<nbclauses<<","<<nbliterals<<")\n";
}

int main(){
	mk_no6hole_givenhullstructure_outside({7,1,0},18);
	mk_no6hole_givenhullstructure_outside({7,1,0},19);
	mk_no6hole_givenhullstructure_outside({7,1,0},20);
} 



