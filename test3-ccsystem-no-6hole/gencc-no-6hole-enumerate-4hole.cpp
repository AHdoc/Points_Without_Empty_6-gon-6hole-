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

const int MAXN=50;

int nbhullstructure;
int nbvar,nbclauses,nbliterals;

int lvl[MAXN+1],known[MAXN+1][MAXN+1][MAXN+1],known2[MAXN+1][MAXN+1][MAXN+1][MAXN+1];
int idx[MAXN+1][MAXN+1][MAXN+1];
int idx2[MAXN+1][MAXN+1][MAXN+1][MAXN+1];

unordered_set<vector<int>,VectorHash> clauses;

void new_clauses(vector<int> d){
	sort(d.begin(),d.end());
	if(clauses.find(d)==clauses.end()){
		++nbclauses;
		nbliterals+=d.size();
		clauses.insert(d);
	}
}

void define_var_and(int nv,vector<int> cond){
	vector<int> tmp;
	for(int x:cond) tmp.push_back(-x);
	tmp.push_back(nv);
	new_clauses(tmp);
	
	for(int x:cond)	new_clauses({-nv,x});
}

int get(int i,int j,int k){
	assert(i!=j && i!=k && j!=k);
	if(i<j && i<k){
		if(j<k) return idx[i][j][k];
		else return -idx[i][k][j];
	}else
		return get(j,k,i);
}

void cc_system(int n){
	for(int i=1;i<=n;i++)
		for(int j=i+1;j<=n;j++)
			for(int k=j+1;k<=n;k++){
				++nbvar;
				idx[i][j][k]=idx[j][k][i]=idx[k][i][j]=nbvar;
				idx[i][k][j]=idx[j][i][k]=idx[k][j][i]=-nbvar;
			}
	// Interiority: If tqr and ptr and pqt, then pqr.
	for(int p=1;p<=n;p++)
		for(int q=p+1;q<=n;q++)
			for(int r=p+1;r<=n;r++)
				for(int t=1;t<=n;t++){
					set<int> pts={p,q,r,t};
					if(pts.size()!=4) continue;
					new_clauses({-get(t,q,r),-get(p,t,r),-get(p,q,t),get(p,q,r)});
				}
	// Transitivity: If tsp and tsq and tsr, and tpq and tqr, then tpr.
	for(int p=1;p<=n;p++)
		for(int q=1;q<=n;q++)
			for(int r=1;r<=n;r++)
				for(int s=1;s<=n;s++)
					for(int t=1;t<=n;t++){
						set<int> pts={p,q,r,s,t};
						if(pts.size()!=5) continue;
						new_clauses({-get(t,s,p),-get(t,s,q),-get(t,s,r),-get(t,p,q),-get(t,q,r),get(t,p,r)});
					}
}

void var_pt_inside_triangle(int n){
	for(int p=1;p<=n;p++)
		for(int q=1;q<=n;q++)
			for(int r=1;r<=n;r++)
				for(int s=1;s<=n;s++)
					idx2[p][q][r][s]=-1;
	for(int p=1;p<=n;p++)
		for(int q=1;q<=n;q++)
			for(int r=1;r<=n;r++)
				for(int s=1;s<=n;s++){
					set<int> pts={p,q,r,s};
					if(pts.size()!=4) continue;
					if(idx2[p][q][r][s]!=-1) continue;
					idx2[p][q][r][s]=idx2[q][r][p][s]=idx2[r][p][q][s]=++nbvar;
					define_var_and(idx2[p][q][r][s],{get(p,q,s),get(q,r,s),get(r,p,s)});
				}
}

void mk_no6hole(int n){
	nbvar=0;
	nbclauses=0;
	nbliterals=0;
	clauses.clear();
	
	cc_system(n);
	var_pt_inside_triangle(n);
	
	// Restriction: No 6-hole.
	for(int p=1;p<=n;p++)
		for(int q=p+1;q<=n;q++)
			for(int r=p+1;r<=n;r++)
				for(int s=p+1;s<=n;s++){
					set<int> pts={p,q,r,s};
					if(pts.size()!=4) continue;
					
					vector<int> left_empty_cond;
					vector<int> right_empty_cond;
					for(int t=1;t<=n;t++){
						if(pts.find(t)!=pts.end()) continue;
						
						define_var_and(++nbvar,{get(p,q,t),get(r,s,t),get(p,s,t)});
						left_empty_cond.push_back(-nbvar);
						
						define_var_and(++nbvar,{get(p,q,t),get(r,s,t),get(r,q,t)});
						right_empty_cond.push_back(-nbvar);
					}
					int left_empty=++nbvar;
					define_var_and(left_empty,left_empty_cond);
					int right_empty=++nbvar;
					define_var_and(right_empty,right_empty_cond);
					
					vector<int> tmp={-get(p,q,r),-get(q,r,s),-get(r,s,p),-get(s,p,q),left_empty,right_empty};
					for(int t=1;t<=n;t++){
						if(pts.find(t)!=pts.end()) continue;
						tmp.push_back(idx2[p][q][r][t]);
						tmp.push_back(idx2[p][r][s][t]);
					}
					new_clauses(tmp);
				}
		
	string s=to_string(n)+"pts-no-6hole.sat";
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
	for(int n=20;n<=30;n++) mk_no6hole(n);
} 

