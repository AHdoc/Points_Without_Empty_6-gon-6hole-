#include<cstdio>
#include<cstring>
#include<iostream>
#include<vector>
#include<algorithm>

#include<assert.h>

using namespace std;

const int MAXN=50;

int n;
bool va[MAXN+1][MAXN+1][MAXN+1];
bool inhull[MAXN+1];

void findhull(int lvl){
	vector<int> pts;
	
	int cnt=0;
	for(int i=1;i<=n;i++)
		if(!inhull[i])
			++cnt;
	if(cnt==0) return;
	if(cnt<=2){
		pts.clear();
		for(int i=1;i<=n;i++)
			if(!inhull[i])
				pts.push_back(i);
	}else{
		vector<pair<int,int>> hulledges;
		for(int i=1;i<=n;i++)
			for(int j=1;j<=n;j++){
				if(i==j) continue;
				if(inhull[i]) continue;
				if(inhull[j]) continue;
				
				bool hulledge=true;
				for(int k=1;k<=n;k++)
					if(k!=i && k!=j && !inhull[k]){
						if(!va[i][j][k]){
							hulledge=false;
							break;
						}
					}
				if(hulledge) hulledges.push_back(make_pair(i,j));
			}
		int nn=hulledges.size();
		int prev=-1,now=hulledges[0].first,nxt;
		pts.clear();
		pts.push_back(now);
		for(int i=2;i<=nn;i++){
			for(auto pii:hulledges){
				if(pii.first==now && pii.second!=prev){
					nxt=pii.second;
					break;
				}
				if(pii.second==now && pii.first!=prev){
					nxt=pii.first;
					break;
				}
			}
			prev=now;
			now=nxt;
			pts.push_back(now);
		}
	}
	
	cerr<<"lvl = "<<lvl<<" : ";
	for(int x:pts){
		cerr<<x<<" ";
		inhull[x]=true;
	}
	cerr<<"\n";
	
	findhull(lvl+1);
} 

void findhull(){
	for(int i=1;i<=n;i++) inhull[i]=false;
	findhull(1);
}

int main(){
	n=22;
	freopen("22.txt","r",stdin);
	
	int nbvar=0;
	for(int i=1;i<=n;i++)
		for(int j=i+1;j<=n;j++)
			for(int k=j+1;k<=n;k++){
				++nbvar;
				int o;
				cin>>o;
				
				assert(abs(o)==nbvar);
				if(o>0){
					va[i][j][k]=va[j][k][i]=va[k][i][j]=true;
					va[i][k][j]=va[j][i][k]=va[k][j][i]=false;
				}else{
					va[i][j][k]=va[j][k][i]=va[k][i][j]=false;
					va[i][k][j]=va[j][i][k]=va[k][j][i]=true;
				}
			}
	
	findhull();
}
