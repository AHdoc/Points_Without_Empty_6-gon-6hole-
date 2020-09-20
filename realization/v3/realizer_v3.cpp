#include<iostream>
#include<cstdio>
#include<cstring>
#include<cmath>
#include<ctime>
#include<random>
#include<vector>
#include<set>
#include<algorithm>
#include<assert.h>
using namespace std;

typedef long long LL;

LL random(LL a,LL b){
	assert(a<=b);
	LL x=rand()*32768;
	x=(x+rand())*32768;
	x=(x+rand())*32768;
	x=(x+rand())%(b-a+1);
	x+=a;
	return x;
}

const int MAXN=30;

#define x first
#define y second

int n;
LL tot_achievement,achievement[MAXN+1];
LL radius[MAXN+1],lvl[MAXN+2];

void dfs(int i,vector<pair<LL,LL>> pt){
	++tot_achievement; ++achievement[i-1];
	if(tot_achievement%1000000==0){
		cerr<<"achieve = "<<tot_achievement<<"     ";
		for(int j=11;j<=30;j++){
			if(achievement[j]==0) break;
			cerr<<j<<":"<<achievement[j]<<" ";
		}
		cerr<<"\n";
	}
	
	if(i==n+1){
		for(int i=1;i<=n;i++) cerr<<"Point #"<<i<<": ("<<pt[i-1].x<<","<<pt[i-1].y<<")\n";
		exit(0);
	}
	
	
}

void Realizer(string pat){
	n=0;
	LL radius0=1LL;
	for(int i=pat.size()-1,j=1,k;i>=0;i--,j=k){
		k=j+pat[i]-'0';
		for(int l=j;l<k;l++){
			lvl[l]=pat.size()-1-i;
			radius[l]=radius0;
		}
		radius0*=10LL;
		n=k-1;
	}
	lvl[0]=0; lvl[n+1]=pat.size();
	cerr<<"pat="<<pat<<"   n="<<n<<"\n"; 
	for(int i=1;i<=n;i++) cerr<<"i="<<i<<":   radius="<<radius[i]<<"   lvl="<<lvl[i]<<"\n";
	cerr<<"lvl[0]="<<lvl[0]<<"   lvl[n+1]="<<lvl[n+1]<<"\n";
	
	tot_achievement=0; for(int i=1;i<=30;i++) achievement[i]=0;
	for(;;) dfs(2,{{0LL,0LL}});
}

int main(){
	Realizer("3477710");
}

