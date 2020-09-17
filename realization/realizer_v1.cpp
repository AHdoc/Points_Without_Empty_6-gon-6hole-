#include<iostream>
#include<cstdio>
#include<cstring>
#include<cmath>
#include<ctime>
#include<random>
#include<vector>
#include<algorithm>

using namespace std;

double eps=1e-9;

default_random_engine gen;

int sign(double x){
	if(fabs(x)<eps) return 0;
	return x>0?1:-1; 
}

struct Tpoint{
	double x,y;
	Tpoint(){}
	Tpoint(double _x,double _y){x=_x; y=_y;}
	Tpoint operator -(const Tpoint &b)const{return Tpoint(x-b.x,y-b.y);}
	double operator ^(const Tpoint &b)const{return x*b.y-y*b.x;}
	double operator *(const Tpoint &b)const{return x*b.x+y*b.y;}
};

struct Tline{
	Tpoint s,e;
	double k;
	Tline(){}
	Tline(Tpoint _s,Tpoint _e){
		s=_s; e=_e;
		k=atan2(e.y-s.y,e.x-s.x);
	}
	Tpoint operator &(const Tline &b)const{
		Tpoint res=s;
		double t=((s-b.s)^(b.s-b.e))/((s-e)^(b.s-b.e));
		res.x+=(e.x-s.x)*t;
		res.y+=(e.y-s.y)*t;
		return res;
	}
};

vector<Tpoint> HPI(vector<Tline> line){
	sort(line.begin(),line.end(),[](Tline a,Tline b){return a.k<b.k;});
	
	//for(Tline l:line) cerr<<"("<<l.s.x<<","<<l.s.y<<")->("<<l.e.x<<","<<l.e.y<<") --> "; cerr<<"\n";
	
	int head=0,tail=1;
	vector<Tline> Q;
	Q.resize(line.size());
	Q[0]=line[0];
	Q[1]=line[1];
	vector<Tpoint> res;
	for(int i=2;i<line.size();i++){
		if(fabs((Q[tail].e-Q[tail].s)^(Q[tail-1].e-Q[tail-1].s))<eps || fabs((Q[head].e-Q[head].s)^(Q[head+1].e-Q[head+1].s))<eps)
			return res;
		while(head<tail && (((Q[tail]&Q[tail-1])-line[i].s)^(line[i].e-line[i].s))>eps)
			--tail;
		while(head<tail && (((Q[head]&Q[head+1])-line[i].s)^(line[i].e-line[i].s))>eps)
			++head;
		Q[++tail]=line[i];
	}
	while(head<tail && (((Q[tail]&Q[tail-1])-Q[head].s)^(Q[head].e-Q[head].s))>eps)
		--tail;
	while(head<tail && (((Q[head]&Q[head-1])-Q[tail].s)^(Q[tail].e-Q[tail].e))>eps)
		++head;
	if(tail<=head+1)
		return res;
	for(int i=head;i<tail;i++)
		res.push_back(Q[i]&Q[i+1]);
	if(head<tail-1)
		res.push_back(Q[head]&Q[tail]);
	return res;
}

void HPI_test(){
	cout.precision(6);
	vector<Tline> line0={
		Tline(Tpoint(0,0),Tpoint(1+0.01,0-0.01)),
		Tline(Tpoint(1+0.01,0-0.01),Tpoint(1-0.01,1+0.01)),
		Tline(Tpoint(1-0.01,1+0.01),Tpoint(0-0.01,1-0.01))};
	vector<Tline> line1={
		Tline(Tpoint(0,0),Tpoint(1+0.01,0-0.01)),
		Tline(Tpoint(1+0.01,0-0.01),Tpoint(1-0.01,1+0.01)),
		Tline(Tpoint(1-0.01,1+0.01),Tpoint(0-0.01,1-0.01)),
		Tline(Tpoint(0-0.01,1-0.01),Tpoint(0,0))};
	vector<Tline> line2={
		Tline(Tpoint(-1,0),Tpoint(2,0)),
		Tline(Tpoint(1,-1),Tpoint(1,2)),
		Tline(Tpoint(2,1),Tpoint(-1,1)),
		Tline(Tpoint(0,2),Tpoint(0,-1)),
		Tline(Tpoint(0.5,1),Tpoint(0,0.5))};
	for(Tpoint x:HPI(line2)){
		cout<<fixed<<"("<<x.x<<","<<x.y<<")-";
	}
	cout<<"\n";
}

double Realnum_Inside01(){
	int x=rand()*32768+rand();
	return (double)(x+1)/(32768*32768+1); 
}

int Weighted_Random(vector<double> weights){
	discrete_distribution<int> dist(weights.begin(),weights.end());
	return dist(gen);
}

Tpoint PointPicking_Triangle(vector<Tpoint> tri){
	double a=Realnum_Inside01(),b=Realnum_Inside01();
	if(a+b<1){
		Tpoint ret {tri[0].x+a*(tri[1].x-tri[0].x)+b*(tri[2].x-tri[0].x),tri[0].y+a*(tri[1].y-tri[0].y)+b*(tri[2].y-tri[0].y)};
		return ret;
	}else
		return PointPicking_Triangle(tri);
}

Tpoint PointPicking_Polygon(vector<Tpoint> pol){
	vector<double> areas;
	for(int i=2;i<pol.size();i++) // Triangle (0,i-1,i)
		areas.push_back((pol[i-1]-pol[0])^(pol[i]-pol[0]));
	int j=2+Weighted_Random(areas);
	return PointPicking_Triangle({pol[0],pol[j-1],pol[j]}); 
}

/**************************************************************/

const int MAXN=30;
const double MAXR=1e5;

int n;
int idx[MAXN+1][MAXN+1][MAXN+1];
bool var[MAXN*MAXN*MAXN+1];

bool query(int i,int j,int k){
	if(idx[i][j][k]>0) return var[idx[i][j][k]];
	else return (!var[-idx[i][j][k]]);
}

bool Solve(Tpoint A,Tpoint B,Tpoint C,Tpoint D,Tpoint P1,Tpoint P2){
	vector<Tpoint> pt;
	pt.push_back(P1);
	pt.push_back(P2);
	for(int i=3;i<=n;i++){
		vector<Tline> line;
		line.push_back(Tline(A,B));
		line.push_back(Tline(B,C));
		line.push_back(Tline(C,D));
		line.push_back(Tline(D,A));
		for(int j=1;j<i;j++)
			for(int k=j+1;k<i;k++)
				if(query(j,k,i))
					line.push_back(Tline(pt[j-1],pt[k-1]));
				else
					line.push_back(Tline(pt[k-1],pt[j-1]));
		vector<Tpoint> res=HPI(line);
		if(res.empty()) return false;
		pt.push_back(PointPicking_Polygon(res));
	}
	
	cout.precision(17);
	for(int i=1;i<=n;i++)
		cerr<<fixed<<"Point #"<<i<<": "<<pt[i-1].x<<" "<<pt[i-1].y<"\n";
	return true;
}

void Realizer(){
	cin>>n;
	int tot=0;
	for(int i=1;i<=n;i++)
		for(int j=i+1;j<=n;j++)
			for(int k=j+1;k<=n;k++){
				int x=++tot;
				idx[i][j][k]=idx[j][k][i]=idx[k][i][j]=x;
				idx[i][k][j]=idx[j][i][k]=idx[k][j][i]=-x;
			}
	int x;
	while(cin>>x){
		if(abs(x)<=tot){
			if(x>0) var[x]=true;
			else var[abs(x)]=false;
		}
	}
	
	Tpoint A(-MAXR+eps,-MAXR-eps);
	Tpoint B(MAXR+eps,-MAXR+eps);
	Tpoint C(MAXR-eps,MAXR+eps);
	Tpoint D(-MAXR-eps,MAXR-eps);
	Tpoint P1(0,0);
	Tpoint P2(100,1);
	
	for(int i=1;;i++){
		if(Solve(A,B,C,D,P1,P2)) break;
		if(i%100==0) cerr<<"i = "<<i<<"\n";
	}
}

int main(){
	//HPI_test();
	Realizer();
}


/*
hullst: 88510
*/

