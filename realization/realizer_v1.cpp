#include<iostream>
#include<cstdio>
#include<cstring>
#include<cmath>
#include<ctime>
#include<random>
#include<vector>
#include<set>
#include<algorithm>

using namespace std;

double eps=1e-9;

default_random_engine gen;

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

double Area_Polygon(vector<Tpoint> pol){
	double area;
	for(int i=2;i<pol.size();i++) // Triangle (0,i-1,i)
		area+=(pol[i-1]-pol[0])^(pol[i]-pol[0]);
	return area;
}

void Print_Q(vector<Tline> Q,int head,int tail){
	cerr<<head<<"---"<<tail<<"   ";
	for(int i=head;i<=tail;i++){
		Tline l=Q[i];
		cerr<<fixed<<"("<<l.s.x<<","<<l.s.y<<")-("<<l.e.x<<","<<l.e.y<<")   ";
	} 
	cerr<<"\n";
}

vector<Tpoint> HPI(vector<Tline> _line){
	sort(_line.begin(),_line.end(),[](Tline a,Tline b){return a.k<b.k;});
	vector<Tline> line;
	line.push_back(_line[0]);
	for(int i=1;i<_line.size();i++)
		if(fabs(_line[i].k-_line[i-1].k)>eps)
			line.push_back(_line[i]);
		else if(((line[line.size()-1].e-_line[i].s)^(_line[i].e-_line[i].s))>-eps){
			line.pop_back();
			line.push_back(_line[i]);
		}
	sort(line.begin(),line.end(),[](Tline a,Tline b){return a.k<b.k;});
	
	//for(Tline l:line) cout<<fixed<<l.k<<"   ("<<l.s.x<<","<<l.s.y<<")->("<<l.e.x<<","<<l.e.y<<")\n";
	
	int head=0,tail=1;
	vector<Tline> Q;
	Q.resize(line.size());
	Q[0]=line[0];
	Q[1]=line[1];
	vector<Tpoint> res;
	//cerr<<"---\n";
	for(int i=2;i<line.size();i++){
		if(fabs((Q[tail].e-Q[tail].s)^(Q[tail-1].e-Q[tail-1].s))<eps || fabs((Q[head].e-Q[head].s)^(Q[head+1].e-Q[head+1].s))<eps)
			return res;
		while(head<tail && (((Q[tail]&Q[tail-1])-line[i].s)^(line[i].e-line[i].s))>-eps)
			--tail;
		while(head<tail && (((Q[head]&Q[head+1])-line[i].s)^(line[i].e-line[i].s))>-eps)
			++head;
		Q[++tail]=line[i];
		//cerr<<"i="<<i<<"  "; Print_Q(Q,head,tail);
	}
	while(head<tail && (((Q[tail]&Q[tail-1])-Q[head].s)^(Q[head].e-Q[head].s))>-eps)
		--tail;
	while(head<tail && (((Q[head]&Q[head+1])-Q[tail].s)^(Q[tail].e-Q[tail].s))>-eps)
		++head;
	//Print_Q(Q,head,tail);
	if(tail<=head+1);
	else{
		for(int i=head;i<tail;i++)
			res.push_back(Q[i]&Q[i+1]);
		if(head<tail-1)
			res.push_back(Q[head]&Q[tail]);
	}
	
	double area=Area_Polygon(res);
	if(area<eps) res.clear();
	
	bool exist_nan=false;
	for(Tpoint x:res)
		if(isnan(x.x) || isnan(x.y))
			exist_nan=true;
	
	if(exist_nan){
		cout<<"_line:\n";
		for(Tline QQ:_line)
			cout<<fixed<<"("<<QQ.s.x<<","<<QQ.s.y<<") --- ("<<QQ.e.x<<","<<QQ.e.y<<")\n";
		cout<<"line:\n";
		for(Tline QQ:line)
			cout<<fixed<<"("<<QQ.s.x<<","<<QQ.s.y<<") --- ("<<QQ.e.x<<","<<QQ.e.y<<")\n";
		cout<<"Q[]:\n";
		for(int i=head;i<=tail;i++)
			cout<<fixed<<"("<<Q[i].s.x<<","<<Q[i].s.y<<") --- ("<<Q[i].e.x<<","<<Q[i].e.y<<")\n";
		cout<<"res:\n";
		for(Tpoint x:res)
			cout<<fixed<<"("<<x.x<<","<<x.y<<") --- ";
		cout<<"\n";
		exit(1); 
	}
	
	return res;
}

void HPI_test(){
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
	int x=(rand()*32768+rand())%100;
	return (double)(x+1)/101; 
}

int Weighted_Random(vector<double> weights){
	discrete_distribution<int> dist(weights.begin(),weights.end());
	return dist(gen);
}

Tpoint PointPicking_Triangle(vector<Tpoint> tri){
	double a=Realnum_Inside01(),b=Realnum_Inside01();
	if(a+b<1){
		Tpoint ret {tri[0].x+a*(tri[1].x-tri[0].x)+b*(tri[2].x-tri[0].x),tri[0].y+a*(tri[1].y-tri[0].y)+b*(tri[2].y-tri[0].y)};
		if(isnan(ret.x) || isnan(ret.y)){
			cerr<<fixed<<"("<<tri[0].x<<","<<tri[0].y<<")\n";
			cerr<<fixed<<"("<<tri[1].x<<","<<tri[1].y<<")\n";
			cerr<<fixed<<"("<<tri[2].x<<","<<tri[2].y<<")\n";
			cerr<<fixed<<"a="<<a<<"   b="<<b<<"\n";
			exit(1);
		}
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

int n,tot_achievement,achievement[MAXN+1];
int idx[MAXN+1][MAXN+1][MAXN+1];
bool var[MAXN*MAXN*MAXN+1];

bool query(int i,int j,int k){
	if(idx[i][j][k]>0) return var[idx[i][j][k]];
	else return (!var[-idx[i][j][k]]);
}

void dfs(int i,int n,Tpoint A,Tpoint B,Tpoint C,Tpoint D,vector<Tpoint> pt){
	++tot_achievement;
	++achievement[i-1];
	if(tot_achievement%1000000==0){
		cerr<<"achieve = "<<tot_achievement<<"     ";
		for(int j=1;j<=30;j++){
			if(achievement[j]==0) break;
			cerr<<j<<":"<<achievement[j]<<" ";
		}
		cerr<<"\n";
	}
	
	if(i==n+1){
		for(int i=1;i<=n;i++)
			cerr<<fixed<<"Point #"<<i<<": Tpoint("<<pt[i-1].x<<","<<pt[i-1].y<<")\n";
		// Check
		for(int i=1;i<=n;i++)
			for(int j=1;j<=n;j++)
				for(int k=1;k<=n;k++){
					set<int> ijk={i,j,k};
					if(ijk.size()!=3) continue;
					double crossprod=(pt[j-1]-pt[i-1])^(pt[k-1]-pt[i-1]);
					
					if(crossprod>0 && query(i,j,k));
					else if(crossprod<0 && (!query(i,j,k)));
					else{
						cerr<<fixed<<"point i = ("<<pt[i-1].x<<","<<pt[i-1].y<<")\n";
						cerr<<fixed<<"point j = ("<<pt[j-1].x<<","<<pt[j-1].y<<")\n";
						cerr<<fixed<<"point k = ("<<pt[k-1].x<<","<<pt[k-1].y<<")\n";
						cerr<<fixed<<"crossprod = "<<crossprod<<"\n";
						cerr<<"query(i,j,k) = "<<query(i,j,k)<<"\n";
						exit(1);
					}
				}
		exit(0);
	}
	
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
	if(res.size()<3)
		return;
	double area=Area_Polygon(res);
		
	for(int t=1;t<=((i==9 || i==11 || i==17)?100:1);){
		Tpoint newP=PointPicking_Polygon(res);
		if(fabs(A.y-newP.y)<eps) continue;
		if(fabs(B.y-newP.y)<eps) continue;
		if(fabs(C.y-newP.y)<eps) continue;
		if(fabs(D.y-newP.y)<eps) continue;
		for(int j=1;j<i;j++)
			if(fabs(pt[j-1].y-newP.y)<eps) continue;
		for(int j=1;j<i;j++)
			for(int k=j+1;k<i;k++)
				if(fabs((pt[k-1]-pt[j-1])^(newP-pt[j-1]))<eps)
					continue;
		
		vector<Tpoint> pt2=pt;
		pt2.push_back(newP);
		dfs(i+1,n,A,B,C,D,pt2);
		++t;
	}
}

void Realizer(){
	freopen("88510.out","r",stdin);
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
	while(cin>>x && x!=0){
		if(abs(x)<=tot){
			if(x>0) var[x]=true;
			else var[abs(x)]=false;
		}
	}
	fclose(stdin);
	
	Tpoint A(-MAXR,-MAXR);
	Tpoint B(MAXR,-MAXR+1);
	Tpoint C(MAXR,MAXR+1);
	Tpoint D(-MAXR,MAXR);
	Tpoint P1(-MAXR+1,MAXR*0.1);
	Tpoint P2(-MAXR+1,-MAXR*0.1);
	
	tot_achievement=0;
	for(int i=1;i<=30;i++) achievement[i]=0; achievement[0]=1; achievement[1]=1; achievement[2]=1;
	for(int i=1;;i++)
		dfs(3,n,A,B,C,D,{P1,P2});
}

int main(){
	cerr.precision(17);
	cout.precision(17);
	
	//HPI_test();
	Realizer();
}


/*
hullst: 88510

	vector<Tpoint> pt={
		Tpoint(-99999.00000000000000000,10000.00000000000000000),
		Tpoint(-99999.00000000000000000,-10000.00000000000000000),
		Tpoint(-6930.69306930692528113,12871.75247524752012396),
		Tpoint(-7146.35820017645164626,53115.87199974375835154),
		Tpoint(-28336.02431590079868329,72148.12528910698893014),
		Tpoint(-36850.53627836801751982,73675.08070030220551416),
		Tpoint(-77491.28025763612822630,66111.91715940831636544),
		Tpoint(-85736.82115335357957520,52753.32605842949124053)};
	//double minx=pt[0].x,maxx=pt[0].x,miny=pt[0].y,maxy=pt[0].y;
	//for(Tpoint p:pt){
	//	minx=fmin(minx,p.x); maxx=fmax(maxx,p.x);
	//	miny=fmin(miny,p.y); maxy=fmax(maxy,p.y);
	//}
	//for(int i=0;i<pt.size();i++){
	//	pt[i].x=-MAXR+2.0*MAXR*(pt[i].x-minx)/(maxx-minx);
	//	pt[i].y=-MAXR+2.0*MAXR*(pt[i].y-miny)/(maxy-miny);
	//}
	for(int i=1;i<=pt.size();i++)
		for(int j=1;j<=pt.size();j++)
			for(int k=1;k<=pt.size();k++){
				set<int> ijk={i,j,k};
				if(ijk.size()!=3) continue;
				double crossprod=(pt[j-1]-pt[i-1])^(pt[k-1]-pt[i-1]);
				
				if(crossprod>0 && query(i,j,k));
				else if(crossprod<0 && (!query(i,j,k)));
				else{
					cerr<<fixed<<"point i = ("<<pt[i-1].x<<","<<pt[i-1].y<<")\n";
					cerr<<fixed<<"point j = ("<<pt[j-1].x<<","<<pt[j-1].y<<")\n";
					cerr<<fixed<<"point k = ("<<pt[k-1].x<<","<<pt[k-1].y<<")\n";
					cerr<<fixed<<"crossprod = "<<crossprod<<"\n";
					cerr<<"query(i,j,k) = "<<query(i,j,k)<<"\n";
					exit(1);
				}
			}
	vector<Tline> line;
	line.push_back(Tline(A,B));
	line.push_back(Tline(B,C));
	line.push_back(Tline(C,D));
	line.push_back(Tline(D,A));
	for(int j=1;j<=pt.size();j++)
		for(int k=j+1;k<=pt.size();k++)
			if(query(j,k,pt.size()+1))
				line.push_back(Tline(pt[j-1],pt[k-1]));
			else
				line.push_back(Tline(pt[k-1],pt[j-1]));
	vector<Tpoint> res=HPI(line);
*/

