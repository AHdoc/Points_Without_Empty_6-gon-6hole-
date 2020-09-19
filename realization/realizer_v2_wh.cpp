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

double eps=1e-3;

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

namespace geo {
	
	#define For(i,n) for(int i=1;i<=n;i++)
	#define Fork(i,k,n) for(int i=k;i<=n;i++)
	#define Rep(i,n) for(int i=0;i<n;i++)
	#define ForD(i,n) for(int i=n;i;i--)
	#define ForkD(i,k,n) for(int i=n;i>=k;i--)
	#define RepD(i,n) for(int i=n;i>=0;i--)
	#define Forp(x) for(int p=Pre[x];p;p=Next[p])
	#define Forpiter(x) for(int &p=iter[x];p;p=Next[p])  
	#define Lson (o<<1)
	#define Rson ((o<<1)+1)
	#define MEM(a) memset(a,0,sizeof(a));
	#define MEMI(a) memset(a,127,sizeof(a));
	#define MEMi(a) memset(a,128,sizeof(a));
	#define INF (2139062143)
	#define F (100000007)
	#define pb push_back
	#define mp make_pair 
	#define fi first
	#define se second
	#define vi vector<int> 
	#define pi pair<int,int>
	#define SI(a) ((a).size())
	#define ALL(x) (x).begin(),(x).end()
	typedef long long ll;
	typedef long double ld;
	typedef unsigned long long ull;
	ll mul(ll a,ll b){return (a*b)%F;}
	ll add(ll a,ll b){return (a+b)%F;}
	ll sub(ll a,ll b){return (a-b+llabs(a-b)/F*F+F)%F;}
	void upd(ll &a,ll b){a=(a%F+b%F)%F;}
	int read()
	{
		int x=0,f=1; char ch=getchar();
		while(!isdigit(ch)) {if (ch=='-') f=-1; ch=getchar();}
		while(isdigit(ch)) { x=x*10+ch-'0'; ch=getchar();}
		return x*f;
	} 
	ll sqr(ll a){return a*a;}
	ld sqr(ld a){return a*a;}
	double sqr(double a){return a*a;}
	const double eps=1e-10;
	int dcmp(double x) {
		if (fabs(x)<eps) return 0; else return x<0 ? -1 : 1; 
	}
	ld PI = 3.141592653589793238462643383;
	class P{
	public:
		double x,y;
		P(double x=0,double y=0):x(x),y(y){}
		friend ld dis2(P A,P B){return sqr(A.x-B.x)+sqr(A.y-B.y);	}
		friend ld Dot(P A,P B) {return A.x*B.x+A.y*B.y; }
		friend ld Length(P A) {return sqrt(Dot(A,A)); }
		friend ld Angle(P A,P B) {
			if (dcmp(Dot(A,A))==0||dcmp(Dot(B,B))==0||dcmp(Dot(A-B,A-B))==0) return 0;
			return acos(max((ld)-1.0, min((ld)1.0, Dot(A,B) / Length(A) / Length(B) )) ); 
		}
			
		friend P operator- (P A,P B) { return P(A.x-B.x,A.y-B.y); }
		friend P operator+ (P A,P B) { return P(A.x+B.x,A.y+B.y); }
		friend P operator* (P A,double p) { return P(A.x*p,A.y*p); }
		friend P operator/ (P A,double p) { return P(A.x/p,A.y/p); }
		friend bool operator< (const P& a,const P& b) {return dcmp(a.x-b.x)<0 ||(dcmp(a.x-b.x)==0&& dcmp(a.y-b.y)<0 );}
		
	}; 
	P read_point() {
		P a;
		scanf("%lf%lf",&a.x,&a.y);
		return a;	
	} 
	bool operator==(const P& a,const P& b) {
		return dcmp(a.x-b.x)==0 && dcmp(a.y-b.y) == 0;
	} 
	typedef P V;
	
	double Cross(V A,V B) {return A.x*B.y - A.y*B.x;}
	double Area2(P A,P B,P C) {return Cross(B-A,C-A);}
	V Rotate(V A,double rad) {
		return V(A.x*cos(rad)-A.y*sin(rad),A.x*sin(rad)+A.y*cos(rad));
	} 
	// A 不是 0向量 
	V Normal(V A) { 
		double L = Length(A);
		return V(-A.y/L , A.x/L); 
	}
	
	P GetLineIntersection(P p,V v,P Q,V w){
		V u = p-Q;
		double t = Cross(w,u)/Cross(v,w);
		return p+v*t;
	}
	P GetLineIntersectionB(P p,V v,P Q,V w){
		return GetLineIntersection(p,v-p,Q,w-Q);
	}
	
	double DistanceToLine(P p,P A,P B) {
		V v1 = B-A, v2 = p-A;
		return fabs(Cross(v1,v2))/Length(v1);
	}
	double DistanceToSegment(P p,P A,P B) {
		if (A==B) return Length(p-A);
		V v1 = B-A, v2 = p-A, v3 = p - B;
		if (dcmp(Dot(v1,v2))<0) return Length(v2);
		else if (dcmp(Dot(v1,v3))>0 ) return Length(v3);
		else return fabs(Cross(v1,v2) ) / Length(v1);
	}
	P GetLineProjection(P p,P A,P B) {
		V v=B-A;
		return A+v*(Dot(v,p-A)/Dot(v,v));
	}
	//规范相交-线段相交且交点不在端点 
	bool SegmentProperIntersection(P a1,P a2,P b1,P b2) { 
		double  c1 = Cross(a2-a1,b1-a1) , c2 = Cross(a2-a1,b2-a1),
				c3 = Cross(b2-b1,a1-b1) , c4 = Cross(b2-b1,a2-b1);
		return dcmp(c1)*dcmp(c2)<0 && dcmp(c3)*dcmp(c4)<0;
	}
	//点在线段上（不包含端点） 
	bool OnSegment(P p,P a1,P a2) {
		return dcmp(Cross(a1-p,a2-p)) == 0 && dcmp(Dot(a1-p,a2-p))<0;
	}
	double PolygonArea(P *p,int n) {
		double area=0;
		For(i,n-2) area+=Cross(p[i]-p[0],p[i+1]-p[0]);
		return area/2;
	} 
	struct Line{
		P p;
		V v;
		double ang;
		Line(){}
		Line(P p,V v):p(p),v(v) {ang=atan2(v.y,v.x); }
		bool operator<(const Line & L) const {
			return ang<L.ang;
		}
		P point(double a) {
			return p+v*a;
		}
	};
	
	//线上不算 
	bool OnLeft(Line L,P p) {
		return Cross(L.v,p-L.p)>0;
	} 
	P GetIntersection(Line a,Line b) {
		V u=a.p-b.p;
		double t = Cross(b.v,u) / Cross(a.v, b.v);
		return a.p + a.v*t;
	}
	int HalfplaneIntersection(Line *L, int n, P* poly) {
		sort(L,L+n);
		int fi,la;
		P *p = new P[n];
		Line *q = new Line[n];
		q[fi=la=0 ] = L[0];
		For(i,n-1) {
			while (fi < la && !OnLeft(L[i],p[la-1])) la--;
			while (fi < la && !OnLeft(L[i],p[fi])) fi++;
			q[++la] = L[i];
			if (fabs(Cross(q[la].v, q[la-1].v))<eps) {
				la--;
				if (OnLeft(q[la],L[i].p)) q[la] = L[i];
			}
			if (fi<la) p[la-1] = GetIntersection(q[la-1],q[la]);
		}
		while(fi < la && !OnLeft(q[fi],p[la-1])) la--;
		if (la-fi<=1) return 0; 
		p[la] = GetIntersection(q[la],q[fi]);
		
		int m=0;
		Fork(i,fi,la) poly[m++]=p[i];
		return m;
	} 
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
	
	geo::Line l[SI(line)];
	int n=SI(line);
	geo::P p[2*n];
	for(int i=0;i<n;i++) {
		l[i]=geo::Line(geo::P(line[i].s.x,line[i].s.y),geo::P(line[i].e.x-line[i].s.x,line[i].e.y-line[i].s.y));
	}	
	int m=geo::HalfplaneIntersection(l,n,p);
	vector<Tpoint> res;
	for(int i=0;i<m;i++) {
		res.push_back(Tpoint(p[i].x,p[i].y) );
	}	
	
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
//		for(int i=head;i<=tail;i++)
//			cout<<fixed<<"("<<Q[i].s.x<<","<<Q[i].s.y<<") --- ("<<Q[i].e.x<<","<<Q[i].e.y<<")\n";
		for(Tpoint x:res)
			cout<<fixed<<"("<<x.x<<","<<x.y<<") --- ";
		cout<<"res:\n";
		for(Tpoint x:res)
			cout<<fixed<<"("<<x.x<<","<<x.y<<") --- ";
		cout<<"\n";
		exit(1); 
	}
	
	return res;
}

double Realnum_Inside01(){
	int x=(rand()*32768+rand())%100;
	return (double)(x+1)/101; 
}

int Weighted_Random(vector<double> _weights){
	vector<double> weights;
	for(int i=0;i<_weights.size();i++)
		if(_weights[i]<eps) weights.push_back(0.0);
		else weights.push_back(_weights[i]);
	
	int j=0;
	while(weights[j]<eps) ++j;
	cerr<<fixed<<"j="<<j<<"   area="<<weights[j]<<"\n";
	return j;
	
	//discrete_distribution<int> dist(weights.begin(),weights.end());
	//return dist(gen);
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
		for(int j=10;j<=30;j++){
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
	Realizer();
}

