#include<iostream>
#include<cstdio>
#include<cstring>
#include<cmath>
#include<ctime>
#include<random>
#include<vector>
#include<set>
#include<algorithm>
#include<queue>
#include<assert.h>
using namespace std;

typedef long long LL;

LL absLL(LL x){return (x<0?-x:x);}
LL gcd(LL a,LL b){return (b==0?a:gcd(b,a%b));}
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

int n;
LL tot_achievement,achievement[MAXN+1];
LL radius[MAXN+1],lvl[MAXN+2];

#define x first
#define y second

bool on_the_left(pair<LL,LL> a,pair<LL,LL> b,pair<LL,LL> c){return (b.x-a.x)*(c.y-a.y)-(c.x-a.x)*(b.y-a.y)>0;}

namespace geo{
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
	ll dcmp(ll x){
		if(x>0) return 1;if(x<0) return -1;return 0;
	}
	class P{
		public:
			ll x,y;
			P(ll x=0,ll y=0):x(x),y(y){}
			
			friend ll dis2(P A,P B){return sqr(A.x-B.x)+sqr(A.y-B.y);	}
			friend ll Dot(P A,P B) {return A.x*B.x+A.y*B.y; }
				
			friend P operator- (P A,P B) { return P(A.x-B.x,A.y-B.y); }
			P(P A,P B):x(B.x-A.x),y(B.y-A.y){}
			friend P operator+ (P A,P B) { return P(A.x+B.x,A.y+B.y); }
			friend P operator* (P A,double p) { return P(A.x*p,A.y*p); }
			friend P operator/ (P A,double p) { return P(A.x/p,A.y/p); }
			friend bool operator< (const P& a,const P& b) {return dcmp(a.x-b.x)<0 ||(dcmp(a.x-b.x)==0&& dcmp(a.y-b.y)<0 );}
		}; 
	P read_point() {
		P a;
		scanf("%lld%lld",&a.x,&a.y);
		return a;	
	} 
	bool operator==(const P& a,const P& b) {
		return dcmp(a.x-b.x)==0 && dcmp(a.y-b.y) == 0;
	} 
	typedef P V;
	
	ll Cross(V A,V B) {return A.x*B.y - A.y*B.x;}
	ll Area2(P A,P B,P C) {return Cross(B-A,C-A);}
	
	bool OnLeft(P A,P B,P C) {
		return Cross(B-A,C-A)>0;
	} 
	P _p;
	int cmp(P A,P B) //1:a>b 0:a<=b
	{
		ll tmp=Cross(V(_p,A),V(_p,B));
		if (tmp>0) return 1;
		else if (tmp==0) return (-(dis2(_p,A)-dis2(_p,B))>0)?1:0;
		else return 0;
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
	bool OnLeft(Line L,P p) {
		return Cross(L.v,p-L.p)>0;
	} 
	
}

void proceed(int i,int j,vector<geo::P> &vp,vector< queue<int> > &q,vector<vector<int> > &vg ) {
	while(!q[i].empty() && geo::OnLeft(geo::Line(q[i].front(),vp[i]-q[i].front()),vp[j])){
		proceed(q[i].front(),j,vp,q,vg);
		q[i].pop();
	}
	vg[i].pb(j); //add_edge(i,j)
	q[j].push(i);
}
int L[MAXN];
vector<vector<int> > get_reverese_graph(vector<vector<int> > vg) {
	vector<vector<int> > rG;
	rG.resize(vg.size());
	Rep(i,vg.size()) {
		Rep(j,vg[i].size()) {
			rG[vg[i][j]].pb(i);
		}
	}
	return rG;
}
void treat(geo::P p,vector<geo::P> vp, vector<int> &i, vector<int> &o) {
	int imax=SI(i),omax=SI(o);
	int l=omax,m=0;
	RepD(j,imax-1) {
		L[i[j]]=m+1;
		while(l>0 && geo::OnLeft(vp[i[j]],p,vp[o[l]])) {
			if(L[o[l]]>m) {
				m=L[o[l]];
				L[i[j]]=m+1;
			}
			l--;
		}
	}	
}
void maxchain(vector<geo::P> &vp, vector<vector<int> > &vg,vector<vector<int> > &rG){
	for(int i=n-1;i>=1;i--) treat(vp[i],vp,rG[i],vg[i]);
}
bool find6hole(vector<pair<LL,LL>> pt,pair<LL,LL> p){
	vector<geo::P> vp;
	for(auto p:pt) {
		vp.pb(geo::P(p.x,p.y));
	}
	geo::P _p=geo::P(p.x,p.y);
	int n=pt.size();
	geo::_p =_p;
	sort(vp.begin(),vp.end(),geo::cmp);
	vector<queue<int> > q;
	q.resize(n);
	vector<vector<int> > vg;
	Rep(i,n-1) {
		proceed(i,i+1,vp,q,vg);
	}
	vector<vector<int> > rG=get_reverese_graph(vg);
	maxchain(vp,vg,rG);
}



bool check(int i,int ii,vector<pair<LL,LL>> pt,pair<LL,LL> p){
	set<pair<LL,LL>> ang;
	for(int j=1;j<i;j++){
		LL x=pt[j-1].x-p.x,y=pt[j-1].y-p.y;
		if(x==0 || y==0) return false;
		LL g=gcd(absLL(x),absLL(y)); x/=g; y/=g;
		if(x<0){x=-x; y=-y;}
		else if(x==0) y=-y;
		ang.insert({x,y}); 
	}
	if(ang.size()!=i-1) return false;
	if(lvl[i]==lvl[i-1]+1){
		for(int j=1;j<i;j++)
			if(p.x<=pt[j-1].x) return false;
	}else if(lvl[i]+1==lvl[i+1]){
		for(int j=1;j<i;j++)
			if((j!=i-1 && !on_the_left(p,pt[j-1],pt[i-2])) || (j!=ii && !on_the_left(p,pt[ii-1],pt[j-1]))) return false;
	}else{
		for(int j=1;j<i;j++)
			if(j!=i-1 && !on_the_left(p,pt[j-1],pt[i-2])) return false;
	}
	return find6hole(pt,p);
}

void dfs(int i,int ii,vector<pair<LL,LL>> pt){ // points numbered from ii to i are of the same level
	++tot_achievement; ++achievement[i-1];
<<<<<<< HEAD
	if(tot_achievement%1000000==0){
		cerr<<"achieve = "<<tot_achievement<<"     ";
		for(int j=11;j<=30;j++){
=======
	if(tot_achievement%100000==0){
		cerr<<"achieve = "<<tot_achievement<<"     ";
		for(int j=1;j<=30;j++){
>>>>>>> 7250bc7c6d1f6c838f92613db5709e212bdf54b5
			if(achievement[j]==0) break;
			cerr<<j<<":"<<achievement[j]<<" ";
		}
		cerr<<"\n";
	}
	
	if(i==n+1){
		for(int i=1;i<=n;i++) cerr<<"Point #"<<i<<": ("<<pt[i-1].x<<","<<pt[i-1].y<<")\n";
		exit(0);
	}
	
	for(int t=1;t<=100;t++){
		pair<LL,LL> p={random(-radius[i],radius[i]),random(-radius[i],radius[i])}; 
		if(check(i,ii,pt,p)){
			vector<pair<LL,LL>> pt2=pt;
			pt2.push_back(p);
			dfs(i+1,lvl[i]==lvl[i+1]?ii:i+1,pt2);
		}
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
	for(;;) dfs(2,lvl[1]==lvl[2]?1:2,{{0LL,0LL}}); // must fix the kernel point
}

int main(){
<<<<<<< HEAD
	Realizer("3477710");
=======
	//Realizer("333330");
	Realizer("3333330");
	//Realizer("8730");
	//Realizer("88510");
	//Realizer("3477710");
	
	//check({{0,0},{59,-35},{-99,81},{-77,6},{16,-87},{96,-82}});
	//cout<<ahdoc::find6hole({{0,0},{59,-35},{-99,81},{-77,6},{16,-87},{96,-82}},{92,-73})<<"\n";
	//check({{0,0},{59,-35},{-99,81},{-77,6},{16,-87},{96,-82},{92,-73}});
>>>>>>> 7250bc7c6d1f6c838f92613db5709e212bdf54b5
}

