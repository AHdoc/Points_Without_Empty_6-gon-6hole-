#include<iostream>
#include<cstdio>
#include<cstring>
#include<cmath>
#include<ctime>
#include<random>
#include<vector>
#include<set>
#include<algorithm>
#include<unordered_map>
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

struct pair_hash
{
	template <class T1, class T2>
	std::size_t operator() (const std::pair<T1, T2> &pair) const
	{
		return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
	}
};
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
	int Quadrant(P a)
	{
	    if(a.x>0&&a.y>=0) return 1;
	    if(a.x<=0&&a.y>0) return 2;
	    if(a.x<0&&a.y<=0) return 3;
	    return 4;
	}

	int cmp(P A,P B) //1:a>b 0:a<=b
	{
		if(Quadrant(A)!=Quadrant(B))
			return Quadrant(A)<Quadrant(B);
		ll tmp=Cross(V(_p,A),V(_p,B));
		if (tmp>0) return 1;
		else if (tmp==0) return (-(dis2(_p,A)-dis2(_p,B))>0)?1:0;
		else return 0;
	}
	void PolarSort(vector<P> &v,P _p2) {
		_p=_p2;
		sort(ALL(v),cmp);
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
	bool OnRight(Line L,P p) {
		return Cross(L.v,p-L.p)<0;
	} 
	vector<P> v;
	bool proceed(int i,int j,P _p,vector<P> &vp,vector< queue<int> > &q,vector<queue<pair<int,int> > > &q2,
	vector<vector<int> > &vg ,vector<pair<int,int> > &C,vector<vector< pair<int,int>> >  &Ce) {
		while(!q[i].empty() && geo::OnLeft(geo::Line(q[i].front(),vp[i]-q[i].front()),vp[j])){
		 //if k can see i && i can see j && turn)left, then k can see j
			if (proceed(q[i].front(),j,_p,vp,q,q2,vg,C,Ce)) return 1; // add k-j and p-k-j 
				
				auto now=q2[i].front();
				C[i].fi=max(C[i].fi,now.fi);
				C[i].se=max(C[i].se,now.se);
				
				q2[i].pop();
				q[i].pop();
		}
		vg[i].pb(j); //add_edge(i,j)
		
		if(C[i].se!=-1 && geo::OnLeft(geo::Line(vp[C[i].se],vp[i]-vp[C[i].se]),_p ) ) {
			return 1;
		} 
		Ce[i].push_back(mp(C[i].fi,C[i].se));
		q[j].push(i);
		return 0;
	}

	bool find6hole(vector<P> vp,P _p){
	
		int n=vp.size();
		if (n<6) return 0;
		PolarSort(vp,_p);
		
		{
			vector<queue<int> > q;
			vector<queue<pair<int,int> > > q2;
			vector<vector<int> > vg;
			vector<pair<int,int> > C; //C_2,C_3
			vector<vector< pair<int,int>  > > Ce;
			q.resize(n); q2.resize(n);
			C.resize(n); vg.resize(n);
			Ce.resize(n);
			Rep(i,n) C[i]=make_pair(-1,-1);
				
			Rep(i,n-1) {
				if(proceed(i,i+1,_p,vp,q,q2,vg,C,Ce))return 1;
			}
		}
		int id=-1; 
		while(id<n&&OnLeft(_p,vp[0],vp[id])) ++id;
		if(id>=n&&id==-1) return 0;
		
		rotate(vp.begin(),vp.begin()+id,vp.end());

		{
			vector<queue<int> > q;
			vector<queue<pair<int,int> > > q2;
			vector<vector<int> > vg;
			vector<pair<int,int> > C; //C_2,C_3
			vector<vector< pair<int,int>  > > Ce;
			q.resize(n); q2.resize(n);
			C.resize(n); vg.resize(n);
			Ce.resize(n);
			Rep(i,n) C[i]=make_pair(-1,-1);
				
			Rep(i,n-1) {
				if(proceed(i,i+1,_p,vp,q,q2,vg,C,Ce))return 1;
			}
		}
		return 0;
	}
}

bool find6hole(vector<pair<LL,LL>> pt,pair<LL,LL> p){
	vector<geo::P> vp;
	for(auto p:pt) {
		vp.pb(geo::P(p.x,p.y));
	}
	geo::P _p=geo::P(p.x,p.y);
	
	bool b=geo::find6hole(vp,_p);
	cout<<b<<endl;
	return b;
	
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
	if(tot_achievement%100000==0){
		cerr<<"achieve = "<<tot_achievement<<"     ";
		for(int j=1;j<=30;j++){
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
	Realizer("3477710");
}

