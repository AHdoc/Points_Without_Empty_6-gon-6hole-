namespace ahdoc{
	typedef long long LL;
	
	#define x first
	#define y second
	#define c1 first.first
	#define c2 first.second
	#define c3 second
	#define MP3(a,b,c) make_pair(make_pair(a,b),c)
	
	bool on_the_left(pair<LL,LL> a,pair<LL,LL> b,pair<LL,LL> c){return (b.x-a.x)*(c.y-a.y)-(c.x-a.x)*(b.y-a.y)>0;}
	
	vector<pair<LL,LL>> pt;
	vector<queue<pair<pair<int,int>,int>>> Q;
	vector<pair<pair<int,int>,int>> C;
	
	bool proceed(int i,int j){
		//cerr<<">   i="<<i<<"   j="<<j<<"\n";
		while(!Q[i].empty() && on_the_left(pt[Q[i].front().c1],pt[i],pt[j])){
			if(proceed(Q[i].front().c1,j)) return true;
			auto tmp=Q[i].front();
			if(tmp.c1>C[i].c1) C[i].c1=tmp.c1;
			if(tmp.c2>C[i].c2) C[i].c2=tmp.c2;
			if(tmp.c3>C[i].c3) C[i].c3=tmp.c3;
			Q[i].pop();
		}
		if(C[i].c3!=-1 && on_the_left(pt[C[i].c3],pt[j],make_pair(0LL,0LL))) return true;
		Q[j].push(MP3(i,C[i].c1,C[i].c2));
		//cerr<<"<   i="<<i<<"   j="<<j<<"\n";
		//for(int k=0;k<pt.size();k++){
		//	cerr<<"       Q["<<k<<"]={ ";
		//	for(auto x:Q[k]) cerr<<"("<<x.c1<<","<<x.c2<<","<<x.c3<<") ";
		//	cerr<<"}\n";
		//	cerr<<"       C["<<k<<"]=("<<C[k].c1<<","<<C[k].c2<<","<<C[k].c3<<")\n";
		//}
		return false;
	}
	
	bool solve(vector<pair<LL,LL>> _pt){
		pt=_pt;
		int n=pt.size();
		//cerr<<"n="<<n<<"\n";
		//for(int i=0;i<n;i++) cerr<<"   "<<i<<": ("<<pt[i].x<<","<<pt[i].y<<")\n";
		Q.clear(); C.clear(); Q.resize(n); C.resize(n);
		for(int i=0;i<n;i++)
			C[i].c1=C[i].c2=C[i].c3=-1;
		for(int i=0;i+1<n;i++)
			if(proceed(i,i+1))
				return true;
		return false; 
	}
	
	bool find6hole(vector<pair<LL,LL>> pt,pair<LL,LL> p){
		int n=pt.size();
		for(int i=0;i<n;i++){
			pt[i].x-=p.x;
			pt[i].y-=p.y;
		}
		sort(pt.begin(),pt.end(),[](pair<LL,LL> p1,pair<LL,LL> p2){return atan2(p1.y,p1.x)<atan2(p2.y,p2.x);});
		if(solve(pt)) return true;
		for(int i=0;i<n;i++){
			pt[i].x*=-1;
			pt[i].y*=-1;
		}
		sort(pt.begin(),pt.end(),[](pair<LL,LL> p1,pair<LL,LL> p2){return atan2(p1.y,p1.x)<atan2(p2.y,p2.x);});
		if(solve(pt)) return true;
		return false;
	}
}
