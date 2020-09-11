#include<iostream>
#include<cstdio>
#include<cstring>
using namespace std;

int cnt;

void check(int n,string s,bool smallest_hulls,bool hulls_with_obstacle,int startswith){
	if(s[0]-'0'!=startswith);
	else{
		if(hulls_with_obstacle){
			if(s[s.size()-1]=='1' || s[s.size()-1]=='2');
			else cout<<++cnt<<" "<<n<<" "<<s<<"\n";
		}
		if(smallest_hulls) cout<<++cnt<<" "<<n<<" "<<s<<"0\n";
	}
}

void dfs(int tot,int n,string s,bool smallest_hulls,bool hulls_with_obstacle,int startswith){
	if(n<=2){
		if(n>0) s.push_back('0'+n);
		check(tot,s,smallest_hulls,hulls_with_obstacle,startswith);
	}else{
		//for(int i=3;i<=8 && i<=n;i++){
		for(int i=min(8,n);i>=3;i--){
			string t=s;
			t.push_back('0'+i);
			dfs(tot,n-i,t,smallest_hulls,hulls_with_obstacle,startswith);
		}
	}
}

int main(){
	freopen("n_g_24_le_30_startswith(8).txt","w",stdout);
	int st=8;
	cnt=0;
	for(int n=25;n<=30;n++){
		if(n<=27) dfs(n,n,"",true,true,st);
		else dfs(n,n,"",true,false,st);
	}
	//fclose(stdout);
	
	//for(int n=24;n<=30;n++){
	//	freopen(("n_"+to_string(n)+"_smallest_hulls.txt").c_str(),"w",stdout);
	//	cnt=0;
	//	dfs(n,n,"",true,false);
	//	fclose(stdout);
	//} 
}

