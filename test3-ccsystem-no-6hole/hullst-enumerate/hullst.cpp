#include<iostream>
#include<cstdio>
#include<cstring>
using namespace std;

int cnt;

void check(int n,string s,bool smallest_hulls,bool hulls_with_obstacle){
	if(hulls_with_obstacle){
		if(s[s.size()-1]=='1' || s[s.size()-1]=='2');
		else cout<<++cnt<<" "<<n<<" "<<s<<"\n";
	}
	if(smallest_hulls) cout<<++cnt<<" "<<n<<" "<<s<<"0\n";
}

void dfs(int tot,int n,string s,bool smallest_hulls,bool hulls_with_obstacle){
	if(n<=2){
		if(n>0) s.push_back('0'+n);
		check(tot,s,smallest_hulls,hulls_with_obstacle);
	}else{
		//for(int i=3;i<=8 && i<=n;i++){
		for(int i=min(8,n);i>=3;i--){
			string t=s;
			t.push_back('0'+i);
			dfs(tot,n-i,t,smallest_hulls,hulls_with_obstacle);
		}
	}
}

int main(){
<<<<<<< HEAD
	//freopen("n_le_24.txt","w",stdout);
	//cnt=0;
	//for(int n=1;n<=24;n++){
	//	dfs(n,n,"",true,true);
	//}
	//fclose(stdout);
	
	for(int n=24;n<=30;n++){
		freopen(("n_"+to_string(n)+"_smallest_hulls.txt").c_str(),"w",stdout);
		cnt=0;
		dfs(n,n,"",true,false);
		fclose(stdout);
	} 
=======
	freopen("n_le_24.txt","w",stdout);
	cnt=0;
	for(int n=1;n<=24;n++){
		dfs(n,n,"",true,true);
	}
	fclose(stdout);
	
	//for(int n=20;n<=30;n++){
	//	freopen(("n_"+to_string(n)+"_smallest_hulls.txt").c_str(),"w",stdout);
	//	cnt=0;
	//	dfs(n,n,"",true,false);
	//	fclose(stdout);
	//} 
>>>>>>> 865821a6b84c745c90c8bf7efc41dc687dd6af28
}

