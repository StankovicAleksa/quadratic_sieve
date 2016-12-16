#include <cstdio>
#include <vector>
using namespace std;
vector <int> vec;
vector <int> ::iterator it;
int main(){
	vec.push_back(2);
	it=vec.begin();
	printf("%d\n",*it);
	return 0;
}
