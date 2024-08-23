#include "solution.hpp"

int main(){
unsigned M,J;
    cin >> M;
    cin >> J; 
    cout << M <<" " << J << endl;
    Solution * a = new Solution(M,J);
    a->isViable();  
}