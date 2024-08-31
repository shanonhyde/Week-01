#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

int main()
{
	double test1, test2, test3, test4, test5, average;
	
	cout << "Enter five test scores: ";
    cin >> test1;
    cin >> test2;
    cin >> test3;
    cin >> test4;
    cin >> test5;

    average = (test1+test2+test3+test4+test5)/5.0;
	cout<<"Average=" << average;
	cout<<showpoint<<setprecision(1)<<fixed;
	cout << "\nAverage: " << average << endl;
	return 0;
	
}