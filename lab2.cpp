#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
using namespace std;

int number_of_digits(float number)
{
    int numberInt = (int) number;
    return floor(log10(numberInt)+4);
}

int main()
{
  int sharesBought;
  float pricePerShare;
  int percentCom;
  float annualPercent;
  int years;

  float amountPaid;
  float amountCom;
  float total;
  float sharesWorth;

  string amountPaidOutput;
  string amountComOutput;
  string totalOutput;
  string sharesWorthOutput;

  char receipt;

  ifstream rfile;
  ofstream wfile;

  rfile.open("sample.txt");

  ofstream writeFile("receipt.txt");

if (rfile.is_open()) 
{
    string line;
  
        getline(rfile, line);
        cout << line << endl;
        cin >> sharesBought;

        getline(rfile, line);
        cout << line << endl;
        cin >> pricePerShare;

        getline(rfile, line);
        cout << line << endl;
        cin >> percentCom;

        getline(rfile, line);
        cout << line << endl;
        cin >> annualPercent;

        getline(rfile, line);
        cout << line << endl;
        cin >> years;

        amountPaid = sharesBought * pricePerShare;
        amountCom = amountPaid * ((float)percentCom/100);
        total = amountPaid + amountCom;
        sharesWorth =  amountPaid * pow((1+(annualPercent/100)), years);

        //cout << precision(2) << fixed << amountPaid;
        // to_string << fixed << setprecision(2);
        
        amountPaidOutput = "The amount paid for the stock alone (without the commission): " + to_string(amountPaid);

        amountComOutput = "The amount of the commission: " + to_string(amountCom);

        totalOutput = "The total amount of paid (the payment for stock plus the commission): " + to_string(total);

        sharesWorthOutput = "After X years, your shares will be worth: " + to_string(sharesWorth);

        cout << amountPaidOutput << endl;
        cout << amountComOutput << endl;
        cout << totalOutput << endl;
        cout << sharesWorthOutput << endl;
        cout << "Would you like a receipt? Y/N: " << endl;

        cin >> receipt;

        if(receipt == 'Y')
        {
          writeFile << "Gabrielle Nibert" << endl;
          writeFile << "-------------------------------" << endl;

          int setWidth = number_of_digits(sharesWorth);
        
          writeFile << "Total Stock:          $ " << setprecision(2) << fixed << setw(setWidth) << right << amountPaid << endl;
          writeFile << "Commission:           $ " << setprecision(2) << fixed << setw(setWidth) << right << amountCom << endl;
          writeFile << "Total amount:         $ " << setprecision(2) << fixed << setw(setWidth) << right << total << endl;
          writeFile << "Net worth in X years: $ " << setprecision(2) << fixed << setw(setWidth) << right << sharesWorth << endl;
        } else {
          cout << "Have a nice day" << endl;
        }

        

    writeFile.close();
    rfile.close();
}

  return 0;
}