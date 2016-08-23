//Code containing formula for the most significant coefficients

using namespace std;

void do_Terms( vector<FunctionBase*> &functions, bool debug = false){

   if(debug) cout << "Start adding terms" << endl;

      Function *fcph1a = new Function("fcph1a","cos(1*z)*cos(1*t)",4.,4.); functions.push_back(fcph1a);

   if(debug) cout << "Done adding terms" << endl;

}

