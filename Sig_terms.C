//Code containing formula for the most significant coefficients

using namespace std;

void do_Terms( vector<FunctionBase*> &functions, bool debug = false){

  Function *fdph1 = new Function("fdph1","cos(1*z)",2.,2.,"Phi"); functions.push_back(fdph1);
  Function *fdph2 = new Function("fdph2","cos(2*z)",2.,2.,"Phi"); functions.push_back(fdph2);
  Function *fdph3 = new Function("fdph3","cos(3*z)",2.,2.,"Phi"); functions.push_back(fdph3);

  Function *fsph1 = new Function("fsph1","cos(1*t)",2.,2.,"Phi"); functions.push_back(fsph1);
  Function *fsph2 = new Function("fsph2","cos(2*t)",2.,2.,"Phi"); functions.push_back(fsph2);
  Function *fsph3 = new Function("fsph3","cos(3*t)",2.,2.,"Phi"); functions.push_back(fsph3);

  //if(debug){}

}

