//Code containing formula for the most significant coefficients

using namespace std;

void do_Terms( vector<FunctionBase*> &functions){

      Function *fcph1a = new Function("fcph1a","cos(1*z)*cos(1*t)",4.,4.,"Phi"); functions.push_back(fcph1a);
      Function *fcph1b = new Function("fcph1b","sin(1*z)*sin(1*t)",4.,4.,"Phi"); functions.push_back(fcph1b);
      Function *ft22 = new Function("ft22","(cos(2*x)*cos(2*y))",4.,4.,"Theta"); functions.push_back(ft22);
      Function *ft24 = new Function("ft24","(cos(2*x)*cos(4*y) + cos(4*x)*cos(2*y))",2.,4.,"Theta"); functions.push_back(ft24);

      Function *fdp1 = new Function("fdp1","cos(1*([0]))",2.,2.,"Mom"); functions.push_back(fdp1);
      Function *fdp2 = new Function("fdp2","cos(2*([0]))",2.,2.,"Mom"); functions.push_back(fdp2);
      Function *fdp3 = new Function("fdp3","cos(3*([0]))",2.,2.,"Mom"); functions.push_back(fdp3);
      Function *fdp4 = new Function("fdp4","cos(4*([0]))",2.,2.,"Mom"); functions.push_back(fdp4);
      Function *fdp5 = new Function("fdp5","cos(5*([0]))",2.,2.,"Mom"); functions.push_back(fdp5);
      Function *fdp6 = new Function("fdp6","cos(6*([0]))",2.,2.,"Mom"); functions.push_back(fdp6);

      Function *fsp1 = new Function("fsp1","cos(1*([1]))",2.,2.,"Mom"); functions.push_back(fsp1);
      Function *fsp2 = new Function("fsp2","cos(2*([1]))",2.,2.,"Mom"); functions.push_back(fsp2);
      Function *fsp3 = new Function("fsp3","cos(3*([1]))",2.,2.,"Mom"); functions.push_back(fsp3);

      Function *fdsp11 = new Function("fdsp11","cos(1*([0]))*cos(1*([1]))",4.,4.,"Mom"); functions.push_back(fdsp11);
      Function *fdsp22 = new Function("fdsp22","cos(2*([0]))*cos(2*([1]))",4.,4.,"Mom"); functions.push_back(fdsp22);
      Function *fdsp33 = new Function("fdsp33","cos(3*([0]))*cos(3*([1]))",4.,4.,"Mom"); functions.push_back(fdsp33);

}

