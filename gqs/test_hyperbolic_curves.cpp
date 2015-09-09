

#include "gq.hpp"


double tol = 1e-6;


#define FLAG_ERROR throw 1

inline void flag_error() { FLAG_ERROR; }

#define EQUAL_TEST( TEST, TYPE ) if( !(TEST) ) {	\
    printf( "Equality Test Failed at line %d: \n", line);	\
  printf( "  Expected value: %" #TYPE "\n", A ); \
  printf( "  Actual value:   %" #TYPE "\n", B ); \
  printf( "\n" ); \
  flag_error();	  \
  }

#define CHECK_EQUAL( EXP, ACT ) check_equal ( (EXP), (ACT), __LINE__ ) 
#define CHECK_REAL_EQUAL( EXP, ACT, EPS ) check_equal ( (EXP), (ACT), (EPS), __LINE__ )
#define CHECK_CUBITVECTORS_EQUAL( EXP, ACT ) check_CubitVectors_equal( (EXP), (ACT), #EXP, #ACT, __LINE__ )

bool check_equal(int A, int B, int line) { EQUAL_TEST( A == B, d) };
bool check_equal(double A, double B, int line) { EQUAL_TEST( A == B, f) };
bool check_equal(double A, double B, double eps, int line) { EQUAL_TEST( (A-B < eps), f ) };
bool check_CubitVectors_equal( CubitVector A, CubitVector B, const char* Aname, const char* Bname, int line ) 
{
if ( (A[0] - B[0] > tol) || ( A[1] - B[1] > tol ) || (A[2] - B[2] > tol ) )
  {

std::cout << "=====" << std::endl << "ERROR" << std::endl << "=====" << std::endl;
std::cout << "CubitVector " << Aname << " does not equal " << Bname
          << " (from line " << line << ")" << std::endl;
CHECK_REAL_EQUAL( A[0], B[0], tol);
CHECK_REAL_EQUAL( A[1], B[1], tol);
CHECK_REAL_EQUAL( A[2], B[2], tol);

}
}



int main()
{

  
  //hyperbolic function coefficients
  double A = 1;
  double B = -2;

  DLIList<RefEdge*> edges;
  //generate the desired curves
  hyperbolic_curves(A,B,edges);

  
  BasicTopologyEntity* bte1 = dynamic_cast<BasicTopologyEntity*>(edges[0]);
  BasicTopologyEntity* bte2 = dynamic_cast<BasicTopologyEntity*>(edges[1]);
    
  CubitBox edge1_box = bte1->bounding_box();
  CubitBox edge2_box = bte2->bounding_box();

  //check the bounds of each curve

  ////// Curve 1 Check ///////
  CubitVector expected_startpoint(1.4142135,2.0,0);
  CubitVector actual_startpoint = edges[0]->start_vertex()->center_point();
  CHECK_CUBITVECTORS_EQUAL(expected_startpoint, actual_startpoint);
  
  CubitVector expected_midpoint(1.0,0.0,0.0);
  CubitVector actual_midpoint;
  edges[0]->mid_point(actual_midpoint);
  CHECK_CUBITVECTORS_EQUAL(expected_midpoint, actual_midpoint );

  CubitVector expected_endpoint(1.4142135,-2.0,0);
  CubitVector actual_endpoint = edges[0]->end_vertex()->center_point();
  CHECK_CUBITVECTORS_EQUAL(expected_endpoint, actual_endpoint);

  ////// Curve 2 Check ///////
  expected_startpoint.set(-1.4142135,2.0,0);
  actual_startpoint = edges[1]->start_vertex()->center_point();
  CHECK_CUBITVECTORS_EQUAL(expected_startpoint, actual_startpoint);
  
  expected_midpoint.set(1.0,0.0,0.0);
  actual_midpoint;
  edges[0]->mid_point(actual_midpoint);
  CHECK_CUBITVECTORS_EQUAL(expected_midpoint, actual_midpoint );

  expected_endpoint.set(-1.4142135,-2.0,0);
  actual_endpoint = edges[1]->end_vertex()->center_point();
  CHECK_CUBITVECTORS_EQUAL(expected_endpoint, actual_endpoint);

  
  
  if (false)
    {
      DLIList<RefEntity*> exp_bodies;
      int exp_ents;
      CubitString cubit_version("12.2");
      CubitCompat_export_solid_model(exp_bodies, "test_file.sat", "ACIS_SAT", exp_ents, cubit_version);
    }



}

