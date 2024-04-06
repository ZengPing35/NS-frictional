function r=Gauss_quadrature_2D_test_time(coefficient_function_name,Gauss_coefficient_local,Gauss_point_local,vertices,test_basis_type,test_basis_index,test_derivative_degree_x,test_derivative_degree_y)
%a coefficient function * a test FE basis function (or its derivatives).
%coefficient_function_name: the coefficient function of the integrand.
%Gauss_coefficient_local,Gauss_point_local:the Gauss coefficients and Gauss points on the triangular element T.
%vertices: the coordinates of all vertices of the triangular element T.
%test_basis_type:the type of the test FE basis function.
%test_basis_type=1:2D linear FE.  
%test_basis_type=11:2D P1b FE.
%test_basis_index: the index of test FE basis function to specify which test FE basis function we want to use.
%test_derivative_degree_x:the derivative degree of the test FE basis function with respect to x.
%test_derivative_degree_y:the derivative degree of the test FE basis function with respect to y.
%Gpn: the number of the Gauss points of the Gauss quadrature we are using.

Gpn=length(Gauss_coefficient_local);
r=0;
for i=1:Gpn
      r=r+Gauss_coefficient_local(i)*feval(coefficient_function_name,Gauss_point_local(i,1),Gauss_point_local(i,2))*triangular_local_basis(Gauss_point_local(i,1),Gauss_point_local(i,2),vertices,test_basis_type,test_basis_index,test_derivative_degree_x,test_derivative_degree_y);
end