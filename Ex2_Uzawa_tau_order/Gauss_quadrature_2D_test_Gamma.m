function r=Gauss_quadrature_2D_test_Gamma(Gauss_coefficient_local,Gauss_point_local,vertices,test_basis_type,test_basis_index,test_derivative_degree_x,test_derivative_degree_y)

Gpn=length(Gauss_coefficient_local);
r=0;
for i=1:Gpn
      r=r + Gauss_coefficient_local(i)*triangular_local_basis(Gauss_point_local(i,1),Gauss_point_local(i,2),vertices,test_basis_type,test_basis_index,test_derivative_degree_x,test_derivative_degree_y);
end