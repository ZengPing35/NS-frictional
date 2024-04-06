function r=Gauss_quadrature_2D_trial_test_triangle_nonlinear(coefficient_function_name,Gauss_coefficient_local,Gauss_point_local,uh_local,vertices,basis_index_coe,basis_type_coe,derivative_degree_x_coe,derivative_degree_y_coe,trial_basis_type,trial_basis_index,trial_derivative_degree_x,trial_derivative_degree_y,test_basis_type,test_basis_index,test_derivative_degree_x,test_derivative_degree_y)
%Gpn: the number of the Gauss points of the Gauss quadrature we are using.

Gpn=length(Gauss_coefficient_local);
r=0;
for i=1:Gpn    
     r=r+Gauss_coefficient_local(i)*feval(coefficient_function_name,Gauss_point_local(i,1),Gauss_point_local(i,2),uh_local,vertices,basis_index_coe,basis_type_coe,derivative_degree_x_coe,derivative_degree_y_coe)*triangular_local_basis(Gauss_point_local(i,1),Gauss_point_local(i,2),vertices,trial_basis_type,trial_basis_index,trial_derivative_degree_x,trial_derivative_degree_y)*triangular_local_basis(Gauss_point_local(i,1),Gauss_point_local(i,2),vertices,test_basis_type,test_basis_index,test_derivative_degree_x,test_derivative_degree_y);
end