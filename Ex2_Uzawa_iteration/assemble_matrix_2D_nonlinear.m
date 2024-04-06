function result=assemble_matrix_2D_nonlinear(coefficient_function_name,uh,basis_index_coe,basis_type_coe,derivative_degree_x_coe,derivative_degree_y_coe,P_partition,T_partition,T_basis_trial,T_basis_test,number_of_elements,matrix_size,trial_basis_index,test_basis_index,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,trial_basis_index_length,trial_basis_type,trial_derivative_degree_x,trial_derivative_degree_y,test_basis_index_length,test_basis_type,test_derivative_degree_x,test_derivative_degree_y)
%Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle: the Gauss coefficients and Gauss points on the reference triangular element.
%vertices: the coordinates of all vertices of a triangular element.
%Gauss_coefficient_local_triangle,Gauss_point_local_triangle: the Gauss coefficients and Gauss points on the local triangular element.

result=sparse(matrix_size(1),matrix_size(2));
for n=1:number_of_elements
   
    vertices=P_partition(:,T_partition(:,n));
    [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);
    uh_local=uh(T_basis_trial(:,n));
    
    for j=1:trial_basis_index_length
        alpha=trial_basis_index(j);
        for i=1:test_basis_index_length   
            beta=test_basis_index(i);   
            temp=Gauss_quadrature_2D_trial_test_triangle_nonlinear(coefficient_function_name,Gauss_coefficient_local_triangle,Gauss_point_local_triangle,uh_local,vertices,basis_index_coe,basis_type_coe,derivative_degree_x_coe,derivative_degree_y_coe,trial_basis_type,alpha,trial_derivative_degree_x,trial_derivative_degree_y,test_basis_type,beta,test_derivative_degree_x,test_derivative_degree_y);
            result(T_basis_test(beta,n),T_basis_trial(alpha,n))=result(T_basis_test(beta,n),T_basis_trial(alpha,n))+temp;
       end
    end

end