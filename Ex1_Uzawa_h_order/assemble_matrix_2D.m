function result=assemble_matrix_2D(coefficient_function_name,P_partition,T_partition,T_basis_trial,T_basis_test,number_of_elements,matrix_size,trial_basis_index,test_basis_index,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,trial_basis_index_length,trial_basis_type,trial_derivative_degree_x,trial_derivative_degree_y,test_basis_index_length,test_basis_type,test_derivative_degree_x,test_derivative_degree_y)

result=sparse(matrix_size(1),matrix_size(2));
% trial_basis_index_length=length(trial_basis_index);
% test_basis_index_length=length(test_basis_index);
%Go through all elements.
%On each element, compute the volume integrals for all possible combinations of trial and test FE basis functions.
%Assemble the values of those volume integrals into the matrix.
for n=1:number_of_elements
   
    vertices=P_partition(:,T_partition(:,n));
    [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);
   
    for j=1:trial_basis_index_length
        alpha=trial_basis_index(j);
        for i=1:test_basis_index_length   
            beta=test_basis_index(i);   
            temp=Gauss_quadrature_2D_trial_test_triangle(coefficient_function_name,Gauss_coefficient_local_triangle,Gauss_point_local_triangle,vertices,trial_basis_type,alpha,trial_derivative_degree_x,trial_derivative_degree_y,test_basis_type,beta,test_derivative_degree_x,test_derivative_degree_y);
            result(T_basis_test(beta,n),T_basis_trial(alpha,n))=result(T_basis_test(beta,n),T_basis_trial(alpha,n))+temp;
       end
    end

end