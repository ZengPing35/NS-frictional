function result=assemble_vector_2D_time(coefficient_function_name,N2,P_partition,T_partition,T_basis_test,number_of_elements,vector_size,test_basis_index,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,test_basis_type,test_derivative_degree_x,test_derivative_degree_y)
%coefficient_function_name: the coefficient function of the integrand.
%P_partition: store the coordinates of all the grid points of the partition,not FE.
%T_partition: store the global indices of the grid points of every element for the partition,not FE.
%T_basis: store the global indices of the nodes of every element for FE,not the partition.
%T_basis_test: T_basis for the test basis function.
%h_partition: the step size of the partition
%number_of_element: the number of the local triangluar elements of the partition.
%test_basis_index: a vector which stores the indices of the local test basis functions we need for the matrix.
%Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle: the Gauss coefficients and Gauss points on the reference triangular element.
%test_basis_type:the type of the test FE basis function.
%test_basis_index: the index of test FE basis function to specify which test FE basis function we want to use.
%test_derivative_degree_x:the derivative degree of the test FE basis function with respect to x.
%test_derivative_degree_y:the derivative degree of the test FE basis function with respect to y.
%vertices: the coordinates of all vertices of a triangular element.
%Gauss_coefficient_local_triangle,Gauss_point_local_triangle: the Gauss coefficients and Gauss points on the local triangular element.

result=sparse(vector_size,1);
test_basis_index_length=length(test_basis_index);
%Go through all elements.
%On each element, compute the volume integrals for all test FE basis functions.
%Assemble the values of those volume integrals into the vector.
for n=1:number_of_elements
    if mod(n,2)==0
        m = n-1;
    else
        m = n;
    end
    if mod(fix(m/N2),2)~=0   %%uper half of the region
        vertices=P_partition(:,T_partition(:,n));
        [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);

        for i=1:test_basis_index_length   
            beta=test_basis_index(i);    
            temp=Gauss_quadrature_2D_test_time(coefficient_function_name,Gauss_coefficient_local_triangle,Gauss_point_local_triangle,vertices,test_basis_type,beta,test_derivative_degree_x,test_derivative_degree_y);
            result(T_basis_test(beta,n),1)=result(T_basis_test(beta,n),1)+temp;
        end
    end
end