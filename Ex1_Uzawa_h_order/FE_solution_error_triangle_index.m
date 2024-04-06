function result=FE_solution_error_triangle_index(uh,uh_exact,left,right,bottom,top,h_partition,basis_index,basis_type,derivative_degree_x,derivative_degree_y,Gauss_point_number,N2_partition_fine)
%Numerically compute a norm error of FE solution on the whole domain [left,right]*[bottom,top].
%uh: the values of the FE solution at all nodes of FE in the whole domain. These values are in 1D index of nodes of FE.
%accurate_function: the accurate function in the error.
%h_partition: the stepsize of the partition.
%basis_type: the type of the FE.
%basis_type=0:2D constant FE.
%basis_type=1:2D Lagrange linear FE.
%basis_type=11:2D P1b FE.
%derivative_degree_x:the derivative degree of the FE solution with respect to x.
%derivative_degree_y:the derivative degree of the FE solution with respect to y.
%Gauss_point_number: the number of the Gauss points of the Gauss quadrature we want to use.
%N1_partition,N2_partition:The N1 and N2 for the partition,not the FE basis functions.
%N1 is the number of the sub-intervals of the partition in x-direction.
%N2 is the number of the sub-intervals of the partition in y-direction.
%vertices: the coordinates of the vertices of a triangular element.
%uh_local: the values of the FE solution at the nodes of FE in a triangular
%element.
%% Parameters
N1_partition=(right-left)/h_partition(1);
N2_partition=(top-bottom)/h_partition(2);
number_of_elements=2*N1_partition*N2_partition;
basis_type = 1;
basis_index = [1 2 3];  

[P_partition,T_partition]=generate_P_T_triangle(left,right,bottom,top,h_partition,1);
P_basis=P_partition;
T_basis=T_partition;

[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(Gauss_point_number);

%Go through all elements and accumulate the error on them.
result1=0;
result2=0;
h = N2_partition_fine/N2_partition;
num = 3;
nodes_fine_mesh = zeros(num, 1);

for n=1:number_of_elements
    vertices=P_partition(:,T_partition(:,n));
    nodes = T_partition(:,n);
    uh_local=uh(nodes);
    for k = 1: num
        i = fix(nodes(k)/(N2_partition+1));
        j = rem(nodes(k), (N2_partition+1));
        if j == 0
            m_fine_mesh = (N2_partition_fine + 1) * (h * (i - 1) + 1);
        else 
            m_fine_mesh = (N2_partition_fine + 1) * h * (i) + (j - 1) * h + 1;
        end
        nodes_fine_mesh(k) = m_fine_mesh;
    end
    uh_local_exact=uh_exact(nodes_fine_mesh);
    result1=result1+Gauss_quad_2D_error_triangle_time(uh_local,uh_local_exact,vertices,basis_index,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type,derivative_degree_x,derivative_degree_y,1);
    result2=result2+Gauss_quad_2D_error_triangle_time(uh_local,uh_local_exact,vertices,basis_index,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type,derivative_degree_x,derivative_degree_y,0);
end
result1=sqrt(result1);
result2=sqrt(result2);
result=result1/result2;