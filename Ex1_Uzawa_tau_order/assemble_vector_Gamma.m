function b = assemble_vector_Gamma(h_partition, lamda_old, t,vector_size, boundary_nodes, P_partition, T_partition, T_basis_test, Gauss_coefficient_reference_triangle, Gauss_point_reference_triangle, test_basis_index, test_basis_type, test_derivative_degree_x, test_derivative_degree_y)
%% Initialization
b = sparse(vector_size,1);
[s,p] = find(boundary_nodes(1,:)==-2);
nbn = size(p,2);
lamda = zeros(1,nbn);
lamda(1,:) = boundary_nodes(3,p);
D = h_partition(1) * h_partition(2);

for k = 1: nbn
    node = lamda(1,k);
    x = P_partition(1,node);
    y = P_partition(2,node);
    g = function_g(x,y);
    [u,v] = find(lamda_old(1,:) == node);
    L = lamda_old(2,v);
    results = (D * L * g);
    b(node,1)=b(node, 1)+results;
end




% %% Initialization
% b = sparse(vector_size,1);
% [s,p] = find(boundary_nodes(1,:)==-2);
% nbn = size(p,2);
% lamda = zeros(1,nbn);
% lamda(1,:) = boundary_nodes(3,p);
% test_basis_index_length=length(test_basis_index);
% D = h_partition(1) * h_partition(2);
% g = function_g(t);
% 
% for k = 1: nbn
%     node = lamda(1,k);
%     [u,v] = find(lamda_old(1,:) == node);
%     L = lamda_old(2,v);
%     [s,m] = find(T_partition == node);
%     num_m = size(m,1);
%     for j = 1: num_m
%         vertices=P_partition(:,T_partition(:, m(j)));
%         [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);
% 
%         for i=1:test_basis_index_length   
%             beta=test_basis_index(i);    
% %             temp=Gauss_quadrature_2D_test_Gamma(Gauss_coefficient_local_triangle,Gauss_point_local_triangle,vertices,test_basis_type,beta,test_derivative_degree_x,test_derivative_degree_y);
%             temp=1;
%             r = T_basis_test(beta,m(j));
%             if r==node
%                 results = (D * L * g * temp);
%                 b(r,1)=b(r, 1)+results;
%             end
%         end
%     end
% end
