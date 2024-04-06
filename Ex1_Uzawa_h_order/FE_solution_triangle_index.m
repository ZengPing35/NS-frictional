function result=FE_solution_triangle_index(x,y,uh_local,vertices,basis_index,basis_type,derivative_degree_x,derivative_degree_y)
%Evaluate a finite element solution at (x,y) which is in a triangular element T.
%uh_local: the values of numerical solution at all the nodes of FE basis functions in the triangular element T.
%vertices: the coordinates of all vertices of the triangular element T.

result=0;
basis_index_length=length(basis_index);

for k=1:basis_index_length   
    i=basis_index(k);
    result=result+uh_local(i)*triangular_local_basis(x,y,vertices,basis_type,i,derivative_degree_x,derivative_degree_y);
end