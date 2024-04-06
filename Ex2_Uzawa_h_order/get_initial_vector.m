function r=get_initial_vector(initial_function_name,P_basis)
%Evaluate the initial function at all nodes.
%initial_function_name: the name of the initial funtion.

number_of_nodes=size(P_basis,2);
r=zeros(number_of_nodes,1);
for i=1:number_of_nodes
    r(i)=feval(initial_function_name,P_basis(1,i),P_basis(2,i));
end