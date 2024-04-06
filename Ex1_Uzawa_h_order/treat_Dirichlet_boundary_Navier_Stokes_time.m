function [A,b]=treat_Dirichlet_boundary_Navier_Stokes_time(Dirichlet_boundary_function_name_u1,Dirichlet_boundary_function_name_u2,Dirichlet_boundary_function_name_u1_Gamma,A,b,boundary_nodes,P_basis_u,number_of_FE_nodes_u)
%Deal with Dirichlet boundary nodes.
%boundary_nodes(1,k): specifiy the type of the kth boundary node for the normal direction(or u1).
%boundary_nodes(1,k)=-1: Dirichlet boundary node in the normal direction(or u1);
%boundary_nodes(1,k)=-2: The nodes on Gamma in the normal direction(or u1);
%boundary_nodes(1,k)=-3: Robin boundary node in the normal direction(or u1). 
%boundary_nodes(2,k): specifiy the type of the kth boundary node for the tangential direction(or u2).
%boundary_nodes(2,k)=-1: Dirichlet boundary node in the tangential direction(or u2);
%boundary_nodes(2,k)=-2: The nodes on Gamma in the tangential direction(or u2);
%boundary_nodes(2,k)=-3: Robin boundary node in the tangential direction(or u2).
%The intersection node between Dirichlet boundary and other boundaries is a Dirichlet boundary node.
%boundary_nodes(3,k): global index of the kth boundary node among all nodes of FE. 
%                     That is, the index of FE is used here.
%number_of_FE_nodes_u: the number of the FE nodes for u. 
%nbn: the total number of all the boundary nodes of FE.

nbn=size(boundary_nodes,2);
%Check all boundary nodes of FE.
for k=1:nbn

    if boundary_nodes(1,k)==-1 
        i=boundary_nodes(3,k);
        A(i,:)=0;
        A(i,i)=1;
        b(i,1)=feval(Dirichlet_boundary_function_name_u1,P_basis_u(1,i),P_basis_u(2,i));
    end    
    if boundary_nodes(2,k)==-1 
        i=boundary_nodes(3,k);
        u2_index=number_of_FE_nodes_u+i;
        A(u2_index,:)=0;
        A(u2_index,u2_index)=1;
        b(u2_index,1)=feval(Dirichlet_boundary_function_name_u2,P_basis_u(1,i),P_basis_u(2,i));
    end    
    
    
    if boundary_nodes(1,k)==-2 
        i=boundary_nodes(3,k);
        A(i,:)=0;
        A(i,i)=1;
        b(i,1)=feval(Dirichlet_boundary_function_name_u1_Gamma,P_basis_u(1,i),P_basis_u(2,i));
    end  
    
%     if boundary_nodes(2,k)==-2
%         i=boundary_nodes(3,k);
%         u2_index=number_of_FE_nodes_u+i;
%         A(u2_index,:)=0;
%         A(u2_index,u2_index)=1;
%         b(u2_index,1)=feval(Dirichlet_boundary_function_name_u2_Gamma,P_basis_u(1,i),P_basis_u(2,i));
%     end
end

