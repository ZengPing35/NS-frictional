function Lamda = function_update_Lamda(P_partition, lamda_old, u2h_old, rho, t, N1_partition, boundary_nodes)

%% Initialization
vector_size = N1_partition - 1;
Lamda = zeros(2, vector_size);
[s,p] = find(boundary_nodes(1,:)==-2);
Lamda(1,:) = boundary_nodes(3,p);
g = function_g(t);

for k = 1: vector_size
    node = Lamda(1,k);
%     x = P_partition(1,node);
%     y = P_partition(2,node);
%     g = function_g(x,y);
    [u,v] = find(lamda_old(1,:) == node);
    L = lamda_old(2,v);
    uh2 = u2h_old(node);
    Temp = L - rho * g * uh2;
    Lamda(2,k) = max(-1, min(1,Temp));
end