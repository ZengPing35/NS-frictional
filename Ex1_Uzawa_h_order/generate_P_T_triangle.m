function [P,T]=generate_P_T_triangle(left,right,bottom,top,h_partition,basis_type)
%h_partition: is the step size of the partition.
%basis_type: the type of the FE.
%basis_type_p=1:2D Lagrange linear FE.
%basis_type_p=11:2D P1b FE.
%basis_type_u=2:2D Lagrange quadratic FE.
%P stores the coordinates of all nodes for the type of finite element specified by "basis_type".
%P(i,j) is the ith coordinate of the jth node.
%T stores the global indices of the nodes of every element for the type of finite element specified by "basis_type".
%T(i,j) stores the global index of the ith node in th jth element.
%T_partition: store the global indices of the grid points of every element for the partition,not FE.
%P_basis: store the coordinates of all the nodes for the FE,not the partition.
%T_basis: store the global indices of the nodes of every element for FE,not the partition.
%N1 is the number of the sub-intervals of the partition in x-direction.
%N2 is the number of the sub-intervals of the partition in y-direction.
%tnp:total number of all the nodes of FE,including inner nodes and boundary nodes.

h=h_partition;
N1=(right-left)/h(1);
N2=(top-bottom)/h(2);
number_of_elements=2*N1*N2;

if basis_type==0
    P=0;
    T=zeros(1,number_of_elements);
    for n=1:number_of_elements
        T(1,n)=n;
    end
    
elseif basis_type==1

   tnp=(N1+1)*(N2+1);
   P=zeros(2,tnp);
   T=zeros(3,number_of_elements);
   Q=zeros(N1+1,N2+1);

   for j=1:tnp
      if mod(j,N2+1)==0
         P(1,j)=left+(j/(N2+1)-1)*h(1);
         P(2,j)=top;
      else
         P(1,j)=left+fix(j/(N2+1))*h(1);
         P(2,j)=bottom+(mod(j,N2+1)-1)*h(2);
      end
   end

   for i=1:N1+1
      for j=1:N2+1
         Q(i,j)=(i-1)*(N2+1)+j;
      end
   end

%Go through all rectangles in the partition. 
%For the nth rectangle, store the information of its two triangular elements whose element indices are 2n-1 and 2n.
   for n=1:N1*N2
   
      if mod(n,N2)==0
         row=N2;
         column=n/N2;
      else
         row=mod(n,N2);
         column=fix(n/N2)+1;
      end
   
      T(1,2*n-1)=Q(column,row);
      T(2,2*n-1)=Q(column+1,row);
      T(3,2*n-1)=Q(column,row+1);  
  
      T(1,2*n)=Q(column,row+1);
      T(2,2*n)=Q(column+1,row);
      T(3,2*n)=Q(column+1,row+1);  
    
   end

elseif basis_type==11
   tnp1=(N1+1)*(N2+1);
   tnp=(N1+1)*(N2+1)+2*N1*N2;
   P=zeros(2,tnp);
   T=zeros(4,number_of_elements);
   Q=zeros(N1+1,N2+1);

   for j=1:tnp1
      if mod(j,N2+1)==0
         P(1,j)=left+(j/(N2+1)-1)*h(1);
         P(2,j)=top;
      else
         P(1,j)=left+fix(j/(N2+1))*h(1);
         P(2,j)=bottom+(mod(j,N2+1)-1)*h(2);
      end
   end
   
   for i=1:N1+1
      for j=1:N2+1
         Q(i,j)=(i-1)*(N2+1)+j;
      end
   end

%Go through all rectangles in the partition. 
%For the nth rectangle, store the information of its two triangular elements whose element indices are 2n-1 and 2n.
   for n=1:N1*N2
   
      if mod(n,N2)==0
         row=N2;
         column=n/N2;
      else
         row=mod(n,N2);
         column=fix(n/N2)+1;
      end
   
      T(1,2*n-1)=Q(column,row);
      T(2,2*n-1)=Q(column+1,row);
      T(3,2*n-1)=Q(column,row+1);  
      T(4,2*n-1)=tnp1+1+2*(n-1);
  
      T(1,2*n)=Q(column,row+1);
      T(2,2*n)=Q(column+1,row);
      T(3,2*n)=Q(column+1,row+1);  
      T(4,2*n)=tnp1+2*n;
    
   end   
   for j=tnp1+1:tnp
       P(1,j)=(P(1,T(1,j-tnp1))+P(1,T(2,j-tnp1))+P(1,T(3,j-tnp1)))/3;
       P(2,j)=(P(2,T(1,j-tnp1))+P(2,T(2,j-tnp1))+P(2,T(3,j-tnp1)))/3;
   end
   
   
   
elseif basis_type==2

   dh=h/2;
   dN1=N1*2;
   dN2=N2*2;
   tnp=(dN1+1)*(dN2+1);
   P=zeros(2,tnp);
   T=zeros(6,number_of_elements);
   Q=zeros(dN1+1,dN2+1);

   for j=1:tnp
      if mod(j,dN2+1)==0
         P(1,j)=left+(j/(dN2+1)-1)*dh(1);
         P(2,j)=top;
      else
         P(1,j)=left+fix(j/(dN2+1))*dh(1);
         P(2,j)=bottom+(mod(j,dN2+1)-1)*dh(2);
      end
   end

   for i=1:dN1+1
      for j=1:dN2+1
         Q(i,j)=(i-1)*(dN2+1)+j;
      end
   end

%Go through all rectangles in the partition. 
%For the nth rectangle, store the information of its two triangular elements whose element indices are 2n-1 and 2n.
   for n=1:N1*N2
   
      if mod(n,N2)==0
         row=N2;
         column=n/N2;
      else
         row=mod(n,N2);
         column=fix(n/N2)+1;
      end
   
      T(1,2*n-1)=Q(2*column-1,2*row-1);
      T(2,2*n-1)=Q(2*column+1,2*row-1); 
      T(3,2*n-1)=Q(2*column-1,2*row+1);
      T(4,2*n-1)=Q(2*column,2*row-1);
      T(5,2*n-1)=Q(2*column,2*row);
      T(6,2*n-1)=Q(2*column-1,2*row);


      T(1,2*n)=Q(2*column-1,2*row+1);
      T(2,2*n)=Q(2*column+1,2*row-1);
      T(3,2*n)=Q(2*column+1,2*row+1);
      T(4,2*n)=Q(2*column,2*row);
      T(5,2*n)=Q(2*column+1,2*row);
      T(6,2*n)=Q(2*column,2*row+1); 

   end    
end
