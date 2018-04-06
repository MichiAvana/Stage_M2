% % % %
%CALCUL DES VALEURS PROPRES
% UTILISANT LES ELEMENTS FINIS 
% % % % 

clear all
close all

%k in [0,2*pi]
k = [-0*pi:0.1*pi:2*pi];
K = length(k);

%N number of points
N = 1000;

%number of eigen values
nv= 6;

%table of degree of freedom 
ddl=zeros(N,2);
vec1 = [1:1:N];
vec2 = [2:1:N,1];
ddl(:,1) = vec1;
ddl(:,2) = vec2;

A = [ dx/3 dx/6 ; dx/6 dx/3];
B = [-1/2 -1/2 ; 1/2 1/2];
C = [1/dx -1/dx ; -1/dx  1/dx];

EspecF = sparse(nv,K);
   for j=1:K
    M = zeros(N,N);
    ML= zeros(N,N);

    %Laplace matrix
    D = (k(j)^2/2)*A - i*k(j)*B + (1/2)*C;

    for e=1:N
        for l=1:2
            for c=1:2
                M(ddl(e,l),ddl(e,c)) = M(ddl(e,l),ddl(e,c))+D(l,c);
                ML(ddl(e,l),ddl(e,c)) = ML(ddl(e,l),ddl(e,c))+A(l,c);
            end
        end
    end 

    EspecF(1:nv,j) = eigs(M,ML,nv,'SM');  
end
figure()
plot(k,sort((EspecF(1:nv,:)))),title('finite element')

