%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATION OF THE SOLUTION 
% USING FINITE ELEMENTS
% WITH A POTENTIAL 10*cos(4*pi*x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%k in [0,2*pi]
k = [-0*pi:0.1*pi:2*pi];
K = length(k);


%N+1 number of points
N = 1000;

%for errorB=1:10    
dx = 1/(N-1);

%number of eigen values
nv= 6;

%table of degree of freedom 
ddl=zeros(N,2);
vec1 = [1:1:N];
vec2 = [2:1:N,1];
ddl(:,1) = vec1;
ddl(:,2) = vec2;

%x varie entre [0,1] (N+1) point
dx1 = 1/(N);
x = [0:dx1:1];
X = length(x);

 A = [ dx/3 dx/6 ; dx/6 dx/3];
 B = [-1/2 -1/2 ; 1/2 1/2];
 C = [1/dx -1/dx ; -1/dx  1/dx];
        
    EspecF = sparse(nv,K);
    for j=1:K
        M = zeros(N,N);
        ML= zeros(N,N);

        for e=1:N
            ind = e;
            sini = sin(4*pi*x(ind));
            sinip = sin(4*pi*x(ind+1));
            cosi = cos(4*pi*x(ind));
            cosip = cos(4*pi*x(ind+1));
            
%             Pii = sin(4*pi*x(ind))*(-(dx^2/(4*pi)) - 1/(32*pi^3))  + (dx/(8*pi^2))*cos(4*pi*x(ind)) + (1/(32*pi^3))*sin(4*pi*x(ind+1)) ;
% 
%             Pipip = sin(4*pi*x(ind+1))*((dx^2/(4*pi)) - 1/(32*pi^3)) + (dx/(8*pi^2))*cos(4*pi*x(ind+1)) + (1/(32*pi^3))*sin(4*pi*x(ind)) ;
%             
%             Piip = (-dx/(16*pi^2))*(cos(4*pi*x(ind+1)) + cos(4*pi*x(ind))) + (dx/(32*pi^3))*(sin(4*pi*x(ind+1)) - sin(4*pi*x(ind)));
%             
            Pii = sini*(-(dx^2/(4*pi)) + 1/(32*pi^3))  + (dx/(8*pi^2))*cosi - (1/(32*pi^3))*sinip ;

            Pipip = sinip*((dx^2/(4*pi)) - 1/(32*pi^3)) + (dx/(8*pi^2))*cosip + (1/(32*pi^3))*sini ;
            
            Piip = (-dx/(16*pi^2))*(cosip + cosi) + (1/(32*pi^3))*(sinip - sini);          
            
            P = (10/(dx^2))*[Pii Piip ; Piip Pipip];

            %Laplace matrix
            D = (k(j)^2/2)*A - i*k(j)*B + (1/2)*C + P;
            
            for l=1:2
                for c=1:2
                    M(ddl(e,l),ddl(e,c)) = M(ddl(e,l),ddl(e,c))+D(l,c);
                    ML(ddl(e,l),ddl(e,c)) = ML(ddl(e,l),ddl(e,c))+A(l,c);
                end
            end
        end 
        
        
        EspecF(1:nv,j) = eigs(M,ML,nv,'SM');  
    end
%plot(k,EspecF(1:nv,:)),title('finite differences solution');