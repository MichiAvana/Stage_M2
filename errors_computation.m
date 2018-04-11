clear all
close all

%k in [0,2*pi]
k = [-0*pi:0.1*pi:2*pi];
K = length(k);

%N number of points
N = 100;
ind=[]
indice = [100:100:1000];

for errorB=1:10    
    dx = 1/(N);
    ind = [ind dx];

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
    %figure()
    %plot(k,sort((EspecF(1:nv,:)))),title('finite element')


    %exact solution
    if mod(nv,2) == 0
        Niter = (nv/2)-1 ;
    else Niter = ((nv-1)/2);
    end

    for n=-Niter:(Niter+1)
        n+Niter+1;
        E(n+Niter+1,:)= (k.^2 + 4.0*n^2*pi^2 - 4.0*k*n*pi)/2.0;
    end 


    %finite differences 
    EspecL = sparse(N,K);
    for j=1:K
        tL0 =  ((1/dx^2) + (1/2)*k(j)^2 )*ones(1,N);
        tL1 = ((-i*k(j)*(1/(2*dx))-(1/(2*dx^2))))*ones(1,N-1);
        tL2 = ((i*k(j)*(1/(2*dx))-(1/(2*dx^2))))*ones(1,N-1);

        ML = sparse(N,N);
        ML = diag(tL0) + diag(tL1,1) + diag(tL2,-1);
        ML(1,N) = ((i*k(j)*(1/(2*dx))-(1/(2*dx^2))));
        ML(N,1) = ((-i*k(j)*(1/(2*dx))-(1/(2*dx^2))));

        EspecL(:,j) = eig(ML);
    end    
    
    EspecL = sort(EspecL);
    E = sort(E);
    EspecF = sort(EspecF);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCUL DES ERREURS 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %erreurs pour tout k=[0, 2pi] pour tout N avec nv=1er
    err1 = abs((E(1,:) - EspecF(1,:)) / EspecF(1,:));
    err2 = abs((E(1,:) - EspecL(1,:)) / EspecL(1,:)) ;
    diffEFz(errorB) = sqrt(sum(err1.*err1));
    diffEDz(errorB) = sqrt(sum(err2.*err2)); 
    
    
    N= N+100
end




%{
%plot(k,EspecL(1:nv,:)),title('finite differences solution');
EspecL = sort(EspecL);
E = sort(E);
EspecF = sort(EspecF);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCUL DES ERREURS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% POUR K=[0,2PI] ET 1ER VALEUR PROPRES %%%%%
errRel = abs( (EspecF(5,:) - E(5,:) )./ E(5,:))
plot(k,E(5,:), k, EspecF(5,:)),title('E_5, N=1000')
legend('Exact sol','Elements finis','Location','northeast')
xlabel('k')
ylabel('E_5')


figure()
plot(k,errRel),title('E_5, N=1000')
xlabel('k')
ylabel('errors')


%%%%% pour tout k=[0, 2pi] pour tout N avec nv=1er %%%%%%
figure()
plot(ind,diffEDz/K,ind,diffEFz/K), title('k=0 , nev = 1st')
legend('error Exact/Difference F','error Exact/Finite E','Location','northeast')
xlabel('dx')
ylabel('errors')
%}







