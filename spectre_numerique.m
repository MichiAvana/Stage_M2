% % % %
%CALCUL DES VALEURS PROPRES
% UTILISANT LES SCHEMAS NUMERIQUES 
% DIFFERENCES FINI
% % % % 

clear all

%pour k varie entre [0;pi]
k = [0*pi:0.1*pi:2*pi];
K = length(k);
N = 1000;
dx = 1/(N-1);


%SANS POTENTIEL
%n est le nombre de valeurs propres
n=6;

%Leap-Frog scheme
EspecL = sparse(N,K);
for j=1:K
    tL0 =  ((1/dx^2) + (1/2)*k(j)^2 )*ones(1,N);
    tL1 = (-(-i*k(j)*(1/(2*dx))-(1/(2*dx^2))))*ones(1,N-1);
    tL2 = (-(i*k(j)*(1/(2*dx))-(1/(2*dx^2))))*ones(1,N-1);

    ML = sparse(N,N);
    ML = diag(tL0) + diag(tL1,1) + diag(tL2,-1);
    ML(1,N) = (-(i*k(j)*(1/(2*dx))-(1/(2*dx^2))));
    ML(N,1) = (-(-i*k(j)*(1/(2*dx))-(1/(2*dx^2))));
    
    EspecL(:,j) = eig(ML);
end


%AVEC POTENTIEL 
%n est le nombre de valeur propres
%n = 10
%central scheme
EspecLP = sparse(N,K);

%x varie entre [0,1]
x = [0:dx:1];

%definition du potentiel 
V = 10*cos(2*(2*pi*x));
%V = 10*cos(2*(2*pi*x)) + 0.5*30*cos(2*pi*x)
%V = 10*cos(2*(2*pi*x)) + 30*tanh(x).*cos(2*pi*x)
%V = 2

for j=1:K
    %construction de la matrice diagonale
    tL0 =  ((1/dx^2) + (1/2)*k(j)^2 +V).*ones(1,N);
    tL1 = (-(-i*k(j)*(1/(2*dx))-(1/(2*dx^2))))*ones(1,N-1);
    tL2 = (-(i*k(j)*(1/(2*dx))-(1/(2*dx^2))))*ones(1,N-1);

    ML = sparse(N,N);
    ML = diag(tL0) + diag(tL1,1) + diag(tL2,-1);
    ML(1,N) = (-(i*k(j)*(1/(2*dx))-(1/(2*dx^2))));
    ML(N,1) = (-(-i*k(j)*(1/(2*dx))-(1/(2*dx^2))));
    
    %calcul des valeurs propres
    EspecLP(:,j) = eig(ML);
end

%afficher les figure 
subplot(2,1,1),plot(k,sort(EspecL(1:n,:))),title('without potential');
subplot(2,1,2),plot(k,sort(EspecLP(1:n,:))),title('with potential V = 10*cos(2*(2*pi*x))}');

xlabel('k')
ylabel('En(k)')
