% % % %
%CALCUL DES VALEURS PROPRES
% UTILISANT L'EXPRESSION ANALYTIQUE
% FORMULATION EXACTE
% % % % 

%Bandes d'energies sans potentiel 
%analytiquement
clear all
k = [-0*pi:0.1*pi:2*pi];
N=2
%Nbp=10
%N = ((Nbp)/2);
K = length(k);
E = zeros(N,K);
l = 6


% dx = 1/(Nbp-1);
% %x varie entre [0,1]
% x = [0:dx:1];
% 
% %definition du potentiel 
% %V = 10*cos(2*(2*pi*x)) + 30*tanh(x).*cos(2*pi*x);
% V = 10*cos(2*(2*pi*x));


for n=-N:(N+1)
    indice = n+N+1
    E(n+N+1,:)= (k.^2 + 4.0*n^2*pi^2 - 4.0*k*n*pi )/2.0;
    %E(n+N+1,:)= (k.^2 + 4.0*n^2*pi^2 - 4.0*k*n*pi+ 2*V(n+N+1))/2.0;
end 

figure()
plot(k,E);,title('solution exact')
xlabel('k');
ylabel('En(k)');

%sort(E(1:5,1))




