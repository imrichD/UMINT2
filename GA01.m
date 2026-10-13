% Simple GA
% Ivan Sekaj, 2022
% ==================

numgen=500	% number of generations
lpop=50;	% number of chromosomes in population
lstring=10;	% number of genes in a chromosome
M=500;          % maximum of the search space

figure(1); hold on;

% Population initialisation

Space=[ones(1,lstring)*(-M); ones(1,lstring)*M];
Delta=Space(2,:)/100;    % additive mutation step

Pop=genrpop(lpop,Space);

% Main cyklus ---------------------------------

for gen=1:numgen

    Fit=testfn3(Pop);

    evolution(gen)=min(Fit);	% convergence graph of the solution

    % GA
    Best=selbest(Pop,Fit,[1,1]);
    Old=selrand(Pop,Fit,10);
    Work1=selsus(Pop,Fit,18);
    Work2=selsus(Pop,Fit,20);
    Work1=crossov(Work1,1,0);
    Work2=mutx(Work2,0.15,Space);
    Work2=muta(Work2,0.15,Delta,Space);
    Pop=[Best;Old;Work1;Work2];

end;  % gen

% -------------------------------------------------

% Best solution
figure(1);
plot(evolution,'b');
title('evolution');
xlabel('generation');
ylabel('fitness');
