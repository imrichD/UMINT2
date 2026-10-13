clc;
clear;
close all;

% --- path bootstrap ---
thisFile = mfilename("fullpath");
thisDir = fileparts(thisFile);
addpath(fullfile(thisDir, "..", "genetic"));
% -----------------------------------------------

% nastavenia, parametre:
nvars = 10; % pocet premennych
ngen = 500; % pocet generacii
runs = 5; % pocet behov experimentu PGA/SGA
M = 512; % ohranicenie priestoru (Eggholder)
Space = [ones(1,nvars)*(-M); ones(1,nvars)*M];

fitFun = @eggholder;

% migracia: kazdy vonkajsi ostrov ma vlastny nahodny casovac 10-50 generacii
migMin = 10;
migMax = 50;

% architektura PGA - hviezda (5 ostrovov x 30 = 150 jedincov)
% ostrov 1 = stred (C), zbiera migrantov, vyvazeny so sklonom k exploatacii
% ostrovy 2-5 = vonkajsie (O1-O4), explorery, posielaju svojho najlepsieho do C
% migrant nahradi nahodneho jedinca v C, nikdy nie nBest najlepsich z C

nazvy = ["C" "O1" "O2" "O3" "O4"];

selekcie = ["seltourn" "selrand" "selsus" "seltourn" "selsus"];

%           -- pocty jedincov (spolu 30) --   ---- mutacia ----
%           nBest  nOld  nCross  nMut         pmutx pmuta dDiv
ostrovy = [   2     3     13     12           0.08  0.15  200     % C  stred, vyvazeny
              1     6     11     12           0.25  0.40   10     % O1 silny explorer
              1     4     12     13           0.20  0.30   20     % O2 explorer
              1     5     10     14           0.30  0.20   50     % O3 explorer, silna mutacia
              1     8      8     13           0.15  0.40   25 ];  % O4 explorer, vela nahodnych

% SGA - nastavenia stredneho ostrova, preskalovane na 150 jedincov
sgaSel = "seltourn";
sgaPar = [10 15 65 60 0.08 0.15 200];

nIsl = size(ostrovy, 1);
assert(sum(ostrovy(:,1:4), "all") == sum(sgaPar(1:4)), ...
       "PGA a SGA musia mat rovnaky pocet jedincov");

farby = lines(nIsl);

% ulozisko vysledkov
evoPGA = nan(nIsl, ngen, runs); % ostrov x generacia x beh
evoSGA = nan(runs, ngen); % beh x generacia

for r = 1:runs

    fprintf("beh %d/%d\n", r, runs);

    evoPGA(:,:,r) = spusti_pga(ostrovy, selekcie, fitFun, Space, ngen, migMin, migMax);
    evoSGA(r,:) = spusti_sga(sgaPar, sgaSel, fitFun, Space, ngen);

    fprintf("   PGA best = %10.2f   SGA best = %10.2f\n", ...
            min(evoPGA(:,end,r)), evoSGA(r,end));

    kresli_beh(evoPGA(:,:,r), evoSGA(r,:), nazvy, farby, r, ...
               fullfile(thisDir, sprintf("beh%d.png", r)));

end

% posledny obrazok - priemery 5 behov
kresli_priemer(evoPGA, evoSGA, nazvy, farby, fullfile(thisDir, "priemer.png"));

% tabulka k diskusii
vypis(evoPGA, evoSGA, nazvy);

% ========== pomocne funkcie ==========

function evo = spusti_pga(ostrovy, selekcie, fitFun, Space, ngen, migMin, migMax)
% ostrovny PGA s hviezdicovou topologiou, migracia O -> C

nIsl = size(ostrovy, 1);
Pops = cell(nIsl, 1);
Fits = cell(nIsl, 1);

for i = 1:nIsl
    Pops{i} = genrpop(sum(ostrovy(i,1:4)), Space); % inicializacia ostrova
end

% casovace migracie, stred (1) nema vlastny
dalsiaMig = inf(nIsl, 1);
for i = 2:nIsl
    dalsiaMig(i) = randi([migMin migMax]);
end

evo = nan(nIsl, ngen);

for gen = 1:ngen

    for i = 1:nIsl
        Fits{i} = fitFun(Pops{i});
        evo(i,gen) = min(Fits{i}); % konvergencna krivka ostrova
    end

    % migracia pred tvorbou novej generacie
    [Pops, Fits, dalsiaMig] = migracia(Pops, Fits, dalsiaMig, gen, ...
                                       ostrovy(1,1), migMin, migMax);

    % GA krok na kazdom ostrove
    for i = 1:nIsl
        Pops{i} = krok_ga(Pops{i}, Fits{i}, ostrovy(i,:), selekcie(i), Space);
    end

end

end

function [Pops, Fits, dalsiaMig] = migracia(Pops, Fits, dalsiaMig, gen, nBestC, migMin, migMax)
% vonkajsi ostrov posle kopiu svojho najlepsieho do stredu,
% ten nahradi nahodneho jedinca v strede okrem nBestC najlepsich

for i = 2:numel(Pops)

    if gen == dalsiaMig(i)

        [migrantFit, iBest] = min(Fits{i});
        migrant = Pops{i}(iBest,:);

        [~, poradie] = sort(Fits{1}); % vzostupne, najlepsi na zaciatku
        volne = poradie(nBestC+1:end);
        ciel = volne(randi(numel(volne)));

        Pops{1}(ciel,:) = migrant;
        Fits{1}(ciel) = migrantFit;

        dalsiaMig(i) = gen + randi([migMin migMax]); % novy casovac

    end

end

end

function evo = spusti_sga(p, selName, fitFun, Space, ngen)
% jednoduchy jednopopulacny GA

Pop = genrpop(sum(p(1:4)), Space);
evo = nan(1, ngen);

for gen = 1:ngen

    Fit = fitFun(Pop);
    evo(gen) = min(Fit);
    Pop = krok_ga(Pop, Fit, p, selName, Space);

end

end

function Pop = krok_ga(Pop, Fit, p, selName, Space)
% jedna generacia GA01, p = [nBest nOld nCross nMut pmutx pmuta dDiv]

nBest = p(1);
nOld = p(2);
nCross = p(3);
nMut = p(4);
pmutx = p(5);
pmuta = p(6);
Delta = Space(2,:) / p(7); % aditivny krok lokalnej mutacie

Best = selbest(Pop, Fit, ones(1,nBest));

if nOld > 0
    Old = selrand(Pop, Fit, nOld);
else
    Old = [];
end

Work1 = vyber(selName, Pop, Fit, nCross);
Work2 = vyber(selName, Pop, Fit, nMut);

Work1 = crossov(Work1, 1, 0);
Work2 = mutx(Work2, pmutx, Space);
Work2 = muta(Work2, pmuta, Delta, Space);

Pop = [Best; Old; Work1; Work2];

end

function Sel = vyber(name, Pop, Fit, n)
% vyberie selekcnu funkciu podla nazvu

switch name

    case "selsus"
        Sel = selsus(Pop, Fit, n);

    case "seltourn"
        Sel = seltourn(Pop, Fit, n);

    case "selrand"
        Sel = selrand(Pop, Fit, n);

    otherwise
        error("neznama selekcia: %s", name);

end

end

function kresli_beh(evoIsl, evoS, nazvy, farby, r, subor)
% jeden beh: vsetky ostrovy PGA (roznymi farbami) + SGA

f = figure("Name", sprintf("Beh %d", r), "Position", [100 100 900 550]);
hold on;
grid on;

for i = 1:size(evoIsl, 1)

    if i == 1
        sirka = 2; % stredny ostrov hrubsie
    else
        sirka = 1.2;
    end

    plot(evoIsl(i,:), "Color", farby(i,:), "LineWidth", sirka);

end

plot(evoS, "k", "LineWidth", 2);

title(sprintf("Beh %d: ostrovy PGA vs SGA (Eggholder, 10 premennych)", r));
xlabel("generacia");
ylabel("fitness");
axis tight;
legend([nazvy "SGA"], "Location", "northeast");

exportgraphics(f, subor, "Resolution", 300);

end

function kresli_priemer(evoPGA, evoSGA, nazvy, farby, subor)
% priemer 5 behov: kazdy ostrov, najlepsi z PGA a SGA

[nIsl, ngen, runs] = size(evoPGA);

priemerOstrovy = mean(evoPGA, 3); % ostrov x generacia
pgaBest = reshape(min(evoPGA, [], 1), ngen, runs)'; % beh x generacia
priemerPGA = mean(pgaBest, 1);
priemerSGA = mean(evoSGA, 1);

f = figure("Name", "Priemer behov", "Position", [100 100 900 550]);
hold on;
grid on;

for i = 1:nIsl
    plot(priemerOstrovy(i,:), "--", "Color", farby(i,:), "LineWidth", 1);
end

plot(priemerPGA, "r", "LineWidth", 2.5);
plot(priemerSGA, "k", "LineWidth", 2.5);

title(sprintf("Priemer %d behov: PGA vs SGA (Eggholder, 10 premennych)", runs));
xlabel("generacia");
ylabel("fitness");
axis tight;
legend([nazvy "PGA (najlepsi ostrov)" "SGA"], "Location", "northeast");

exportgraphics(f, subor, "Resolution", 300);

end

function vypis(evoPGA, evoSGA, nazvy)
% tabulka vysledkov z poslednej generacie

runs = size(evoSGA, 1);
pgaKon = zeros(runs, 1);

fprintf("\n%-6s %12s %8s %12s\n", "beh", "PGA best", "ostrov", "SGA best");

for r = 1:runs
    [pgaKon(r), iIsl] = min(evoPGA(:,end,r));
    fprintf("%-6d %12.2f %8s %12.2f\n", r, pgaKon(r), nazvy(iIsl), evoSGA(r,end));
end

sgaKon = evoSGA(:,end);

fprintf("\n%-6s %12s %12s %12s\n", "", "najlepsi", "priemer", "najhorsi");
fprintf("%-6s %12.2f %12.2f %12.2f\n", "PGA", min(pgaKon), mean(pgaKon), max(pgaKon));
fprintf("%-6s %12.2f %12.2f %12.2f\n", "SGA", min(sgaKon), mean(sgaKon), max(sgaKon));

end