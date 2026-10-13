% UI2 - cvicenie 5, algoritmus Diferencialna evolucia (DE)
% DE/rand/1/bin, ucelova funkcia: Schwefel 10 premennych (testfn3)
% porovnanie s GA zo zadania 1 (varianta e = GA01), rovnaka velkost populacie 50

clc;
clear;
close all;

% --- path bootstrap ---
thisFile = mfilename("fullpath");
thisDir = fileparts(thisFile);
addpath(fullfile(thisDir, "..", "genetic")); % Sekaj toolbox
addpath(fullfile(thisDir, ".."));            % testfn3.m
outDir = fullfile(thisDir, "vystupy");
if ~exist(outDir, "dir"), mkdir(outDir); end
% -----------------------------------------------

% nastavenia, parametre:
nvars = 10;     % pocet premennych
nPop = 50;      % velkost populacie
maxGen = 500;   % pocet generacii (50 ohodnoteni fitness na generaciu, rovnako ako GA)
runs = 5;       % pocet behov
M = 500;        % ohranicenie priestoru
lb = -M*ones(1,nvars); ub = M*ones(1,nvars);
Space = [lb; ub];
fitFun = @testfn3;
seedy = [2088161188 148145123 795436664 648306378 1359796845]; % rovnake ako v zadani 4
fOpt = -418.9829*nvars;  % globalne minimum Schwefel

% GA zo zadania 1, varianta e (GA01): [nBest nOld nCross nMut pmutx pmuta dDiv]
gaPar = [2 10 18 20 0.15 0.15 100];

% ========== 1. experimentalne hladanie parametrov F, CR ==========
Fset = [0.3 0.5 0.7 0.9];
CRset = [0.1 0.3 0.5 0.7 0.9];
gridMean = zeros(numel(Fset), numel(CRset));

fprintf("Hladanie parametrov DE (priemer %d behov, %d generacii)\n", runs, maxGen);
fprintf("%8s", "F \\ CR"); fprintf("%10.1f", CRset); fprintf("\n");
for a = 1:numel(Fset)
    fprintf("%8.1f", Fset(a));
    for b = 1:numel(CRset)
        res = zeros(1, runs);
        for r = 1:runs
            rng(seedy(r));
            Pop0 = genrpop(nPop, Space);
            e = spusti_de(Pop0, fitFun, lb, ub, maxGen, Fset(a), CRset(b));
            res(r) = e(end);
        end
        gridMean(a,b) = mean(res);
        fprintf("%10.1f", gridMean(a,b));
    end
    fprintf("\n");
end

[~, idx] = min(gridMean(:));
[ia, ib] = ind2sub(size(gridMean), idx);
F = Fset(ia); CR = CRset(ib);
fprintf("\nNajlepsie nastavenie: F = %.1f, CR = %.1f (priemer %.2f)\n\n", F, CR, gridMean(ia,ib));

writetable(array2table(gridMean, "VariableNames", "CR_" + string(CRset), "RowNames", "F_" + string(Fset)), ...
    fullfile(outDir, "tab1_mriezka_F_CR.csv"), "WriteRowNames", true);

% ========== 2. vysledny experiment: DE vs GA, 5 behov ==========
evoDE = nan(runs, maxGen);
evoGA = nan(runs, maxGen);

for r = 1:runs
    rng(seedy(r));
    Pop0 = genrpop(nPop, Space); % rovnaka startovacia populacia pre DE aj GA
    evoDE(r,:) = spusti_de(Pop0, fitFun, lb, ub, maxGen, F, CR);
    evoGA(r,:) = spusti_ga(Pop0, gaPar, fitFun, Space, maxGen);
    fprintf("beh %d/%d: DE = %9.2f   GA = %9.2f\n", r, runs, evoDE(r,end), evoGA(r,end));
end

fprintf("\nglobalne minimum = %.2f\n", fOpt);
fprintf("%-6s %10s %10s %10s\n", "", "najlepsi", "priemer", "najhorsi");
fprintf("%-6s %10.2f %10.2f %10.2f\n", "DE", min(evoDE(:,end)), mean(evoDE(:,end)), max(evoDE(:,end)));
fprintf("%-6s %10.2f %10.2f %10.2f\n", "GA", min(evoGA(:,end)), mean(evoGA(:,end)), max(evoGA(:,end)));

T = table((1:runs)', evoDE(:,end), evoGA(:,end), 'VariableNames', ["beh" "DE" "GA"]);
writetable(T, fullfile(outDir, "tab2_vysledky.csv"));

% spolocna mierka y
yAll = [evoDE(:); evoGA(:)];
yLim = [min(yAll) max(yAll)] + [-1 1]*0.03*(max(yAll) - min(yAll));

% ========== GRAF 0: mriezka F x CR ==========
f0 = figure("Name", "Mriezka F x CR", "Position", [100 100 800 500]);
imagesc(CRset, Fset, gridMean); colorbar; colormap(parula);
set(gca, "YDir", "normal", "XTick", CRset, "YTick", Fset);
for a = 1:numel(Fset)
    for b = 1:numel(CRset)
        text(CRset(b), Fset(a), sprintf("%.0f", gridMean(a,b)), "HorizontalAlignment", "center", "Color", "w", "FontWeight", "bold");
    end
end
xlabel("CR"); ylabel("F");
title(sprintf("DE: priemerna konecna fitness (%d behov) pre rozne F a CR", runs));
exportgraphics(f0, fullfile(outDir, "0_mriezka_F_CR.png"), "Resolution", 300);

% ========== GRAF 1: DE, 5 behov + priemer ==========
f1 = kresli_behy(evoDE, sprintf("DE/rand/1/bin (F = %.1f, CR = %.1f) - Schwefel 10 premennych, populacia 50", F, CR), "generacia", yLim);
exportgraphics(f1, fullfile(outDir, "1_DE.png"), "Resolution", 300);

% ========== GRAF 2: GA, 5 behov + priemer ==========
f2 = kresli_behy(evoGA, "GA (zadanie 1, GA01) - Schwefel 10 premennych, populacia 50", "generacia", yLim);
exportgraphics(f2, fullfile(outDir, "2_GA.png"), "Resolution", 300);

% ========== GRAF 3: DE vs GA ==========
f3 = figure("Name", "DE vs GA", "Position", [150 150 900 550]);
hold on; grid on;
hD = plot(evoDE', "Color", [0.55 0.70 1.00], "LineWidth", 0.8);
hG = plot(evoGA', "Color", [1.00 0.65 0.60], "LineWidth", 0.8);
hDm = plot(mean(evoDE, 1), "Color", [0 0.25 0.8], "LineWidth", 2.5);
hGm = plot(mean(evoGA, 1), "Color", [0.8 0.1 0], "LineWidth", 2.5);
yline(fOpt, "k--");
title("Porovnanie: DE vs GA (populacia 50, 5 behov)");
xlabel("generacia"); ylabel("fitness"); xlim([1 maxGen]); ylim(yLim);
legend([hD(1) hG(1) hDm hGm], ["DE behy" "GA behy" "DE priemer" "GA priemer"]);
exportgraphics(f3, fullfile(outDir, "3_porovnanie.png"), "Resolution", 300);

save(fullfile(outDir, "vysledky.mat"), "evoDE", "evoGA", "gridMean", "Fset", "CRset", "F", "CR");

% ========== funkcie ==========

function f = kresli_behy(evo, nazov, xlab, yLim)
% 5 behov + priemer do jedneho grafu
f = figure("Name", nazov, "Position", [100 100 900 550]);
hold on; grid on;
plot(evo', "LineWidth", 1);                    % jednotlive behy
plot(mean(evo, 1), "k", "LineWidth", 2.5);     % priemer
title(nazov);
xlabel(xlab); ylabel("fitness"); xlim([1 size(evo,2)]); ylim(yLim);
legend(["beh " + (1:size(evo,1)), "priemer"]);
end

function evolution = spusti_de(X, fitFun, lb, ub, maxGen, F, CR)
% DE/rand/1/bin, X = pociatocna populacia

[nPop, nvars] = size(X);
Fit = fitFun(X);
evolution = nan(1, maxGen);

for gen = 1:maxGen
    U = X;
    for i = 1:nPop
        % mutacia: 3 nahodni jedinci, rozni navzajom aj od i
        r = randperm(nPop - 1, 3);
        r(r >= i) = r(r >= i) + 1;
        V = X(r(1),:) + F*(X(r(2),:) - X(r(3),:));

        % binomicke krizenie, aspon 1 gen z mutanta
        mask = rand(1, nvars) < CR;
        mask(randi(nvars)) = true;
        U(i,mask) = V(mask);
    end
    U = max(min(U, ub), lb);       % orezanie na hranice priestoru

    % selekcia 1:1 - potomok nahradi rodica, ak je lepsi alebo rovnaky
    FitU = fitFun(U);
    lepsi = FitU <= Fit;
    X(lepsi,:) = U(lepsi,:);
    Fit(lepsi) = FitU(lepsi);

    evolution(gen) = min(Fit);
end

end

function evo = spusti_ga(Pop, p, fitFun, Space, ngen)
% GA zo zadania 1 (varianta e = GA01), p = [nBest nOld nCross nMut pmutx pmuta dDiv]

Delta = Space(2,:) / p(7); % aditivny krok lokalnej mutacie
evo = nan(1, ngen);

for gen = 1:ngen
    Fit = fitFun(Pop);
    evo(gen) = min(Fit); % konvergencna krivka

    Best = selbest(Pop, Fit, ones(1,p(1)));
    Old = selrand(Pop, Fit, p(2));
    Work1 = selsus(Pop, Fit, p(3));
    Work2 = selsus(Pop, Fit, p(4));
    Work1 = crossov(Work1, 1, 0);
    Work2 = mutx(Work2, p(5), Space);
    Work2 = muta(Work2, p(6), Delta, Space);
    Pop = [Best; Old; Work1; Work2];
end

end
