% UI2 - cvicenie 4, algoritmus PSO
% PSO prevzate z Clanok.m (Mech), ucelova funkcia: Schwefel 10 premennych (testfn3)
% porovnanie s GA zo zadania 1 (varianta e = GA01), rovnaka velkost populacie 50

clc;
clear;
close all;

% ===== PREPINAC =====
nScout = 0;   % 0 = zakladne PSO, 5 = PSO s prieskumnikmi (najhorsie castice skocia do gBest s 1 nahodnym genom)
% ====================

% --- path bootstrap ---
thisFile = mfilename("fullpath");
thisDir = fileparts(thisFile);
addpath(fullfile(thisDir, "..", "genetic")); % Sekaj toolbox
addpath(fullfile(thisDir, ".."));            % testfn3.m
if nScout == 0
    outDir = fullfile(thisDir, "vystupy_bez_scout");
    verzia = "zakladne PSO";
else
    outDir = fullfile(thisDir, "vystupy_scout");
    verzia = sprintf("PSO so %d prieskumnikmi", nScout);
end
if ~exist(outDir, "dir"), mkdir(outDir); end
% -----------------------------------------------

% nastavenia, parametre:
nvars = 10;     % pocet premennych
nPart = 50;     % pocet castic (= velkost populacie GA)
maxIter = 500;  % pocet iteracii / generacii (ako v zadani 1)
runs = 5;       % pocet behov
M = 500;        % ohranicenie priestoru
lb = -M*ones(1,nvars); ub = M*ones(1,nvars);
Space = [lb; ub];
fitFun = @testfn3;
seedy = [2088161188 148145123 795436664 648306378 1359796845]; % pre opakovatelnost

% PSO (ako v Clanok.m, len w linearne klesa)
wMax = 0.9;               % zotrvacnost na zaciatku - prehladavanie
wMin = 0.4;               % zotrvacnost na konci - doladovanie
c1 = 1.5;                 % kognitivna zlozka (vlastne najlepsie)
c2 = 1.5;                 % socialna zlozka (najlepsie roja)
vmax = 0.2*(ub - lb);     % obmedzenie rychlosti

% GA zo zadania 1, varianta e (GA01): [nBest nOld nCross nMut pmutx pmuta dDiv]
gaPar = [2 10 18 20 0.15 0.15 100];

% ========== behy ==========
evoPSO = nan(runs, maxIter);
evoGA = nan(runs, maxIter);

for r = 1:runs

    rng(seedy(r));
    Pop0 = genrpop(nPart, Space); % rovnaka startovacia populacia pre PSO aj GA

    if r == 1
        % z 1. behu si ulozime aj historiu poloh castic (na vizualizaciu roja)
        [evoPSO(r,:), Xhist, gHist] = spusti_pso(Pop0, fitFun, lb, ub, maxIter, wMax, wMin, c1, c2, vmax, nScout);
    else
        evoPSO(r,:) = spusti_pso(Pop0, fitFun, lb, ub, maxIter, wMax, wMin, c1, c2, vmax, nScout);
    end
    evoGA(r,:) = spusti_ga(Pop0, gaPar, fitFun, Space, maxIter);

    fprintf("beh %d/%d: PSO = %9.2f   GA = %9.2f\n", r, runs, evoPSO(r,end), evoGA(r,end));

end

% tabulka vysledkov (globalne minimum Schwefel = -418.98*10 = -4189.8)
fprintf("\n%s\n", verzia);
fprintf("%-6s %10s %10s %10s\n", "", "najlepsi", "priemer", "najhorsi");
fprintf("%-6s %10.2f %10.2f %10.2f\n", "PSO", min(evoPSO(:,end)), mean(evoPSO(:,end)), max(evoPSO(:,end)));
fprintf("%-6s %10.2f %10.2f %10.2f\n", "GA",  min(evoGA(:,end)),  mean(evoGA(:,end)),  max(evoGA(:,end)));

% spolocna mierka y pre vsetky 3 grafy
yAll = [evoPSO(:); evoGA(:)];
yLim = [min(yAll) max(yAll)] + [-1 1]*0.03*(max(yAll) - min(yAll));

% ========== GRAF 1: PSO, 5 behov + priemer ==========
f1 = kresli_behy(evoPSO, "PSO - " + verzia + ", Schwefel 10 premennych, 50 castic", "iteracia", yLim);
exportgraphics(f1, fullfile(outDir, "1_PSO.png"), "Resolution", 300);

% ========== GRAF 2: GA, 5 behov + priemer ==========
f2 = kresli_behy(evoGA, "GA (zadanie 1, GA01) - Schwefel 10 premennych, populacia 50", "generacia", yLim);
exportgraphics(f2, fullfile(outDir, "2_GA.png"), "Resolution", 300);

% ========== GRAF 3: PSO vs GA ==========
f3 = figure("Name", "PSO vs GA", "Position", [150 150 900 550]);
hold on; grid on;
hP = plot(evoPSO', "Color", [0.55 0.70 1.00], "LineWidth", 0.8);
hG = plot(evoGA', "Color", [1.00 0.65 0.60], "LineWidth", 0.8);
hPm = plot(mean(evoPSO, 1), "Color", [0 0.25 0.8], "LineWidth", 2.5);
hGm = plot(mean(evoGA, 1), "Color", [0.8 0.1 0], "LineWidth", 2.5);
title("Porovnanie: " + verzia + " vs GA (populacia 50, 5 behov)");
xlabel("iteracia / generacia"); ylabel("fitness"); xlim([1 maxIter]); ylim(yLim);
legend([hP(1) hG(1) hPm hGm], ["PSO behy" "GA behy" "PSO priemer" "GA priemer"]);
exportgraphics(f3, fullfile(outDir, "3_porovnanie.png"), "Resolution", 300);

% ========== GRAF 4: pohyb roja (beh 1) ==========
% priestor ma 10 rozmerov -> zobrazime premenne x1, x2 nad plochou Schwefel
% Schwefel je separabilna, preto vyska bodu = presny prispevok x1, x2 k fitness
[A, B] = meshgrid(linspace(-M, M, 200));
Z = -A.*sin(sqrt(abs(A))) - B.*sin(sqrt(abs(B)));
schw2 = @(x1, x2) -x1.*sin(sqrt(abs(x1))) - x2.*sin(sqrt(abs(x2)));
xOpt = 420.9687;

snaps = [0 20 100 maxIter];
f4 = figure("Name", "Pohyb roja", "Position", [50 50 1200 900], "Color", "w");
for s = 1:numel(snaps)
    subplot(2, 2, s);
    surf(A, B, Z, "EdgeColor", "none"); hold on;
    shading interp; colormap(parula);
    P = Xhist(:,:,snaps(s)+1);
    g = gHist(snaps(s)+1,:);
    % castice su nad plochou (zTop), ciara ukazuje ich miesto na ploche
    zTop = max(Z(:)) + 300;
    zP = schw2(P(:,1), P(:,2));
    plot3([P(:,1) P(:,1)]', [P(:,2) P(:,2)]', [zP zTop*ones(size(zP))]', "-", "Color", [0.3 0.3 0.3]);
    scatter3(P(:,1), P(:,2), zTop*ones(size(zP)), 30, "k", "filled");
    plot3(g(1), g(2), zTop + 40, "rp", "MarkerSize", 18, "MarkerFaceColor", "r");
    plot3([xOpt xOpt], [xOpt xOpt], [schw2(xOpt, xOpt) zTop], "m--", "LineWidth", 1.5);
    plot3(xOpt, xOpt, zTop, "mo", "MarkerSize", 12, "LineWidth", 2.5);
    axis([-M M -M M min(Z(:)) zTop+80]); view(-35, 40);
    xlabel("x_1"); ylabel("x_2"); zlabel("f(x_1) + f(x_2)");
    title(sprintf("iteracia %d   (gBest F = %.1f)", snaps(s), evoPSOpos(evoPSO(1,:), snaps(s), Xhist, fitFun)));
end
sgtitle(verzia + ": pohyb roja v rovine x_1, x_2 (beh 1)   cierna = castice, cervena = gBest, fialovy kruh = globalne minimum");
print(f4, "-dpng", "-r200", fullfile(outDir, "4_pohyb_roja.png")); % print namiesto exportgraphics (spravne vykresli plochu)

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

function F = evoPSOpos(evo, it, Xhist, fitFun)
% fitness gBest v iteracii it (it = 0 -> najlepsia castica startovacej populacie)
if it == 0
    F = min(fitFun(Xhist(:,:,1)));
else
    F = evo(it);
end
end

function [evolution, Xhist, gHist] = spusti_pso(X, fitFun, lb, ub, maxIter, wMax, wMin, c1, c2, vmax, nScout)
% PSO z Clanok.m, X = pociatocne polohy castic

[nPart, nvars] = size(X);
V = zeros(nPart, nvars);                % pociatocne rychlosti
Fit = fitFun(X);
pBest = X;  pBestF = Fit;               % najlepsie miesto kazdej castice
[gBestF, ib] = min(Fit);
gBest = X(ib,:);                        % najlepsie miesto roja
evolution = nan(1, maxIter);

Xhist = zeros(nPart, nvars, maxIter+1); % historia poloh castic
gHist = zeros(maxIter+1, nvars);        % historia gBest
Xhist(:,:,1) = X; gHist(1,:) = gBest;

for it = 1:maxIter
    w = wMax - (wMax - wMin)*(it - 1)/(maxIter - 1);     % klesajuca zotrvacnost
    r1 = rand(nPart,nvars); r2 = rand(nPart,nvars);
    V = w*V + c1*r1.*(pBest - X) + c2*r2.*(gBest - X);   % aktualizacia rychlosti
    V = max(min(V, vmax), -vmax);                        % obmedzenie rychlosti
    X = max(min(X + V, ub), lb);                         % nova poloha v hraniciach

    if nScout > 0
        % agresivny prieskum (Schwefel je separabilna): nScout najhorsich castic
        % skoci do gBest a jeden nahodny gen dostane nahodnu hodnotu z celeho priestoru
        [~, poradie] = sort(pBestF, "descend");
        sc = poradie(1:nScout);
        X(sc,:) = repmat(gBest, nScout, 1);
        for s = sc'
            j = randi(nvars);
            X(s,j) = lb(j) + rand*(ub(j) - lb(j));
        end
        V(sc,:) = 0;
    end

    Fit = fitFun(X);
    lepsi = Fit < pBestF;                                % lepsie vlastne minimum
    pBest(lepsi,:) = X(lepsi,:); pBestF(lepsi) = Fit(lepsi);

    [minF, ib] = min(pBestF);
    if minF < gBestF                                     % lepsie minimum roja
        gBestF = minF; gBest = pBest(ib,:);
    end
    evolution(it) = gBestF;
    Xhist(:,:,it+1) = X; gHist(it+1,:) = gBest;
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
