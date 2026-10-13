clc;
clear;
close all;

% --- path bootstrap ---
thisFile = mfilename("fullpath");
thisDir = fileparts(thisFile);
addpath(fullfile(thisDir, "..", "genetic"));
vystupDir = fullfile(thisDir, "vystupy");
if ~exist(vystupDir, "dir")
    mkdir(vystupDir);
end
% -----------------------------------------------

% nastavenia, parametre (spolocne pre PGA aj SGA):
nvars = 10; % pocet premennych
ngen = 1000; % pocet generacii = ukoncovacia podmienka
runs = 5; % pocet parovych behov PGA/SGA
M = 512; % ohranicenie priestoru (Eggholder)
Space = [ones(1,nvars)*(-M); ones(1,nvars)*M];

fitFun = @eggholder;

% migracia: kazdy vonkajsi ostrov ma vlastny nahodny casovac 10-50 generacii
migMin = 10;
migMax = 50;

% architektura PGA - hierarchia (3 vrstvy), 9 rovnakych ostrovov x 40 = 360 jedincov
% ostrov 1 = vrchol (C), zbiera migrantov zo strednej vrstvy
% ostrovy 2-5 = stredna vrstva (O1-O4), posielaju kopiu svojho najlepsieho do C
% ostrovy 6-9 = spodna vrstva (O5-O8), posielaju kopiu svojho najlepsieho do vsetkych ostrovov strednej vrstvy
% migrant nahradi najhorsieho jedinca v cielovom ostrove

nazvy = ["C" "O1" "O2" "O3" "O4" "O5" "O6" "O7" "O8"];

% topologia: ciele{i} = ostrovy, do ktorych posiela ostrov i
ciele = {[], 1, 1, 1, 1, 2:5, 2:5, 2:5, 2:5};

selekcie = repmat("seltourn", 1, 9);

%                  nBest nOld nCross nMut  pmutx pmuta dDiv
ostrovy = repmat([   3    3    17    17    0.08  0.30   50 ], 9, 1);

% SGA - rovnake nastavenia preskalovane na 360 jedincov
sgaSel = "seltourn";
sgaPar = [27 27 153 153 0.08 0.30 50];

nIsl = size(ostrovy, 1);
assert(sum(ostrovy(:,1:4), "all") == sum(sgaPar(1:4)), ...
       "PGA a SGA musia mat rovnaky pocet jedincov");

% farby ostrovov (C modra, ziadna cervena ani cierna - tie su pre PGA best a SGA)
farby = [0.000 0.447 0.741    % C
         0.929 0.694 0.125    % O1
         0.494 0.184 0.556    % O2
         0.466 0.674 0.188    % O3
         0.301 0.745 0.933    % O4
         0.850 0.325 0.098    % O5
         0.550 0.350 0.150    % O6
         0.950 0.450 0.750    % O7
         0.500 0.500 0.500];  % O8

% nahodne seedy - nic neriadia, len sa zaznamenaju kvoli zopakovatelnosti
rng("shuffle");
seedy = randi(2^31 - 1, 1, runs);
%seedy = [2088161188 148145123 795436664 648306378 1359796845]; % pevne seedy - zakomentuj pre nahodne

% ulozisko vysledkov
evoPGA = nan(nIsl, ngen, runs); % ostrov x generacia x beh
evoSGA = nan(runs, ngen); % beh x generacia
migLog = cell(runs, 1); % [generacia zdroj ciel uspesna] pre kazdy beh
bestXP = nan(runs, nvars); bestFP = nan(runs, 1); % najlepsi jedinec PGA
bestXS = nan(runs, nvars); bestFS = nan(runs, 1); % najlepsi jedinec SGA

for r = 1:runs

    fprintf("beh %d/%d   seed = %d\n", r, runs, seedy(r));

    rng(seedy(r));

    % parovy start: SGA zacina presne zo zjednotenia populacii ostrovov
    Pops0 = cell(nIsl, 1);
    for i = 1:nIsl
        Pops0{i} = genrpop(sum(ostrovy(i,1:4)), Space);
    end

    [evoPGA(:,:,r), migLog{r}, bestXP(r,:), bestFP(r)] = ...
        spusti_pga(Pops0, ostrovy, selekcie, fitFun, Space, ngen, migMin, migMax, ciele);

    [evoSGA(r,:), bestXS(r,:), bestFS(r)] = ...
        spusti_sga(vertcat(Pops0{:}), sgaPar, sgaSel, fitFun, Space, ngen);

    fprintf("   PGA best = %10.2f   SGA best = %10.2f\n", bestFP(r), bestFS(r));

end

% ========== grafy ==========

% spolocna mierka pre vsetky grafy
yVsetko = [evoPGA(:); evoSGA(:)];
rozsah = max(yVsetko) - min(yVsetko);
yLim = [min(yVsetko) - 0.03*rozsah, max(yVsetko) + 0.03*rozsah];

popisy = [nazvy "PGA najlepsi" "SGA" "uspesna migracia do C"];

% samostatny graf pre kazdy beh (export)
for r = 1:runs

    f = figure("Visible", "off", "Position", [100 100 1000 560]);
    ax = axes(f);
    h = kresli_beh(ax, evoPGA(:,:,r), evoSGA(r,:), migLog{r}, farby, yLim, ...
                   sprintf("Beh %d (seed %d): ostrovy PGA vs SGA", r, seedy(r)));
    legend(ax, h, popisy, "Location", "eastoutside");
    exportgraphics(f, fullfile(vystupDir, sprintf("beh%d.png", r)), "Resolution", 300);
    close(f);

end

% suhrnne grafy (export)
f = figure("Visible", "off", "Position", [100 100 1000 560]);
kresli_priemer(axes(f), evoPGA, evoSGA, yLim);
exportgraphics(f, fullfile(vystupDir, "priemer.png"), "Resolution", 300);
close(f);

f = figure("Visible", "off", "Position", [100 100 1000 560]);
kresli_priemer_ostrovy(axes(f), evoPGA, farby, nazvy, yLim);
exportgraphics(f, fullfile(vystupDir, "priemer_ostrovy.png"), "Resolution", 300);
close(f);

% prehlad - vsetky behy + priemer v jednom okne
f = figure("Name", "Prehlad behov", "Position", [50 50 1700 900]);
tl = tiledlayout(f, 2, 3, "TileSpacing", "compact", "Padding", "compact");

for r = 1:runs
    ax = nexttile(tl);
    h = kresli_beh(ax, evoPGA(:,:,r), evoSGA(r,:), migLog{r}, farby, yLim, ...
                   sprintf("Beh %d (seed %d)", r, seedy(r)));
end

lg = legend(h, popisy);
lg.Layout.Tile = "east";

kresli_priemer(nexttile(tl), evoPGA, evoSGA, yLim);

exportgraphics(f, fullfile(vystupDir, "prehlad.png"), "Resolution", 200);

% ========== tabulky ==========

Tnast = tab_nastavenia(ostrovy, selekcie, sgaPar, sgaSel, nazvy, Space);
Tvys = tab_vysledky(evoPGA, evoSGA, seedy, nazvy);
Tvek = tab_vektor(bestXP, bestFP, bestXS, bestFS);
Tmig = tab_migracia(migLog, nazvy);

fprintf("\nspolocne: nvars = %d, ngen = %d, rozsah <%d, %d>, populacia = %d, migracia kazdych %d-%d gen.\n", ...
        nvars, ngen, -M, M, sum(sgaPar(1:4)), migMin, migMax);

fprintf("\n--- Tabulka 1: nastavenia ostrovov a SGA ---\n");
disp(Tnast);
fprintf("--- Tabulka 2: vysledky behov (posledna generacia) ---\n");
disp(Tvys);
fprintf("--- Tabulka 3: najlepsi najdeny vektor ---\n");
disp(Tvek);
fprintf("--- Tabulka 4: migracia (sucet za %d behov) ---\n", runs);
disp(Tmig);

writetable(Tnast, fullfile(vystupDir, "tab1_nastavenia.csv"), "Delimiter", ";");
writetable(Tvys, fullfile(vystupDir, "tab2_vysledky.csv"), "Delimiter", ";");
writetable(Tvek, fullfile(vystupDir, "tab3_vektor.csv"), "Delimiter", ";");
writetable(Tmig, fullfile(vystupDir, "tab4_migracia.csv"), "Delimiter", ";");

save(fullfile(vystupDir, "vysledky.mat"), "evoPGA", "evoSGA", "seedy", "migLog", ...
     "bestXP", "bestFP", "bestXS", "bestFS", "ostrovy", "selekcie", "sgaPar", "sgaSel", ...
     "nvars", "ngen", "runs", "Space", "migMin", "migMax");

% ========== pomocne funkcie ==========

function [evo, migLog, bestX, bestF] = spusti_pga(Pops, ostrovy, selekcie, fitFun, Space, ngen, migMin, migMax, ciele)
% ostrovny PGA s hierarchickou topologiou, migracia spodna -> stredna vrstva -> C

nIsl = numel(Pops);
Fits = cell(nIsl, 1);

% casovace migracie, stred (1) nema vlastny
dalsiaMig = inf(nIsl, 1);
for i = 2:nIsl
    dalsiaMig(i) = randi([migMin migMax]);
end

evo = nan(nIsl, ngen);
migLog = zeros(0, 4);
bestF = inf;
bestX = [];
najlepsi = inf(nIsl, 1); posledne = zeros(nIsl, 1); zrod = zeros(nIsl, 1); % pre reset

for gen = 1:ngen

    for i = 1:nIsl
        Fits{i} = fitFun(Pops{i});
        [evo(i,gen), iB] = min(Fits{i}); % konvergencna krivka ostrova

        if evo(i,gen) < bestF % globalne najlepsi jedinec PGA
            bestF = evo(i,gen);
            bestX = Pops{i}(iB,:);
        end
    end

    % migracia pred tvorbou novej generacie
    [Pops, Fits, dalsiaMig, zaznam] = migracia(Pops, Fits, dalsiaMig, gen, migMin, migMax, ciele);
    migLog = [migLog; zaznam];

    % reset stagnujuceho vonkajsieho ostrova: za 100 generacii zlepsenie < 20
    for i = 2:nIsl
        if evo(i,gen) < najlepsi(i) - 20
            najlepsi(i) = evo(i,gen);
            posledne(i) = gen;
        elseif gen - posledne(i) >= 100 && gen < 0.9*ngen
            [Pops, Fits] = posli(Pops, Fits, i, ciele{i}, gen); % najlepsi sa pred resetom posle do cielov
            Pops{i} = genrpop(size(Pops{i}, 1), Space);
            Fits{i} = fitFun(Pops{i});
            najlepsi(i) = min(Fits{i});
            posledne(i) = gen;
            zrod(i) = gen;
        end
    end

    % GA krok na kazdom ostrove, krok lokalnej mutacie sa zmensuje od vzniku ostrova
    for i = 1:nIsl
        kf = max(1 - (gen - zrod(i)) / (ngen - zrod(i)), 0.05);
        Pops{i} = krok_ga(Pops{i}, Fits{i}, ostrovy(i,:), selekcie(i), Space, kf);
    end

end

end

function [Pops, Fits, dalsiaMig, zaznam] = migracia(Pops, Fits, dalsiaMig, gen, migMin, migMax, ciele)
% ostrov posle kopiu svojho najlepsieho do vsetkych svojich cielovych ostrovov
% zaznam = [generacia zdroj ciel uspesna], uspesna = migrant lepsi ako doterajsi najlepsi v cieli

zaznam = zeros(0, 4);

for i = 2:numel(Pops)

    if gen == dalsiaMig(i)

        [Pops, Fits, z] = posli(Pops, Fits, i, ciele{i}, gen);
        zaznam = [zaznam; z];

        dalsiaMig(i) = gen + randi([migMin migMax]); % novy casovac

    end

end

end

function [Pops, Fits, zaznam] = posli(Pops, Fits, i, ciele, gen)
% kopia najlepsieho z ostrova i nahradi najhorsieho jedinca v kazdom cielovom ostrove

[migrantFit, iBest] = min(Fits{i});
migrant = Pops{i}(iBest,:);
zaznam = zeros(0, 4);

for c = ciele
    uspesna = migrantFit < min(Fits{c});
    [~, ciel] = max(Fits{c}); % najhorsi v cieli
    Pops{c}(ciel,:) = migrant;
    Fits{c}(ciel) = migrantFit;
    zaznam(end+1,:) = [gen i c uspesna];
end

end

function [evo, bestX, bestF] = spusti_sga(Pop, p, selName, fitFun, Space, ngen)
% jednoduchy jednopopulacny GA, startuje z danej populacie

evo = nan(1, ngen);
bestF = inf;
bestX = [];

for gen = 1:ngen

    Fit = fitFun(Pop);
    [evo(gen), iB] = min(Fit);

    if evo(gen) < bestF
        bestF = evo(gen);
        bestX = Pop(iB,:);
    end

    kf = max(1 - gen/ngen, 0.05); % rovnake zmensovanie kroku ako v PGA
    Pop = krok_ga(Pop, Fit, p, selName, Space, kf);

end

end

function Pop = krok_ga(Pop, Fit, p, selName, Space, kf)
% jedna generacia GA01, p = [nBest nOld nCross nMut pmutx pmuta dDiv]
% kf = faktor kroku lokalnej mutacie (1 na zaciatku, 0.05 na konci)

nBest = p(1);
nOld = p(2);
nCross = p(3);
nMut = p(4);
pmutx = p(5);
pmuta = p(6);
Delta = Space(2,:) / p(7) * kf; % aditivny krok lokalnej mutacie

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

function h = kresli_beh(ax, evoIsl, evoS, migLog, farby, yLim, titulok)
% jeden beh: vsetky ostrovy PGA, globalne najlepsi PGA, SGA a uspesne migracie

hold(ax, "on");
grid(ax, "on");

[nIsl, ngen] = size(evoIsl);
h = gobjects(nIsl + 3, 1);

for i = 2:nIsl
    h(i) = plot(ax, evoIsl(i,:), "Color", farby(i,:), "LineWidth", 0.8);
end

h(1) = plot(ax, evoIsl(1,:), "Color", farby(1,:), "LineWidth", 1.8); % stred hrubsie

h(nIsl+1) = plot(ax, min(evoIsl, [], 1), "Color", [0.85 0 0], "LineWidth", 2.5);
h(nIsl+2) = plot(ax, evoS, "k", "LineWidth", 2.5);

% uspesne migracie - efekt je viditelny v nasledujucej generacii C
usp = migLog(migLog(:,3) == 1 & migLog(:,4) == 1, 1);
gx = min(usp + 1, ngen);
h(nIsl+3) = plot(ax, gx, evoIsl(1,gx), "v", "LineStyle", "none", "MarkerSize", 6, ...
                 "MarkerEdgeColor", "k", "MarkerFaceColor", farby(1,:));

title(ax, titulok);
xlabel(ax, "generacia");
ylabel(ax, "fitness");
xlim(ax, [1 ngen]);
ylim(ax, yLim);

end

function kresli_priemer(ax, evoPGA, evoSGA, yLim)
% priemer globalne najlepsieho PGA vs priemer SGA, pas = +- smerodajna odchylka

[~, ngen, runs] = size(evoPGA);

pgaBest = reshape(min(evoPGA, [], 1), ngen, runs)'; % beh x generacia
mP = mean(pgaBest, 1); sP = std(pgaBest, 0, 1);
mS = mean(evoSGA, 1); sS = std(evoSGA, 0, 1);
x = 1:ngen;

hold(ax, "on");
grid(ax, "on");

fill(ax, [x fliplr(x)], [mP+sP fliplr(mP-sP)], [0.85 0 0], "FaceAlpha", 0.15, "EdgeColor", "none");
fill(ax, [x fliplr(x)], [mS+sS fliplr(mS-sS)], [0 0 0], "FaceAlpha", 0.10, "EdgeColor", "none");
h1 = plot(ax, mP, "Color", [0.85 0 0], "LineWidth", 2.5);
h2 = plot(ax, mS, "k", "LineWidth", 2.5);

title(ax, sprintf("Priemer %d behov: PGA vs SGA (pas = +- smer. odch.)", runs));
xlabel(ax, "generacia");
ylabel(ax, "fitness");
xlim(ax, [1 ngen]);
ylim(ax, yLim);
legend(ax, [h1 h2], ["PGA najlepsi (priemer)" "SGA (priemer)"], "Location", "northeast");

end

function kresli_priemer_ostrovy(ax, evoPGA, farby, nazvy, yLim)
% priemer 5 behov pre kazdy ostrov zvlast + priemer globalne najlepsieho PGA

[nIsl, ngen, runs] = size(evoPGA);

priemerOstrovy = mean(evoPGA, 3); % ostrov x generacia
pgaBest = reshape(min(evoPGA, [], 1), ngen, runs)';

hold(ax, "on");
grid(ax, "on");

h = gobjects(nIsl + 1, 1);
for i = 2:nIsl
    h(i) = plot(ax, priemerOstrovy(i,:), "--", "Color", farby(i,:), "LineWidth", 1);
end
h(1) = plot(ax, priemerOstrovy(1,:), "Color", farby(1,:), "LineWidth", 2);
h(nIsl+1) = plot(ax, mean(pgaBest, 1), "Color", [0.85 0 0], "LineWidth", 2.5);

title(ax, sprintf("Priemer %d behov: jednotlive ostrovy PGA", runs));
xlabel(ax, "generacia");
ylabel(ax, "fitness");
xlim(ax, [1 ngen]);
ylim(ax, yLim);
legend(ax, h, [nazvy "PGA najlepsi"], "Location", "eastoutside");

end

function T = tab_nastavenia(ostrovy, selekcie, sgaPar, sgaSel, nazvy, Space)
% parametre vsetkych ostrovov a SGA

P = [ostrovy; sgaPar];

T = table([nazvy "SGA"]', [selekcie sgaSel]', P(:,1), P(:,2), P(:,3), P(:,4), ...
          sum(P(:,1:4), 2), P(:,5), P(:,6), Space(2,1) ./ P(:,7), ...
          'VariableNames', ["ostrov" "selekcia" "nBest" "nOld" "nCross" "nMut" ...
                            "spolu" "pmutx" "pmuta" "Delta_start"]);

end

function T = tab_vysledky(evoPGA, evoSGA, seedy, nazvy)
% vysledky z poslednej generacie + priemer a smerodajna odchylka

runs = size(evoSGA, 1);

[pgaKon, iIsl] = min(evoPGA(:,end,:), [], 1);
pgaKon = reshape(pgaKon, [], 1);
iIsl = reshape(iIsl, [], 1);
sgaKon = evoSGA(:,end);
rozdiel = pgaKon - sgaKon; % zaporny = PGA lepsi

lepsi = repmat("remiza", runs, 1);
lepsi(rozdiel < 0) = "PGA";
lepsi(rozdiel > 0) = "SGA";

T = table([string(1:runs)'; "priemer"; "smer. odch."], ...
          [seedy(:); NaN; NaN], ...
          [pgaKon; mean(pgaKon); std(pgaKon)], ...
          [reshape(nazvy(iIsl), [], 1); ""; ""], ...
          [sgaKon; mean(sgaKon); std(sgaKon)], ...
          [rozdiel; mean(rozdiel); std(rozdiel)], ...
          [lepsi; ""; ""], ...
          'VariableNames', ["beh" "seed" "PGA_best" "ostrov" "SGA_best" "rozdiel" "lepsi"]);

end

function T = tab_vektor(bestXP, bestFP, bestXS, bestFS)
% najlepsi najdeny vektor PGA a SGA zo vsetkych behov

[fP, rP] = min(bestFP);
[fS, rS] = min(bestFS);

X = [bestXP(rP,:); bestXS(rS,:)];

T = array2table(X, 'VariableNames', cellstr("x" + (1:size(X, 2))));
T = addvars(T, ["PGA"; "SGA"], [rP; rS], [fP; fS], 'Before', 1, ...
            'NewVariableNames', {'algoritmus', 'beh', 'fitness'});

end
function T = tab_migracia(migLog, nazvy)
% pocet migracii a uspesnych migracii podla zdrojoveho ostrova (sucet za behy, vsetky ciele)

vsetky = vertcat(migLog{:});
nOut = numel(nazvy) - 1;
pocet = zeros(nOut, 1);
usp = zeros(nOut, 1);

for i = 2:numel(nazvy)
    pocet(i-1) = sum(vsetky(:,2) == i);
    usp(i-1) = sum(vsetky(:,2) == i & vsetky(:,4) == 1);
end

T = table([reshape(nazvy(2:end), [], 1); "spolu"], ...
          [pocet; sum(pocet)], ...
          [usp; sum(usp)], ...
          round(100 * [usp ./ pocet; sum(usp) / sum(pocet)], 1), ...
          'VariableNames', ["zdroj" "migracie" "uspesne" "uspesnost_pct"]);

end