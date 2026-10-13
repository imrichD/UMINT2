clc;
clear;
close all;

% --- cesty ---
thisDir = fileparts(mfilename("fullpath"));
addpath(fullfile(thisDir, "..", "genetic"));
vystupDir = fullfile(thisDir, "vystupy_cpga"); % novy priecinok, stare vystupy zostanu nedotknute
if ~exist(vystupDir, "dir"), mkdir(vystupDir); end

% ========== nastavenia (spolocne pre PGA, CPGA aj SGA) ==========
nvars = 10;   % pocet premennych
ngen = 1000;  % pocet generacii = ukoncovacia podmienka
runs = 5;     % pocet behov PGA/CPGA/SGA
M = 512;      % ohranicenie priestoru (Eggholder)
Space = [-M*ones(1,nvars); M*ones(1,nvars)];
fitFun = @eggholder;

migMin = 10; migMax = 50;   % nahodny casovac migracie kazdeho vonkajsieho ostrova
stagK = 100; stagEps = 20;  % reset ostrova: za stagK generacii zlepsenie < stagEps

% architektura PGA - hviezda, 9 rovnakych ostrovov x 40 = 360 jedincov
% ostrov 1 = stred (C), zbiera migrantov, nikdy sa neresetuje
% ostrovy 2-9 = vonkajsie (O1-O8), posielaju kopiu svojho najlepsieho do C
% migrant nahradi najhorsieho jedinca v C
nazvy = ["C" "O1" "O2" "O3" "O4" "O5" "O6" "O7" "O8"];
nIsl = numel(nazvy);
sel = @seltourn; % selekcia pre vsetky ostrovy aj SGA

%         nBest nOld nCross nMut  pmutx pmuta dDiv
ostrov = [  3    3    17    17    0.08  0.30   50 ]; % kazdy ostrov (40 jedincov)
sgaPar = [ 27   27   153   153    0.08  0.30   50 ]; % SGA (360 jedincov)

% architektura CPGA - 2D toroidna mriezka, okolie Moore (8 susedov),
% partner = vitaz turnaja 2 nahodnych susedov, krizenie 1-bodove,
% potomok nahradi bunku iba ak je lepsi
% asynchronna aktualizacia: bunky sa v kazdej generacii prechadzaju v nahodnom poradi,
% zmena bunky je hned viditelna pre dalsie bunky
% mutuje sa kazdy potomok, preto nizsie pravdepodobnosti (efektivne ako v PGA/SGA)
%          riadky stlpce pmutx pmuta dDiv
cpgaPar = [ 18    20    0.04  0.15   50 ];
snapGen = [1 25 100 300 ngen]; % generacie pre snimky mriezky
cpgaOkolie = "moore";          % hlavny CPGA podla prednasky (8 susedov)

% varianty CPGA na porovnanie: tvar mriezky x okolie (vsetky maju 360 jedincov)
% varianta 1 = hlavny CPGA (prevezme sa, nespusta sa znova)

varMriezky = [cpgaPar(1:2); 9 40; 6 60];
varOkolia = ["moore" "neumann"];
[im, io] = ndgrid(1:size(varMriezky, 1), 1:numel(varOkolia));
varPar = [varMriezky(im(:),:) repmat(cpgaPar(3:5), numel(im), 1)];
varOk = reshape(varOkolia(io(:)), [], 1);
nVar = size(varPar, 1);
assert(all(prod(varPar(:,1:2), 2) == sum(sgaPar(1:4))), "varianty CPGA musia mat rovnaky pocet jedincov");
assert(varOk(1) == cpgaOkolie, "varianta 1 musi byt hlavny CPGA");

assert(nIsl*sum(ostrov(1:4)) == sum(sgaPar(1:4)), "PGA a SGA musia mat rovnaky pocet jedincov");
assert(prod(cpgaPar(1:2)) == sum(sgaPar(1:4)), "CPGA musi mat rovnaky pocet jedincov");

% farby ostrovov: C modra, O1-O8 (cervena = PGA best, cierna = SGA, tmavozelena = CPGA)
farby = [0 .447 .741; .929 .694 .125; .494 .184 .556; .466 .674 .188; .301 .745 .933;
         .85 .325 .098; .55 .35 .15; .95 .45 .75; .5 .5 .5];

rng("shuffle"); seedy = randi(2^31 - 1, 1, runs); % nahodne seedy
seedy = [2088161188 148145123 795436664 648306378 1359796845]; % pevne seedy - zakomentuj pre nahodne

% ========== experiment ==========
evoPGA = nan(nIsl, ngen, runs); % ostrov x generacia x beh
evoSGA = nan(runs, ngen);       % beh x generacia
evoCPGA = nan(runs, ngen);      % beh x generacia
migLog = cell(runs, 1);         % [generacia zdroj uspesna] pre kazdy beh
snapC = cell(runs, 1);          % snimky fitness mriezky CPGA
evoVar = nan(nVar, ngen, runs); % varianta CPGA x generacia x beh
bestXP = nan(runs, nvars); bestFP = nan(runs, 1); % najlepsi jedinec PGA
bestXS = nan(runs, nvars); bestFS = nan(runs, 1); % najlepsi jedinec SGA
bestXC = nan(runs, nvars); bestFC = nan(runs, 1); % najlepsi jedinec CPGA

for r = 1:runs

    fprintf("beh %d/%d   seed = %d\n", r, runs, seedy(r));
    rng(seedy(r));

    % spolocny start: SGA aj CPGA zacinaju presne zo zjednotenia populacii ostrovov
    Pops0 = cell(nIsl, 1);
    for i = 1:nIsl
        Pops0{i} = genrpop(sum(ostrov(1:4)), Space);
    end
    Pop0 = vertcat(Pops0{:});

    [evoPGA(:,:,r), migLog{r}, bestXP(r,:), bestFP(r)] = ...
        spusti_pga(Pops0, ostrov, sel, fitFun, Space, ngen, migMin, migMax, stagK, stagEps);
    [evoSGA(r,:), bestXS(r,:), bestFS(r)] = ...
        spusti_sga(Pop0, sgaPar, sel, fitFun, Space, ngen);
    [evoCPGA(r,:), bestXC(r,:), bestFC(r), snapC{r}] = ...
        spusti_cpga(Pop0, cpgaPar, cpgaOkolie, fitFun, Space, ngen, snapGen);

    fprintf("   PGA best = %10.2f   CPGA best = %10.2f   SGA best = %10.2f\n", ...
            bestFP(r), bestFC(r), bestFS(r));

    % varianty CPGA - spustaju sa az po hlavnych algoritmoch, aby sa ich vysledky nezmenili
    evoVar(1,:,r) = evoCPGA(r,:);
    for v = 2:nVar
        evoVar(v,:,r) = spusti_cpga(Pop0, varPar(v,:), varOk(v), fitFun, Space, ngen, []);
        fprintf("   CPGA %dx%d %-7s = %10.2f\n", varPar(v,1), varPar(v,2), varOk(v), evoVar(v,end,r));
    end

end

% kontrola: PGA a SGA maju byt zhodne s vysledkami z predchadzajuceho zadania
starySubor = fullfile(thisDir, "vystupy", "vysledky.mat");
if isfile(starySubor)
    S = load(starySubor, "evoPGA", "evoSGA", "seedy");
    if isequal(S.seedy, seedy)
        fprintf("\nkontrola s predoslym zadanim: max |rozdiel| PGA = %g, SGA = %g\n", ...
                max(abs(S.evoPGA(:) - evoPGA(:))), max(abs(S.evoSGA(:) - evoSGA(:))));
    else
        fprintf("\nkontrola s predoslym zadanim: ine seedy, porovnanie preskocene\n");
    end
end

% ========== grafy ==========
yAll = [evoPGA(:); evoSGA(:); evoCPGA(:)];
yLim = [min(yAll) max(yAll)] + [-1 1] * 0.03 * (max(yAll) - min(yAll)); % spolocna mierka
popisy = [nazvy "PGA najlepsi" "SGA" "CPGA" "uspesna migracia do C"];

% samostatny graf pre kazdy beh (export)
for r = 1:runs
    f = figure("Visible", "off", "Position", [100 100 1000 560]);
    ax = axes(f);
    h = kresli_beh(ax, evoPGA(:,:,r), evoSGA(r,:), evoCPGA(r,:), migLog{r}, farby, yLim, ...
                   sprintf("Beh %d (seed %d): ostrovy PGA vs CPGA vs SGA", r, seedy(r)));
    legend(ax, h, popisy, "Location", "eastoutside");
    exportgraphics(f, fullfile(vystupDir, sprintf("beh%d.png", r)), "Resolution", 300);
    close(f);
end

% suhrnne grafy (export) - hlavny graf zadania
f = figure("Visible", "off", "Position", [100 100 1000 560]);
kresli_priemer(axes(f), evoPGA, evoSGA, evoCPGA, yLim);
exportgraphics(f, fullfile(vystupDir, "priemer.png"), "Resolution", 300);
close(f);

f = figure("Visible", "off", "Position", [100 100 1000 560]);
kresli_priemer_ostrovy(axes(f), evoPGA, farby, nazvy, yLim);
exportgraphics(f, fullfile(vystupDir, "priemer_ostrovy.png"), "Resolution", 300);
close(f);

% snimky mriezky CPGA (beh 1) - sirenie dobrych rieseni / niky
f = figure("Visible", "off", "Position", [100 100 1500 330]);
kresli_mriezku(f, snapC{1}, snapGen, sprintf("CPGA beh 1 (seed %d): fitness buniek mriezky", seedy(1)));
exportgraphics(f, fullfile(vystupDir, "cpga_mriezka.png"), "Resolution", 300);
close(f);

% varianty CPGA - vlastna mierka (detail od 50. generacie), inak su rozdiely neviditelne
yV = mean(evoVar(:, min(50, ngen):end, :), 3);
yLimVar = [min(yV(:)) max(yV(:))] + [-1 1] * 0.05 * (max(yV(:)) - min(yV(:)));

f = figure("Visible", "off", "Position", [100 100 1000 560]);
kresli_varianty(axes(f), evoVar, varPar, varOk, yLimVar);
exportgraphics(f, fullfile(vystupDir, "cpga_varianty.png"), "Resolution", 300);
close(f);

% prehlad - vsetky behy + priemer v jednom okne
f = figure("Name", "Prehlad behov", "Position", [50 50 1700 900]);
tl = tiledlayout(f, 2, 3, "TileSpacing", "compact", "Padding", "compact");
for r = 1:runs
    h = kresli_beh(nexttile(tl), evoPGA(:,:,r), evoSGA(r,:), evoCPGA(r,:), migLog{r}, farby, yLim, ...
                   sprintf("Beh %d (seed %d)", r, seedy(r)));
end
lg = legend(h, popisy);
lg.Layout.Tile = "east";
kresli_priemer(nexttile(tl), evoPGA, evoSGA, evoCPGA, yLim);
exportgraphics(f, fullfile(vystupDir, "prehlad.png"), "Resolution", 200);

% varianty CPGA - samostatne okno
f = figure("Name", "Varianty CPGA", "Position", [150 150 1000 560]);
kresli_varianty(axes(f), evoVar, varPar, varOk, yLimVar);

% ========== tabulky ==========
P = [ostrov; sgaPar; nan(1, 4) cpgaPar(3:5)];
Tnast = table(["kazdy ostrov (9x)"; "SGA"; "CPGA (" + cpgaPar(1) + "x" + cpgaPar(2) + ")"], ...
              P(:,1), P(:,2), P(:,3), P(:,4), [sum(P(1:2,1:4), 2); prod(cpgaPar(1:2))], ...
              P(:,5), P(:,6), M ./ P(:,7), ...
              'VariableNames', ["populacia" "nBest" "nOld" "nCross" "nMut" "spolu" "pmutx" "pmuta" "Delta_start"]);
Tvys = tab_vysledky(evoPGA, evoSGA, evoCPGA, seedy, nazvy);
Tvek = tab_vektor(bestXP, bestFP, bestXC, bestFC, bestXS, bestFS);
Tmig = tab_migracia(migLog, nazvy);
Tvar = tab_varianty(evoVar, varPar, varOk);

fprintf("\nspolocne: nvars = %d, ngen = %d, rozsah <%d, %d>, populacia = %d, selekcia PGA/SGA %s\n", ...
        nvars, ngen, -M, M, sum(sgaPar(1:4)), func2str(sel));
fprintf("PGA: migracia kazdych %d-%d gen., reset: zlepsenie < %g za %d gen.\n", migMin, migMax, stagEps, stagK);
fprintf("CPGA: toroidna mriezka %dx%d, okolie Moore (8), partner = turnaj 2 susedov, nahrada ak lepsi, asynchronne\n", ...
        cpgaPar(1), cpgaPar(2));

fprintf("\n--- Tabulka 1: nastavenia ostrovov, SGA a CPGA ---\n"); disp(Tnast);
fprintf("--- Tabulka 2: vysledky behov (posledna generacia) ---\n"); disp(Tvys);
fprintf("--- Tabulka 3: najlepsi najdeny vektor ---\n"); disp(Tvek);
fprintf("--- Tabulka 4: migracia PGA (sucet za %d behov) ---\n", runs); disp(Tmig);
fprintf("--- Tabulka 5: varianty CPGA (posledna generacia, %d behov) ---\n", runs); disp(Tvar);
writetable(Tvar, fullfile(vystupDir, "tab5_varianty.csv"), "Delimiter", ";");

writetable(Tnast, fullfile(vystupDir, "tab1_nastavenia.csv"), "Delimiter", ";");
writetable(Tvys, fullfile(vystupDir, "tab2_vysledky.csv"), "Delimiter", ";");
writetable(Tvek, fullfile(vystupDir, "tab3_vektor.csv"), "Delimiter", ";");
writetable(Tmig, fullfile(vystupDir, "tab4_migracia.csv"), "Delimiter", ";");

save(fullfile(vystupDir, "vysledky.mat"), "evoPGA", "evoSGA", "evoCPGA", "seedy", "migLog", "snapC", ...
     "bestXP", "bestFP", "bestXS", "bestFS", "bestXC", "bestFC", "ostrov", "sgaPar", "cpgaPar", "sel", ...
     "nvars", "ngen", "runs", "Space", "migMin", "migMax", "stagK", "stagEps", "snapGen", ...
     "cpgaOkolie", "evoVar", "varPar", "varOk");

% ========== pomocne funkcie ==========

function [evo, migLog, bestX, bestF] = spusti_pga(Pops, p, sel, fitFun, Space, ngen, migMin, migMax, stagK, stagEps)
% ostrovny PGA: hviezda O -> C, reset stagnujucich vonkajsich ostrovov

nIsl = numel(Pops);
Fits = cell(nIsl, 1);
evo = nan(nIsl, ngen);
migLog = zeros(0, 3);
bestF = inf; bestX = [];
najlepsi = inf(nIsl, 1); posledne = zeros(nIsl, 1); zrod = zeros(nIsl, 1); % pre reset

dalsiaMig = inf(nIsl, 1); % casovace migracie, stred (1) nema vlastny
for i = 2:nIsl
    dalsiaMig(i) = randi([migMin migMax]);
end

for gen = 1:ngen

    % ohodnotenie ostrovov
    for i = 1:nIsl
        Fits{i} = fitFun(Pops{i});
        [evo(i,gen), iB] = min(Fits{i}); % konvergencna krivka ostrova
        if evo(i,gen) < bestF, bestF = evo(i,gen); bestX = Pops{i}(iB,:); end
    end

    % migracia podla casovacov
    for i = 2:nIsl
        if gen == dalsiaMig(i)
            [Pops, Fits, z] = do_stredu(Pops, Fits, i, gen);
            migLog(end+1,:) = z;
            dalsiaMig(i) = gen + randi([migMin migMax]); % novy casovac
        end
    end

    % reset stagnujuceho vonkajsieho ostrova (nie v poslednych 10 % generacii)
    for i = 2:nIsl
        if evo(i,gen) < najlepsi(i) - stagEps
            najlepsi(i) = evo(i,gen); posledne(i) = gen;
        elseif gen - posledne(i) >= stagK && gen < 0.9*ngen
            [Pops, Fits] = do_stredu(Pops, Fits, i, gen); % najlepsi sa pred resetom posle do C
            Pops{i} = genrpop(size(Pops{i}, 1), Space);
            Fits{i} = fitFun(Pops{i});
            najlepsi(i) = min(Fits{i}); posledne(i) = gen; zrod(i) = gen;
        end
    end

    % GA krok, krok lokalnej mutacie sa zmensuje od vzniku ostrova
    for i = 1:nIsl
        kf = max(1 - (gen - zrod(i)) / (ngen - zrod(i)), 0.05);
        Pops{i} = krok_ga(Pops{i}, Fits{i}, p, sel, Space, kf);
    end

end

end

function [Pops, Fits, zaznam] = do_stredu(Pops, Fits, i, gen)
% kopia najlepsieho z ostrova i nahradi najhorsieho jedinca v C
% zaznam = [generacia zdroj uspesna], uspesna = migrant lepsi ako doterajsi najlepsi v C

[fm, iB] = min(Fits{i});
[~, w] = max(Fits{1});
zaznam = [gen i fm < min(Fits{1})];
Pops{1}(w,:) = Pops{i}(iB,:);
Fits{1}(w) = fm;

end

function [evo, bestX, bestF] = spusti_sga(Pop, p, sel, fitFun, Space, ngen)
% jednopopulacny GA, startuje z danej populacie

evo = nan(1, ngen);
bestF = inf; bestX = [];

for gen = 1:ngen
    Fit = fitFun(Pop);
    [evo(gen), iB] = min(Fit);
    if evo(gen) < bestF, bestF = evo(gen); bestX = Pop(iB,:); end
    Pop = krok_ga(Pop, Fit, p, sel, Space, max(1 - gen/ngen, 0.05)); % rovnake zmensovanie kroku
end

end

function [evo, bestX, bestF, snap] = spusti_cpga(Pop, p, okolie, fitFun, Space, ngen, snapGen)
% celularny GA, p = [riadky stlpce pmutx pmuta dDiv], okolie = "moore" (8) / "neumann" (4)
% jedinec k lezi v bunke (riadok, stlpec) mriezky R x S (linearny index po stlpcoch)
% kazda bunka: partner = lepsi z 2 nahodnych susedov (torus),
% 1-bodove krizenie + mutacia, potomok nahradi bunku iba ak je lepsi,
% asynchronne - bunky v nahodnom poradi, zmena je hned viditelna pre susedov

R = p(1); S = p(2); N = R*S;
nvars = size(Pop, 2);
evo = nan(1, ngen);
snap = nan(R, S, numel(snapGen));
bestF = inf; bestX = [];

% tabulka susedov N x nS na toruse
[ii, jj] = ndgrid(1:R, 1:S);
ii = ii(:); jj = jj(:);
if okolie == "neumann"
    d = [-1 0; 0 -1; 0 1; 1 0];                          % von Neumann, 4 susedia
else
    d = [-1 -1; -1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0; 1 1];  % Moore, 8 susedov
end
nS = size(d, 1);
sus = zeros(N, nS);
for k = 1:nS
    sus(:,k) = sub2ind([R S], mod(ii - 1 + d(k,1), R) + 1, mod(jj - 1 + d(k,2), S) + 1);
end

Fit = fitFun(Pop); Fit = Fit(:);

for gen = 1:ngen

    % zaznam (Fit je aktualny po nahradzovani z predoslej generacie)
    [evo(gen), iB] = min(Fit);
    if evo(gen) < bestF, bestF = evo(gen); bestX = Pop(iB,:); end
    s = find(snapGen == gen, 1);
    if ~isempty(s), snap(:,:,s) = reshape(Fit, R, S); end

    % krok lokalnej mutacie sa zmensuje ako v SGA
    Delta = Space(2,:) / p(5) * max(1 - gen/ngen, 0.05);

    % nahodne cisla pre celu generaciu (vyber susedov, bod krizenia)
    n1 = randi(nS, N, 1); n2 = randi(nS, N, 1);
    bod = randi(nvars - 1, N, 1);

    % asynchronna aktualizacia - bunky v nahodnom poradi, Pop a Fit sa menia priebezne
    for k = randperm(N)

        % turnaj 2 nahodnych susedov (podla aktualneho stavu mriezky)
        a = sus(k, n1(k)); b = sus(k, n2(k));
        if Fit(b) < Fit(a), a = b; end

        % 1-bodove krizenie: zaciatok od vlastneho jedinca, zvysok od partnera
        dieta = [Pop(k, 1:bod(k)) Pop(a, bod(k)+1:end)];
        dieta = muta(mutx(dieta, p(3), Space), p(4), Delta, Space);

        % elitarna nahrada v bunke
        fd = fitFun(dieta);
        if fd < Fit(k)
            Pop(k,:) = dieta;
            Fit(k) = fd;
        end

    end

end

end

function Pop = krok_ga(Pop, Fit, p, sel, Space, kf)
% jedna generacia GA01, p = [nBest nOld nCross nMut pmutx pmuta dDiv]

Delta = Space(2,:) / p(7) * kf; % aditivny krok lokalnej mutacie

Best  = selbest(Pop, Fit, ones(1, p(1)));
Old   = selrand(Pop, Fit, p(2));
Work1 = sel(Pop, Fit, p(3));
Work2 = sel(Pop, Fit, p(4));

Work1 = crossov(Work1, 1, 0);
Work2 = muta(mutx(Work2, p(5), Space), p(6), Delta, Space);

Pop = [Best; Old; Work1; Work2];

end

function h = kresli_beh(ax, evoIsl, evoS, evoC, migLog, farby, yLim, titulok)
% jeden beh: vsetky ostrovy PGA, globalne najlepsi PGA, SGA, CPGA a uspesne migracie

hold(ax, "on"); grid(ax, "on");
[nIsl, ngen] = size(evoIsl);
h = gobjects(nIsl + 4, 1);

for i = 2:nIsl
    h(i) = plot(ax, evoIsl(i,:), "Color", farby(i,:), "LineWidth", 0.8);
end
h(1) = plot(ax, evoIsl(1,:), "Color", farby(1,:), "LineWidth", 1.8); % stred hrubsie
h(nIsl+1) = plot(ax, min(evoIsl, [], 1), "Color", [0.85 0 0], "LineWidth", 2.5);
h(nIsl+2) = plot(ax, evoS, "k", "LineWidth", 2.5);
h(nIsl+3) = plot(ax, evoC, "-.", "Color", [0 0.5 0], "LineWidth", 2.5);

% uspesne migracie - efekt je viditelny v nasledujucej generacii C
gx = min(migLog(migLog(:,3) == 1, 1) + 1, ngen);
h(nIsl+4) = plot(ax, gx, evoIsl(1,gx), "v", "LineStyle", "none", "MarkerSize", 6, ...
                 "MarkerEdgeColor", "k", "MarkerFaceColor", farby(1,:));

title(ax, titulok); xlabel(ax, "generacia"); ylabel(ax, "fitness");
xlim(ax, [1 ngen]); ylim(ax, yLim);

end

function kresli_priemer(ax, evoPGA, evoSGA, evoCPGA, yLim)
% priemer globalne najlepsieho PGA vs CPGA vs SGA, pas = +- smerodajna odchylka

[~, ngen, runs] = size(evoPGA);
pgaBest = reshape(min(evoPGA, [], 1), ngen, runs)'; % beh x generacia
mP = mean(pgaBest, 1); sP = std(pgaBest, 0, 1);
mC = mean(evoCPGA, 1); sC = std(evoCPGA, 0, 1);
mS = mean(evoSGA, 1);  sS = std(evoSGA, 0, 1);
x = 1:ngen;

hold(ax, "on"); grid(ax, "on");
fill(ax, [x fliplr(x)], [mP+sP fliplr(mP-sP)], [0.85 0 0], "FaceAlpha", 0.15, "EdgeColor", "none");
fill(ax, [x fliplr(x)], [mC+sC fliplr(mC-sC)], [0 0.5 0], "FaceAlpha", 0.15, "EdgeColor", "none");
fill(ax, [x fliplr(x)], [mS+sS fliplr(mS-sS)], [0 0 0], "FaceAlpha", 0.10, "EdgeColor", "none");
h1 = plot(ax, mP, "Color", [0.85 0 0], "LineWidth", 2.5);
h2 = plot(ax, mC, "-.", "Color", [0 0.5 0], "LineWidth", 2.5);
h3 = plot(ax, mS, "k", "LineWidth", 2.5);

title(ax, sprintf("Priemer %d behov: PGA vs CPGA vs SGA (pas = +- smer. odch.)", runs));
xlabel(ax, "generacia"); ylabel(ax, "fitness");
xlim(ax, [1 ngen]); ylim(ax, yLim);
legend(ax, [h1 h2 h3], ["PGA najlepsi (priemer)" "CPGA (priemer)" "SGA (priemer)"], "Location", "northeast");

end

function kresli_priemer_ostrovy(ax, evoPGA, farby, nazvy, yLim)
% priemer behov pre kazdy ostrov zvlast + priemer globalne najlepsieho PGA

[nIsl, ngen, runs] = size(evoPGA);
priemerOstrovy = mean(evoPGA, 3); % ostrov x generacia
pgaBest = reshape(min(evoPGA, [], 1), ngen, runs)';

hold(ax, "on"); grid(ax, "on");
h = gobjects(nIsl + 1, 1);
for i = 2:nIsl
    h(i) = plot(ax, priemerOstrovy(i,:), "--", "Color", farby(i,:), "LineWidth", 1);
end
h(1) = plot(ax, priemerOstrovy(1,:), "Color", farby(1,:), "LineWidth", 2);
h(nIsl+1) = plot(ax, mean(pgaBest, 1), "Color", [0.85 0 0], "LineWidth", 2.5);

title(ax, sprintf("Priemer %d behov: jednotlive ostrovy PGA", runs));
xlabel(ax, "generacia"); ylabel(ax, "fitness");
xlim(ax, [1 ngen]); ylim(ax, yLim);
legend(ax, h, [nazvy "PGA najlepsi"], "Location", "eastoutside");

end

function kresli_mriezku(f, snap, snapGen, titulok)
% fitness buniek CPGA mriezky vo vybranych generaciach, spolocna farebna skala

tl = tiledlayout(f, 1, numel(snapGen), "TileSpacing", "compact", "Padding", "compact");
cl = [min(snap(:)) max(snap(:))];
for s = 1:numel(snapGen)
    ax = nexttile(tl);
    imagesc(ax, snap(:,:,s), cl);
    axis(ax, "image"); ax.XTick = []; ax.YTick = [];
    title(ax, sprintf("generacia %d", snapGen(s)));
end
cb = colorbar(ax);
cb.Layout.Tile = "east";
cb.Label.String = "fitness";
title(tl, titulok);

end

function kresli_varianty(ax, evoVar, varPar, varOk, yLim)
% priemer behov pre kazdu variantu CPGA, farba = tvar mriezky, ciara = okolie
% plna = Moore, prerusovana = von Neumann, hrubsia = hlavny CPGA

[nVar, ngen, runs] = size(evoVar);
mV = mean(evoVar, 3);
mriezky = unique(varPar(:,1:2), "rows", "stable");
fMriezka = [0 0.5 0; 0.85 0.33 0.1; 0.3 0.3 0.9; 0.6 0.2 0.6];

hold(ax, "on"); grid(ax, "on");
h = gobjects(nVar, 1);
for v = 1:nVar
    [~, iM] = ismember(varPar(v,1:2), mriezky, "rows");
    styl = "-";
    if varOk(v) == "neumann", styl = "--"; end
    h(v) = plot(ax, mV(v,:), styl, "Color", fMriezka(iM,:), "LineWidth", 1.5 + (v == 1));
end

title(ax, sprintf("Priemer %d behov: varianty CPGA (detail)", runs));
xlabel(ax, "generacia"); ylabel(ax, "fitness");
xlim(ax, [1 ngen]); ylim(ax, yLim);
legend(ax, h, varPar(:,1) + "x" + varPar(:,2) + " " + varOk + " (" + round(mV(:,end)) + ")", ...
       "Location", "northeast");

end

function T = tab_varianty(evoVar, varPar, varOk)
% vysledky variant CPGA z poslednej generacie cez vsetky behy

kon = reshape(evoVar(:,end,:), size(evoVar, 1), []); % varianta x beh
T = table(varPar(:,1) + "x" + varPar(:,2), varOk, mean(kon, 2), std(kon, 0, 2), ...
          min(kon, [], 2), max(kon, [], 2), ...
          'VariableNames', ["mriezka" "okolie" "priemer" "smer_odch" "najlepsi_beh" "najhorsi_beh"]);

end

function T = tab_vysledky(evoPGA, evoSGA, evoCPGA, seedy, nazvy)
% vysledky z poslednej generacie + priemer a smerodajna odchylka

runs = size(evoSGA, 1);
[pgaKon, iIsl] = min(evoPGA(:,end,:), [], 1);
pgaKon = reshape(pgaKon, [], 1);
iIsl = reshape(iIsl, [], 1);
cpgaKon = evoCPGA(:,end);
sgaKon = evoSGA(:,end);

alg = ["PGA" "CPGA" "SGA"];
[~, iNaj] = min([pgaKon cpgaKon sgaKon], [], 2);
najlepsi = reshape(alg(iNaj), [], 1);
stat = @(v) [v; mean(v); std(v)];

T = table([string(1:runs)'; "priemer"; "smer. odch."], [seedy(:); NaN; NaN], ...
          stat(pgaKon), [reshape(nazvy(iIsl), [], 1); ""; ""], ...
          stat(cpgaKon), stat(sgaKon), [najlepsi; ""; ""], ...
          'VariableNames', ["beh" "seed" "PGA_best" "ostrov" "CPGA_best" "SGA_best" "najlepsi"]);

end

function T = tab_vektor(bestXP, bestFP, bestXC, bestFC, bestXS, bestFS)
% najlepsi najdeny vektor PGA, CPGA a SGA zo vsetkych behov

[fP, rP] = min(bestFP);
[fC, rC] = min(bestFC);
[fS, rS] = min(bestFS);
T = array2table([bestXP(rP,:); bestXC(rC,:); bestXS(rS,:)], 'VariableNames', cellstr("x" + (1:size(bestXP, 2))));
T = addvars(T, ["PGA"; "CPGA"; "SGA"], [rP; rC; rS], [fP; fC; fS], 'Before', 1, ...
            'NewVariableNames', {'algoritmus', 'beh', 'fitness'});

end

function T = tab_migracia(migLog, nazvy)
% pocet migracii a uspesnych migracii podla zdrojoveho ostrova (sucet za behy)

vsetky = vertcat(migLog{:});
zdroje = (2:numel(nazvy))';
pocet = arrayfun(@(i) sum(vsetky(:,2) == i), zdroje);
usp = arrayfun(@(i) sum(vsetky(:,2) == i & vsetky(:,3) == 1), zdroje);

T = table([reshape(nazvy(2:end), [], 1); "spolu"], [pocet; sum(pocet)], [usp; sum(usp)], ...
          round(100 * [usp ./ pocet; sum(usp) / sum(pocet)], 1), ...
          'VariableNames', ["zdroj" "migracie" "uspesne" "uspesnost_pct"]);

end