clc;
clear;
close all;

% --- path bootstrap ---
thisFile = mfilename("fullpath");
thisDir = fileparts(thisFile);
addpath(fullfile(thisDir, "genetic"));
% -----------------------------------------------

% nastavenia, parametre:
nvars = 10; % pocet premennych
popsize = 50; % velkost populacie
ngen = 500; % pocet generacii
runs = 5; % pocet behov na jednu variantu
M = 500; % ohranicenie priestoru
Space = [ones(1,nvars)*(-M); ones(1,nvars)*M];

fitFun = @testfn3;
doBonus = true; % bonusova uloha - eggholder

% varianty a-i, kazdy riadok je jedna parametrizacia
% tlak: seltourn (velky) > selsus (stredny) > selrand (ziadny)
%       viac nBest tlak zvysuje, nOld (nahodni jedinci) ho znizuje
% diverzita: pmutx, pmuta a velkost kroku Delta = Space(2,:)/dDiv
% varianta e ma hodnoty priamo z GA01, f-i su e s vypnutym operatorom

varianty = ["a" "b" "c" "d" "e" "f" "g" "h" "i"];

popisy = ["a) velky tlak, velka diverzita"
          "b) velky tlak, mala diverzita"
          "c) maly tlak, velka diverzita"
          "d) maly tlak, mala diverzita"
          "e) kompromis (GA01)"
          "f) e bez krizenia"
          "g) e bez mutx"
          "h) e bez muta"
          "i) e bez oboch mutacii"];

%          -- pocty jedincov (spolu 50) --   ---- mutacia ----    -- zapnute --
%          nBest  nOld  nCross  nMut         pmutx pmuta dDiv     cross mutx muta
cisla = [    4     0     20     26           0.25  0.40   10       1     1    1    % a) velky tlak, velka div.
             4     0     40      6           0.05  0.02  500       1     1    1    % b) velky tlak, mala div.
             1     9     20     20           0.25  0.40   10       1     1    1    % c) maly tlak, velka div.
             1     9     34      6           0.05  0.02  500       1     1   1    % d) maly tlak, mala div.
             2    10     18     20           0.15  0.15  100       1     1    1    % e) GA01
             2    10     18     20           0.15  0.15  100       0     1    1    % f) bez krizenia
             2    10     18     20           0.15  0.15  100       1     0    1    % g) bez mutx
             2    10     18     20           0.15  0.15  100       1     1    0    % h) bez muta
             2    10     18     20           0.15  0.15  100       1     0    0 ]; % i) bez mutacii

% selekcia pre pracovne skupiny (Work1, Work2)
selekcie = ["seltourn" "seltourn" "selrand" "selrand" "selsus" ...
            "selsus" "selsus" "selsus" "selsus"];

% ===== SCHWEFEL =====
% ulozisko vysledkov, evo{k}(beh, generacia)
evo = cell(1, numel(varianty));

for k = 1:numel(varianty)

    fprintf("%s\n", popisy(k));
    evo{k} = spusti_ga(cisla(k,:), selekcie(k), fitFun, Space, popsize, ngen, runs);

end

% ULOHA 1 - varianty a-e do jedneho obrazku
kresli(evo, popisy, varianty, 1:5, "Uloha 1: a-e (Schwefel)", ...
       fullfile(thisDir, "uloha1.png"));

% ULOHA 2 - varianty e-i do druheho obrazku
kresli(evo, popisy, varianty, 5:9, "Uloha 2: e-i (Schwefel)", ...
       fullfile(thisDir, "uloha2.png"));

% ULOHA 3 - cisla k diskusii (globalne minimum Schwefel = -418.98*10 = -4189.8)
fprintf("\nSCHWEFEL\n");
vypis(evo, popisy);

% ===== BONUS - EGGHOLDER, to iste pre a-i =====
if doBonus

    SpaceE = [ones(1,nvars)*(-512); ones(1,nvars)*512];
    evoE = cell(1, numel(varianty));

    for k = 1:numel(varianty)

        fprintf("egg %s\n", popisy(k));
        evoE{k} = spusti_ga(cisla(k,:), selekcie(k), @eggholder, SpaceE, popsize, ngen, runs);

    end

    kresli(evoE, popisy, varianty, 1:5, "Bonus 1: a-e (Eggholder)", ...
           fullfile(thisDir, "bonus1.png"));

    kresli(evoE, popisy, varianty, 5:9, "Bonus 2: e-i (Eggholder)", ...
           fullfile(thisDir, "bonus2.png"));

    fprintf("\nEGGHOLDER\n");
    vypis(evoE, popisy);

end

% ========== pomocne funkcie ==========

function evo = spusti_ga(p, selName, fitFun, Space, popsize, ngen, runs)
% GA01 prepisany ako funkcia, cisla sa beru z jedneho riadku "cisla"

% p = [nBest nOld nCross nMut pmutx pmuta dDiv cross mutx muta]
nBest = p(1);
nOld = p(2);
nCross = p(3);
nMut = p(4);
pmutx = p(5);
pmuta = p(6);
Delta = Space(2,:) / p(7); % aditivny krok lokalnej mutacie (GA01: /100)
useCross = p(8);
useMutx = p(9);
useMuta = p(10);

evo = nan(runs, ngen);

for r = 1:runs

    Pop = genrpop(popsize, Space); % inicializacia populacie

    for gen = 1:ngen

        Fit = fitFun(Pop);
        evo(r,gen) = min(Fit); % konvergencna krivka

        % GA
        Best = selbest(Pop, Fit, ones(1,nBest));

        if nOld > 0
            Old = selrand(Pop, Fit, nOld);
        else
            Old = [];
        end

        Work1 = vyber(selName, Pop, Fit, nCross);
        Work2 = vyber(selName, Pop, Fit, nMut);

        if useCross
            Work1 = crossov(Work1, 1, 0);
        end

        if useMutx
            Work2 = mutx(Work2, pmutx, Space);
        end

        if useMuta
            Work2 = muta(Work2, pmuta, Delta, Space);
        end

        Pop = [Best; Old; Work1; Work2];

    end

    fprintf("   beh %d/%d: bestF = %10.2f\n", r, runs, evo(r,end));

end

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

function kresli(evo, popisy, varianty, idx, nazov, subor)
% 5 grafov (behy + priemer) a 6. graf s priemermi, ulozi ako PNG pre Word

f = figure("Name", nazov, "Position", [100 100 1200 700]);

for j = 1:numel(idx)

    k = idx(j);

    subplot(2,3,j);
    hold on;
    grid on;

    plot(evo{k}', "LineWidth", 1); % jednotlive behy
    plot(mean(evo{k},1), "k", "LineWidth", 2); % priemer

    title(popisy(k));
    xlabel("generacia");
    ylabel("fitness");
    axis tight;

end

subplot(2,3,6);
hold on;
grid on;

for k = idx
    plot(mean(evo{k},1), "LineWidth", 1.5);
end

title("priemery");
xlabel("generacia");
ylabel("fitness");
axis tight;
legend(varianty(idx));

sgtitle(nazov);
exportgraphics(f, subor, "Resolution", 300);

end

function vypis(evo, popisy)
% tabulka vysledkov z poslednej generacie

fprintf("%-32s %10s %10s %10s\n", "varianta", "najlepsi", "priemer", "najhorsi");

for k = 1:numel(popisy)

    konec = evo{k}(:,end);
    fprintf("%-32s %10.2f %10.2f %10.2f\n", popisy(k), min(konec), mean(konec), max(konec));

end

end