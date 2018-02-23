close all

global CI_PROT CI_RNA CRO_PROT CRO_RNA
global CI CRO
CI_PROT = 1;
CI_RNA = 2;
CRO_PROT = 3;
CRO_RNA = 4;
CI = 1;
CRO = 2;
[chiCIProt, chiCIRna] = deal(1.2);
[chiCroProt, chiCroRna] = deal(0.8);
[omegaCI, muCI, omegaCro, muCro] = deal(50);
[kCI, kCro] = deal(10);
dt = 0.01;
totalTime = 25;

% Part 1
figureNum = 1
figure(figureNum)
initCons = zeros(3, 4);
initCons(1,:) = [0, 0, 0, 0];
initCons(2,:) = [0, 0, 0, 20];
initCons(3,:) = [0, 50, 0, 0];
titles = ["Initial concentrations = 0", "cro RNA = 20 molecules", "cI RNA = 50 molecules"];
for i = 1 : size(initCons, 1)
    [time, cIProt, cIRna, croProt, croRna] = findConcentration(totalTime, dt, ...
    initCons(i, :), [chiCIProt, chiCIRna, chiCroProt chiCroRna], [omegaCI, omegaCro], ...
    [muCI, muCro], [kCI, kCro]);


    subplot(2, 2, i);
    hold on
    plot(time, cIProt);
    plot(time, cIRna);
    plot(time, croProt);
    plot(time, croRna);
    hold off
    molVsTimeSettings(1, titles(i)) 
end

% Part 2
figureNum = figureNum + 1;
varCons = (0 : 1 : 20);
varCons1 = (0 : 500 : 2000);
for k = 1 : 2
    if (k == 1); cons = varCons; else; cons = varCons1; end
    figure(figureNum);
    hold on
    for i = 1 : size(cons, 2) % CI_RNA
        for j = 1 : size(cons, 2) % CRO_RNA
            initCons = zeros(1, 4);
            initCons(CI_RNA) = cons(i);
            initCons(CRO_RNA) = cons(j);
            [time, cIProt, cIRna, croProt, croRna] = findConcentration(totalTime, dt, initCons, ...
                [chiCIProt, chiCIRna, chiCroProt chiCroRna], [omegaCI, omegaCro], ...
                [muCI, muCro], [kCI, kCro]);
                plot(cIProt, croProt);
        end
    end
    figureNum = figureNum + 1;
end



function [time, cIProt, cIRna, croProt, croRna] = findConcentration(totalTime, dt, initCons, ...
    chis, omegas, mus, ks)
    CI_PROT = 1;
    CI_RNA = 2;
    CRO_PROT = 3;
    CRO_RNA = 4;
    CI = 1;
    CRO = 2;
    [cIProt, cIRna, croProt, croRna] = deal(zeros(1, totalTime/dt));
    cIProt(1) = initCons(CI_PROT);
    cIRna(1) = initCons(CI_RNA);
    croProt(1) = initCons(CRO_PROT);
    croRna(1) = initCons(CRO_RNA);
    time = zeros(1, totalTime/dt);
    % Forward Euler Algorithm
    for i = 1 : totalTime/dt    
        cIProt(i+1) =  cIProt(i) + (dtCIProt(omegas(CI), cIProt(i), cIRna(i), chis(CI_PROT)) * dt);
        cIRna(i+1) =  cIRna(i) + (dtCIRna(mus(CI), croProt(i), cIRna(i), ks(CRO), chis(CI_RNA)) * dt);
        croProt(i+1) = croProt(i) + (dtCroProt(omegas(CRO), croProt(i), croRna(i), chis(CRO_PROT)) * dt);
        croRna(i+1) = croRna(i) + (dtCroRna(mus(CRO), cIProt(i), croRna(i), ks(CI), chis(CRO_RNA)) * dt);
        time(i+1) = time(i) + dt;
    end
end



function y = dtCIProt(omegaCI, cIProt, cIRna, chiCIProt)
    y = (omegaCI * cIRna) - (chiCIProt * cIProt);
end

function y = dtCIRna(muCI, croProt, cIRna, kCro, chiCIRna)
    gen = muCI * (1 - (croProt .^ 2)/(kCro .^ 2 + croProt .^ 2));
    deg = chiCIRna * cIRna;
    y = gen - deg;
end

function y = dtCroProt(omegaCro, croProt, croRna, chiCroProt)
    y = (omegaCro * croRna) - (chiCroProt * croProt);
end

function y = dtCroRna(muCro, cIProt, croRna, kCI, chiCroRna)
    gen = muCro * (1 - (cIProt .^ 2)/(kCI .^ 2 + cIProt .^ 2));
    deg = chiCroRna * croRna;
    y = gen - deg;
end

function molVsTimeSettings(figureNum, t)
    X_LABEL = "Time (s)";
    Y_LABEL = "Molecules";
    figure(figureNum);
    xlabel(X_LABEL);
    ylabel(Y_LABEL);
    labels = ["cI Protein", "cI RNA", "cro Protein", "cro RNA"];
    legend(labels);
    title(t);
end

