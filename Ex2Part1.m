close all

global CI CRO
CI_PROT = 1;
CI_RNA = 2;
CRO_PROT = 3;
CRO_RNA = 4;
CI = 1;
CRO = 2;


% Part 1.
figureNum = 1;
[chiCIProt, chiCIRna] = deal(1.2);
[chiCroProt, chiCroRna] = deal(0.8);
[omegaCI, muCI, omegaCro, muCro] = deal(50);
[kCI, kCro] = deal(10);
steps = 50000;
[initCIProt, initCIRna, initCroProt, initCroRna] = deal(0);
[time, cIProt, cIRna, croProt, croRna] = simulate(initCIProt, initCIRna, initCroProt, ...
     initCroRna, steps);
hold on
plot(time, cIProt);
plot(time, cIRna);
plot(time, croProt);
plot(time, croRna);
molVsTimeSettings(figureNum, "Exercise 2, Question 1")
hold off

% Part 2.
figureNum = figureNum + 1;
figure(figureNum);
% 20 iterations of part 1.
it = 20;
hold on
for i = 1 : it
    [time, cIProt, cIRna, croProt, croRna] = simulate(initCIProt, initCIRna, initCroProt, ...
     initCroRna, steps);
    plot(time, cIProt);
    plot(time, cIRna);
    plot(time, croProt);
    plot(time, croRna);
end
hold off 
molVsTimeSettings(figureNum, "Exercise 1, Question 2")

% Phase plane.
figureNum = figureNum + 1;
figure(figureNum);
hold on
for i = 1 : it
    [time, cIProt, cIRna, croProt, croRna] = simulate(initCIProt, initCIRna, initCroProt, ...
     initCroRna, steps);
    plot(croProt, cIProt);
end
hold off 

% Part 3
figureNum = figureNum + 1;
figure(figureNum);
initCons = deal(zeros(2, 4));
initCons(1,:) = [0, 20, 0, 0]; % cIRna = 20;
initCons(2,:) = [0, 0, 0, 20]; % croRna = 20;
titles = ["cIRna = 20", "croRna = 20"];
for i = 1 : size(initCons, 1)
    figure(figureNum);
    hold on
    cons = initCons(i,:);
    for j = 1 : it
        [time, cIProt, cIRna, croProt, croRna] = simulate(cons(1), cons(2), cons(2), ...
         cons(3), steps);
        plot(time, cIProt);
        plot(time, cIRna);
        plot(time, croProt);
        plot(time, croRna);
    end
    molVsTimeSettings(figureNum, titles(i))
    hold off 
    figureNum = figureNum + 1;
end









function [time, cIProt, cIRna, croProt, croRna] = simulate(initCIProt, initCIRna, initCroProt, ...
     initCroRna, steps)
    [chiCIProt, chiCIRna] = deal(1.2);
    [chiCroProt, chiCroRna] = deal(0.8);
    [omegaCI, muCI, omegaCro, muCro] = deal(50);
    [kCI, kCro] = deal(10);
    [cIProt, cIRna, croProt, croRna] = deal(zeros(1, steps));
    cIProt(1) = initCIProt;
    cIRna(1) = initCIRna;
    croProt(1) = initCroProt;
    croRna(1) = initCroRna;
    time = zeros(1, steps);

    for i = 1 : steps
       % Calculate time step.
       numReactions = 8;
       reactions = zeros(1, numReactions);
       reactions(1) = v1(omegaCI, cIRna(i));
       reactions(2) = v2(cIProt(i), chiCIProt);
       reactions(3) = v3(muCI, croProt(i), kCro);
       reactions(4) = v4(cIRna(i), chiCIRna);
       reactions(5) = v5(omegaCro, croRna(i));
       reactions(6) = v6(croProt(i), chiCroProt);
       reactions(7) = v7(muCro, cIProt(i), kCI);
       reactions(8) = v8(croRna(i), chiCroRna);

       rTot = sum(reactions);
       tStep = exprnd(1 / rTot);
       time(i+1) = time(i) + tStep;
       % Choose reaction.
       prob = rand() * rTot;
       nextReaction = 1;
       for j = 1 : numReactions
           if (prob < sum(reactions(1:j)))
               nextReaction = j;
               break
           end
       end
       % Execute reaction.
       % Update values initially.   
       cIProt(i+1) = cIProt(i);
       cIRna(i+1) = cIRna(i);
       croProt(i+1) = croProt(i);     
       croRna(i+1) = croRna(i);
       switch nextReaction
           case 1
               cIProt(i+1) = cIProt(i) + 1;
           case 2
               cIProt(i+1) = cIProt(i) - 1;
           case 3
               cIRna(i+1) = cIRna(i) + 1;
           case 4
               cIRna(i+1) = cIRna(i) - 1;
           case 5
               croProt(i+1) = croProt(i) + 1;
           case 6
               croProt(i+1) = croProt(i) - 1;
           case 7
               croRna(i+1) = croRna(i) + 1;
           case 8
               croRna(i+1) = croRna(i) - 1;
           otherwise
               disp("SOMETHING WENT WRONG!!!");
       end
    end
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

function gen = v1(omegaCI, cIRna)
    gen = (omegaCI * cIRna);
end

function deg = v2(cIProt, chiCIProt)
    deg = (chiCIProt * cIProt);
end

function gen = v3(muCI, croProt, kCro)
    gen = muCI * (1 - (croProt .^ 2)/(kCro .^ 2 + croProt .^ 2));
end

function deg = v4(cIRna, chiCIRna)
    deg = chiCIRna * cIRna;
end

function gen = v5(omegaCro, croRna)
    gen = (omegaCro * croRna);
end

function deg = v6(croProt, chiCroProt)
    deg = (chiCroProt * croProt);
end

function gen = v7(muCro, cIProt, kCI)
    gen = muCro * (1 - (cIProt .^ 2)/(kCI .^ 2 + cIProt .^ 2));
end

function deg = v8(croRna, chiCroRna)
    deg = chiCroRna * croRna;
end

