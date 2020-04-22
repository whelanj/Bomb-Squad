%% Documentation

% The University of Texas at Austin - Spring 2020
% ME 337G - Nuclear Safety and Security - Dr. HAAS, Derek 
% Team 7 - Bomb Squad - INANC, Ece; Shelby; WHELAN, Jack  

% This code uses the Bateman equation to solve for the concentrations of
% the daughter products of Xe-140, a fission fragment of U-235.

%% Constants 

HL_XE = 14; % Half-life of Xe-140, [s]
HL_CS = 64; % Half-life of Cs-140, [s]
HL_BA = 13*24*60*60; % Half-life of Ba-140, [s]
HL_LA = 40*60*60; % Half-life of La-140, [s]

LL_XE = log(2)/HL_XE; % Decay constant of Xe-140, [1/s]
LL_CS = log(2)/HL_CS; % Decay constant of Cs-140, [1/s]
LL_BA = log(2)/HL_BA; % Decay constant of Ba-140, [1/s]
LL_LA = log(2)/HL_LA; % Decay constant of La-140, [1/s]
LL_CE = 0;

LL = [LL_XE LL_CS LL_BA LL_LA LL_CE];

NN = zeros(101,5);
NN(1,1) = 1;

TT = 2; % Increment time, [s]

%% Compuate Coefficients

COEFF = zeros(5);
for i = 1:4
    NUM_COEFF = i;
    for j = 1:NUM_COEFF
        NUMERATOR = 1;
        DENOMINATOR = 1;
        for k = 1:NUM_COEFF
            NUMERATOR = NUMERATOR*LL(1,k);
            if k ~= j
                DENOMINATOR = DENOMINATOR*(LL(1,k)-LL(1,j));
            end
        end
        COEFF(i,j) = NUMERATOR/DENOMINATOR;
    end
end

%% Compute Radionuclide Density

T = 0;
for i = 2:101
    T = T+TT;
    for j = 1:5
        SUMM = 0;
        for k = 1:j
            SUMM = SUMM+COEFF(j,k)*exp(-LL(1,k)*T);
        end
        NN(i,j) = NN(1,1)*SUMM/(LL(1,j));
    end
end

for i = 2:101
    NN(i,5) = NN(i-1,5)+LL(1,4)*(NN(i,4)-NN(i-1,4));
end

%% Plot 
plot(NN(:,1))
hold on
plot(NN(:,2))
hold on
plot(NN(:,3))
hold on
plot(NN(:,4))
hold on
plot(NN(:,5))
set(gca, 'YScale', 'log')
hold on
grid on
xlabel('Time, [s]')
ylabel('Nuclear Density, [# of nuclei/cm^3]')
title('Nuclear Density of Xe-140 Daughter Products')
legend('Xe-140','Cs-140','Ba-140','La-140','Ce-140')
            
  