function M3_MM_PGOX50_019_24(time, enzymeMat, subConcentrations, refValues)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 132 
% Program Description 
% This program plots the v0s against their respective product concentration
% graphs and analyze the fit of the michalis menton plot determeined.
%
% Function Call
% M3_MM_PGOX50_019_24(time, enzymeMat, subConcentrations, refValues)
%
% Input Arguments
% time = a column matrix that contains the recorded time from 0-7300+ (s)
% enzymeMat = a matrix containing all the tests at different concentrations 
%     and the duplicate tests for an enzyme the user wishes to analyze
%     (μmol/L aka μM)
% subConcentrations = a row vector containing the various substrate
%     concentrations that the tests were conducted at. (units are μM)
% refValues = the PGOX50 reference values for the 10 v0s & Vmax (μM/seconds) &
% the Km(μM)
%
% Output Arguments
%
% Assignment Information
%   Assignment:     M3, Part 2A
%   Team member:    Seena Pourzand, spourzan@purdue.edu
%   Team member:    Sergio Monge, smonge@purdue.edu
%   Team member:    Greg Szymchack, gszymcha@purdue.edu
%   Team member:    Nathan Thorson, njthorso@purdue.edu
%   Team ID:        019-24
%   Academic Integrity:
%     Academic Integrity:
%     [] We worked with one or more peers but our collaboration
%        maintained academic integrity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ____________________
%% INITIALIZATION

%rawData = readmatrix("Data_PGOX50_enzyme.csv");
%enzymeMat = rawData(5:end,2:end);
%enzymeMat = rawData(5:end,2:end);
%subConcentrations = rawData(2,2:end);
%refValues = [0.025 0.049 0.099 0.176 0.329 0.563 0.874 1.192 1.361 1.603 1.806 269.74];

% initializing the reference parameters by retrieving it from the passed in
% vector
refVmax = refValues(11);
refKm = refValues(12);

%% ____________________
%% CALCULATIONS

% calculating the menton model by utilizing the reference parameters and
% substrate concentrations
refMentonModel = ((refVmax) .* subConcentrations) ./ ((refKm) + subConcentrations);

% SSE Calculations
SSE = sum((refValues(1:10)-refMentonModel).^2);



%% ____________________
%% FORMATTED TEXT/FIGURE DISPLAYS
fprintf("The SSE between the Reference Model and the provided reference v0s is %.4f\n",SSE);


% Here we create two vectors that we will use to plot our tangent/v0 lines,
% acting as the corresponding x values and we have them be seperate from
% the time vector so can limit the length of the tangent line and how it appears on the graph. 
scaledLineVector1 = 1:200;
scaledLineVector2 = 1:300;
scaledLineVector3 = 1:800;

figure(1);

subplot(2,5,1);
plot(time,enzymeMat(:,1),"m.");
grid on
hold on
plot(scaledLineVector1, (scaledLineVector1 .* refValues(1)),"-k");
title({'Product Concentration (uM) vs Time (s)', 'at Substrate Concentration(μM)'},subConcentrations(1))
legend("Product Concentration (μM)","Initial Velocity",'location','southeast');
xlabel("Time(seconds)");
ylabel("Product Concentration(μM)");

subplot(2,5,2);
plot(time,enzymeMat(:,2),"m.");
grid on
hold on
plot(scaledLineVector1, (scaledLineVector1 .* refValues(2)),"-k");
title({'Product Concentration (uM) vs Time (s)', 'at Substrate Concentration(μM)'},subConcentrations(2))
legend("Product Concentration (μM)","Initial Velocity",'location','southeast');
xlabel("Time(seconds)");
ylabel("Product Concentration(μM)");

subplot(2,5,3);
plot(time,enzymeMat(:,3),"m.");
grid on
hold on
plot(scaledLineVector1, (scaledLineVector1 .* refValues(3)),"-k");
title({'Product Concentration (uM) vs Time (s)', 'at Substrate Concentration(μM)'},subConcentrations(3))
legend("Product Concentration (μM)","Initial Velocity",'location','southeast');
xlabel("Time(seconds)");
ylabel("Product Concentration(μM)");


subplot(2,5,4);
plot(time,enzymeMat(:,4),"m.");
grid on
hold on
plot(scaledLineVector1, (scaledLineVector1 .* refValues(4)),"-k");
title({'Product Concentration (uM) vs Time (s)', 'at Substrate Concentration(μM)'},subConcentrations(4))
legend("Product Concentration (μM)","Initial Velocity",'location','southeast');
xlabel("Time(seconds)");
ylabel("Product Concentration(μM)");


subplot(2,5,5);
plot(time,enzymeMat(:,5),"m.");
grid on
hold on
plot((1:250), ((1:250) .* refValues(5)),"-k");
title({'Product Concentration (uM) vs Time (s)', 'at Substrate Concentration(μM)'},subConcentrations(5))
legend("Product Concentration (μM)","Initial Velocity",'location','southeast');
xlabel("Time(seconds)");
ylabel("Product Concentration(μM)");


subplot(2,5,6);
plot(time,enzymeMat(:,6),"m.");
grid on
hold on
plot(scaledLineVector2, (scaledLineVector2 .* refValues(6)),"-k");
title({'Product Concentration (uM) vs Time (s)', 'at Substrate Concentration(μM)'},subConcentrations(6))
legend("Product Concentration (μM)","Initial Velocity",'location','southeast');
xlabel("Time(seconds)");
ylabel("Product Concentration(μM)");

subplot(2,5,7);
plot(time,enzymeMat(:,7),"m.");
grid on
hold on
plot(scaledLineVector2, (scaledLineVector2 .* refValues(7)),"-k");
title({'Product Concentration (uM) vs Time (s)', 'at Substrate Concentration(μM)'},subConcentrations(7))
legend("Product Concentration (μM)","Initial Velocity",'location','southeast');
xlabel("Time(seconds)");
ylabel("Product Concentration(μM)");


subplot(2,5,8);
plot(time,enzymeMat(:,8),"m.");
grid on
hold on
plot((1:500), ((1:500) .* refValues(8)),"-k");
title({'Product Concentration (uM) vs Time (s)', 'at Substrate Concentration(μM)'},subConcentrations(8))
legend("Product Concentration (μM)","Initial Velocity",'location','southeast');
xlabel("Time(seconds)");
ylabel("Product Concentration(μM)");


subplot(2,5,9);
plot(time,enzymeMat(:,9),"m.");
grid on
hold on
plot((1:600), ((1:600) .* refValues(9)),"-k");
title({'Product Concentration (uM) vs Time (s)', 'at Substrate Concentration(μM)'},subConcentrations(9))
legend("Product Concentration (μM)","Initial Velocity",'location','southeast');
xlabel("Time(seconds)");
ylabel("Product Concentration(μM)");


subplot(2,5,10);
plot(time,enzymeMat(:,10),"m.");
grid on
hold on
plot(scaledLineVector3, (scaledLineVector3 .* refValues(10)),"-k");
title({'Product Concentration (uM) vs Time (s)', 'at Substrate Concentration(μM)'},subConcentrations(10))
legend("Product Concentration (μM)","Initial Velocity",'location','southeast');
xlabel("Time(seconds)");
ylabel("Product Concentration(μM)");


figure(2)
grid on
plot(subConcentrations,refMentonModel,"-k");
hold on
plot(subConcentrations,refValues(1:10),"rx");
title("Reference Michaelis-Menten(μM/s) & Reference v0s(μM/s) vs Substrate Concentration(μM)");
legend("Reference MM Reaction Velocities(μM/s)","Reference v0s(μM/s)",'location','southeast');
xlabel("Substrate Concentration(μM)");
ylabel("Reaction Velocity(μM/s)");


%% ____________________
%% RESULTS

% No outputs to the command window or workspace, but SSE is computed and
% can be displayed if the user decides and 2 figures are generated


%% ____________________
%% ACADEMIC INTEGRITY STATEMENT
% We have not used source code obtained from any other unauthorized
% source, either modified or unmodified. Neither have we provided
% access to my code to another. The program we are submitting
% is our own original work.


end
