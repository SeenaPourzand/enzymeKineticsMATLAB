function M3_main_019_24()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 132 
% Program Description 
% This Program utilizes our parameter identification program to create both
% a figure containing all of the product concentrations and their initial
% velocities as well as second figure that has the menton model determined
% from the algorithm with the v0s plotted on the same axis. The SSE is also
% computed.
% 
%
% Function Call
% M3_main_019_24()
%
% Input Arguments
%
% Output Arguments
%
% Assignment Information
%   Assignment:     M3, Part 2B
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

rawData = readmatrix("Data_PGOX50_enzyme.csv");
enzymeMat = rawData(5:end,2:end);
time = rawData(5:end,1);
subConcentrations = rawData(2,2:end);

%% ____________________
%% CALCULATIONS

[kM, vMax, v_0] = M3_Algorithm_019_24(time, enzymeMat, subConcentrations);


algoMentonModel = ((vMax) .* subConcentrations) ./ ((kM) + subConcentrations);

% SSE Calculations
SSE = sum((v_0-algoMentonModel).^2);




%% ____________________
%% FORMATTED TEXT/FIGURE DISPLAYS

fprintf("The SSE between the UDF Model and its v0s is %.4f\n",SSE);

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
plot(scaledLineVector1, (scaledLineVector1 .* v_0(1)),"-k");
title({'Product Concentration (uM) vs Time (s)', 'at Substrate Concentration(μM)'},subConcentrations(1))
legend("Product Concentration (μM)","Initial Velocity",'location','southeast');
xlabel("Time(seconds)");
ylabel("Product Concentration(μM)");

subplot(2,5,2);
plot(time,enzymeMat(:,2),"m.");
grid on
hold on
plot(scaledLineVector1, (scaledLineVector1 .* v_0(2)),"-k");
title({'Product Concentration (uM) vs Time (s)', 'at Substrate Concentration(μM)'},subConcentrations(2))
legend("Product Concentration (μM)","Initial Velocity",'location','southeast');
xlabel("Time(seconds)");
ylabel("Product Concentration(μM)");

subplot(2,5,3);
plot(time,enzymeMat(:,3),"m.");
grid on
hold on
plot(scaledLineVector1, (scaledLineVector1 .* v_0(3)),"-k");
title({'Product Concentration (uM) vs Time (s)', 'at Substrate Concentration(μM)'},subConcentrations(3))
legend("Product Concentration (μM)","Initial Velocity",'location','southeast');
xlabel("Time(seconds)");
ylabel("Product Concentration(μM)");


subplot(2,5,4);
plot(time,enzymeMat(:,4),"m.");
grid on
hold on
plot(scaledLineVector1, (scaledLineVector1 .* v_0(4)),"-k");
title({'Product Concentration (uM) vs Time (s)', 'at Substrate Concentration(μM)'},subConcentrations(4))
legend("Product Concentration (μM)","Initial Velocity",'location','southeast');
xlabel("Time(seconds)");
ylabel("Product Concentration(μM)");


subplot(2,5,5);
plot(time,enzymeMat(:,5),"m.");
grid on
hold on
plot((1:250), ((1:250) .* v_0(5)),"-k");
title({'Product Concentration (uM) vs Time (s)', 'at Substrate Concentration(μM)'},subConcentrations(5))
legend("Product Concentration (μM)","Initial Velocity",'location','southeast');
xlabel("Time(seconds)");
ylabel("Product Concentration(μM)");

subplot(2,5,6);
plot(time,enzymeMat(:,6),"m.");
grid on
hold on
plot(scaledLineVector2, (scaledLineVector2 .* v_0(6)),"-k");
title({'Product Concentration (uM) vs Time (s)', 'at Substrate Concentration(μM)'},subConcentrations(6))
legend("Product Concentration (μM)","Initial Velocity",'location','southeast');
xlabel("Time(seconds)");
ylabel("Product Concentration(μM)");

subplot(2,5,7);
plot(time,enzymeMat(:,7),"m.");
grid on
hold on
plot(scaledLineVector2, (scaledLineVector2 .* v_0(7)),"-k");
title({'Product Concentration (uM) vs Time (s)', 'at Substrate Concentration(μM)'},subConcentrations(7))
legend("Product Concentration (μM)","Initial Velocity",'location','southeast');
xlabel("Time(seconds)");
ylabel("Product Concentration(μM)");


subplot(2,5,8);
plot(time,enzymeMat(:,8),"m.");
grid on
hold on
plot((1:500), ((1:500) .* v_0(8)),"-k");
title({'Product Concentration (uM) vs Time (s)', 'at Substrate Concentration(μM)'},subConcentrations(8))
legend("Product Concentration (μM)","Initial Velocity",'location','southeast');
xlabel("Time(seconds)");
ylabel("Product Concentration(μM)");


subplot(2,5,9);
plot(time,enzymeMat(:,9),"m.");
grid on
hold on
plot((1:600), ((1:600) .* v_0(9)),"-k");
title({'Product Concentration (uM) vs Time (s)', 'at Substrate Concentration(μM)'},subConcentrations(9))
legend("Product Concentration (μM)","Initial Velocity",'location','southeast');
xlabel("Time(seconds)");
ylabel("Product Concentration(μM)");


subplot(2,5,10);
plot(time,enzymeMat(:,10),"m.");
grid on
hold on
plot(scaledLineVector3, (scaledLineVector3 .* v_0(10)),"-k");
title({'Product Concentration (uM) vs Time (s)', 'at Substrate Concentration(μM)'},subConcentrations(10))
legend("Product Concentration (μM)","Initial Velocity",'location','southeast');
xlabel("Time(seconds)");
ylabel("Product Concentration(μM)");


figure(2)
grid on
plot(subConcentrations,algoMentonModel,"-k");
hold on
plot(subConcentrations,v_0,"rx");
title("UDF Michaelis-Menten(μM/s) & UDF-Determined v0s(μM/s) vs Substrate Concentration(μM)");
legend("UDF MM Reaction Velocities(μM/s)","UDF-Determined v0s(μM/s)",'location','southeast');
xlabel("Substrate Concentration(μM)");
ylabel("Reaction Velocity(μM/s)");


%% ____________________
%% RESULTS
% M3_main_019_24()

% M3 does not output anything as it has no outputs. It does however create
% two figures and computes the SSE for the parameter identification algorithm

%% ____________________
%% ACADEMIC INTEGRITY STATEMENT
% We have not used source code obtained from any other unauthorized
% source, either modified or unmodified. Neither have we provided
% access to my code to another. The program we are submitting
% is our own original work.


end
