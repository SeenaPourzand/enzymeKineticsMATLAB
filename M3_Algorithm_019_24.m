function [kM, vMax, v_0] = M3_Algorithm_019_24(time, enzymeMat, subConcentrations)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 132 
% Program Description 
% In this program we ask the user for a matrix containing the tests and duplicate tests
% for a given enzyme and calculate the Vmax and Km parameters.
%
% Function Call
% [kM, vMax, v_0] = M2_Algorithm_019_24(time, enzymeMat, subConcentrations)
%
% Input Arguments
% time = a column matrix that contains the recorded time from 0-7300+ (s)
% enzymeMat = a matrix containing all the tests at different concentrations 
%     and the duplicate tests for an enzyme the user wishes to analyze
%     (μmol/L aka μM)
% subConcentrations = a row vector containing the various substrate
%     concentrations that the tests were conducted at. (units are μM)
% 
% Output Arguments
% vMax = maximum velocity reached by the data set matrix (μM/seconds)
% kM = the value that corresponds to half the max veloctiy (μM)
% v_0 = a matrix of ten initial velocities of data curve for ten concentrations (μM/seconds)
%
% Assignment Information
%   Assignment:     M3, UDF
%   Team member:    Seena Pourzand, spourzan@purdue.edu
%   Team member:    Sergio Monge, smonge@purdue.edu
%   Team member:    Greg Szymchack, gszymcha@purdue.edu
%   Team member:    Nathan Thorson, njthorso@purdue.edu
%   Team ID:        019-24
%   Academic Integrity:
%     Academic Integrity:
%     [] We worked with one or more peers but our collaboration
%        maintained academic integrity.
%			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ____________________
%% INITIALIZATION


% We removed the previous hardcoded values for substrate concentration and
% now pass them in as a parameter.


% Here we just preallocate a size to the vectors to improve performace of
% the program, continually changing the size of a vector hurts peformance.
originalV0Vec = zeros(1,10); % will store the v0s we get for the first 10 tests of an enzyme



%% ____________________
%% CALCULATIONS

% Here we changed the limit of the for loop from the previous 20 that we
% had in our M2 algorithm to 10 because the PGO data does not contain
% duplicate tests
for c = 1:1:10
    
    % We create a column vector just containing the information of current
    % index. For example at index 1 we make a column just containing the
    % product concentrations for the 3.75 substrate concentration for the
    % first regular test.
    currCol = enzymeMat(:,c);

    % We then smooth the vector we just made by using the moving mean
    % function that helps reduce noise and variation while retaining the
    % same shape and characteristics of the data.
    smoothed = movmean(currCol,10);

    % We preallocate a vector called averageSlopeSmoothed that will contain
    % the 100 slopes that we will calculate in this for loop below.
    avgSlopeSmoothed = zeros(1,100);
    % In this loop we calculate the slope between two adjacent values and
    % repeat it 100 times.
    for k = 1:1:101

        % find the slope of the next two values and add it to a vector
        avgSlopeSmoothed(k) = (smoothed(k+1) - smoothed(k)) / (time(k+1) - time(k));

    end
    
    % here we average those 100 values to end up with one final slope for
    % this column aka substrate concentration.
    originalV0Vec(c) = mean(avgSlopeSmoothed);


    
end



% Regression Portion

% Here we simply follow the basic instrucutions of the line weaver burke
% method of taking the reciporcal of our data and assigning to vectors to
% use later.
Ydata = 1 ./ originalV0Vec;
Xdata = 1 ./ subConcentrations;


% Here we use polyfit to find a lienar regression line of the lienarized
% data, giving us a value for m and b in the equation y = mx +b
linearizedCoeffs = polyfit(Xdata,Ydata,1);
linMVal = linearizedCoeffs(1);
linBVal = linearizedCoeffs(2);

%  B = 1/Vmax, so Vmax is just 1/B
vMax = 1/ linBVal;


%Using some basic algebra we can rearrange for K(m)
kM = linMVal * vMax;


% We set v_0 to the averaged v0s we determined earlier
v_0 = originalV0Vec;



%% ____________________
%% FORMATTED TEXT/FIGURE DISPLAYS

% This is all of the code we used in our development process to both
% graphically and numerically determine how well our regressions fit.
% This is not essential to the function of UDF but we decided to leave it
% here commented out if the graders wanted to understand how we test our
% fit.


mentonModel = (abs(vMax) .* subConcentrations) ./ (abs(kM) + subConcentrations);


% figure(3);
% grid on
% plot(Xdata,Ydata,"mx");
% hold on
% linModel = (Xdata .* kM .* linBVal) + linBVal;
% plot(Xdata,linModel,"-b");
% 
% figure(2);
% grid on
% plot(subConcentrations,originalV0Vec,"mx");
% hold on
% plot(subConcentrations,mentonModel,"-b");


%% ____________________
%% RESULTS
% Prior to running our function we have to run this so we allocate the
% dataset to some vectors that our function can properly utilize:
%
% raw_data = readmatrix('Data_nextGen_KEtesting_allresults.csv');
%   We input "4:end" because the first 4 rows are filled with text & NaN values 
% enzymeValues = raw_data(4:end,2:end);
% timeVector = raw_data(4:end,1);
%
% Actual function call:

% 
%% ____________________
%% ACADEMIC INTEGRITY STATEMENT
% We have not used source code obtained from any other unauthorized
% source, either modified or unmodified. Neither have we provided
% access to my code to another. The program we are submitting
% is our own original work.







    


end