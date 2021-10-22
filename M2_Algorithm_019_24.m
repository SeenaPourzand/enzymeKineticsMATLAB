function [kM, vMax, v_0] = M2_Algorithm_019_24(time, enzymeMat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 132 
% Program Description 
% In this program we ask the user for a matrix containing the tests and duplicate tests
% for a given enzyme and calculate the Vmax and Km parameters.
%
% Function Call
% [kM, vMax, v_0] = M2_Algorithm_019_24(time, enzymeMat)
%
% Input Arguments
% time = a column matrix that contains the recorded time from 0-7300+ (s)
% enzymeMat = a matrix containing all the tests at different concentrations 
%     and the duplicate tests for an enzyme the user wishes to analyze
%     (μmol/L aka μM)
% 
% Output Arguments
% vMax = maximum velocity reached by the data set matrix (μM/seconds)
% kM = the value that corresponds to half the max veloctiy (μM)
% v_0 = a matrix of ten initial velocities of data curve for ten concentrations (μM/seconds)
%
% Assignment Information
%   Assignment:     M2, UDF
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


% Here we just initialize the substrate concentrations. Since all the
% enzymes and their tests share the same increments in concentration we can
% just initialize it once without having to pass it in. Units are μmol/L
% aka μM
subConcentrations = [3.75 7.5 15 30 65 125 250 500 1000 2000];


% Here we just preallocate a size to the vectors to improve performace of
% the program, continually changing the size of a vector hurts peformance.
originalV0Vec = zeros(1,10); % will store the v0s we get for the first 10 tests of an enzyme
duplicateV0Vec = zeros(1,10);% will store the v0s we get for the 10 duplicate tests of an enzyme
% unts in μM/seconds


%% ____________________
%% CALCULATIONS

% Here we have a for loop that will iterate 20 times, for the number of
% total tests we have avaliable for each enzyme, 10 regular + 10 duplicate
for c = 1:1:20
    
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
    slopeOneSmooth  = mean(avgSlopeSmoothed);


    % This selection structure simply allocates where our values will go.
    % If the index is less or equal to ten then we know it is one of the ten
    % normal tests. If it is greater than ten , we know it is one of the
    % ten duplicate tests
    if (c > 10)
        duplicateV0Vec(c-10) = slopeOneSmooth;
    else
        originalV0Vec(c) = slopeOneSmooth;
        
    end
end

% Here we preallocate an array so we can eventually store our v0s in here.
averagedV0s = zeros(1,10);

for x = 1:1:10
    
    %Here we take the averages of the corresponding v0s that share the same
    %concentration so we have a final vector of 10 v0s.
    averagedV0s(x) = (originalV0Vec(x) + duplicateV0Vec(x)) / 2;
end



% Regression Portion

% Here we simply follow the basic instrucutions of the line weaver burke
% method of taking the reciporcal of our data and assigning to vectors to
% use later.
Ydata = 1 ./ averagedV0s;
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
v_0 = averagedV0s;



%% ____________________
%% FORMATTED TEXT/FIGURE DISPLAYS

% This is all of the code we used in our development process to both
% graphically and numerically determine how well our regressions fit.
% This is not essential to the function of UDF but we decided to leave it
% here commented out if the graders wanted to understand how we test our
% fit.


mentonModel = (abs(vMax) .* subConcentrations) ./ (abs(kM) + subConcentrations);

% Calculate SSE of non-linear model
SSE = sum((averagedV0s-mentonModel).^2);

% Calculate SST of non-linear model
SST = sum((averagedV0s-mean(averagedV0s)).^2);

rSquared = 1 - (SSE/SST);
% 
% figure(1);
% grid on
% plot(Xdata,Ydata,"mx");
% hold on
% linModel = (Xdata .* kM .* linBVal) + linBVal;
% plot(Xdata,linModel,"-b");
% 
% figure(2);
% grid on
% plot(subConcentrations,averagedV0s,"mx");
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

% [kM, vMax, v_0] = M2_Algorithm_019_24(timeVector, enzymeValues(:,1:20))
% rSquared = 0.9993
% kM = 188.6738 (μM)
% vMax = 0.9280 (μM/seconds)
% v_0 = 0.0181    0.0354    0.0668    0.1308    0.2367    0.3746    0.5343  0.6771    0.7973    0.8649 (μM/seconds)
%
%
% [kM, vMax, v_0] = M2_Algorithm_019_24(timeVector, enzymeValues(:,21:40))
% rSquared = 0.9971
% kM = 371.6624 (μM)
% vMax = 0.8471 (μM/seconds)
% v_0 = 0.0085    0.0167    0.0332    0.0634    0.1239    0.2092    0.3393  0.4822    0.6596    0.7045 (μM/seconds)


% [kM, vMax, v_0] = M2_Algorithm_019_24(timeVector, enzymeValues(:,41:60))
% rSquared = 0.9798
% kM = 283.1645 (μM)
% vMax = 1.3836 (μM/seconds)
% v_0 = 0.0180    0.0363    0.0699    0.1341    0.2600    0.4163    0.6263   0.8186    0.9780    1.0882 (μM/seconds)

% [kM, vMax, v_0] = M2_Algorithm_019_24(timeVector, enzymeValues(:,61:80))
% rSquared = 0.9934
% kM = 340.4079 (μM)
% vMax = 1.4948 (μM/seconds)
% v_0 = 0.0163    0.0318    0.0626    0.1215    0.2402    0.4108    0.6567  0.9290    1.1949    1.3586 (μM/seconds)


% [kM, vMax, v_0] = M2_Algorithm_019_24(timeVector, enzymeValues(:,81:100))
% rSquared = 0.9995
% kM = 248.5041 (μM)
% vMax = 1.6485 (μM/seconds)
% v_0 = 0.0245    0.0480    0.0957    0.1701    0.3444    0.5644    0.8433  1.1177    1.3388    1.4853 (μM/seconds)

%% ____________________
%% ACADEMIC INTEGRITY STATEMENT
% We have not used source code obtained from any other unauthorized
% source, either modified or unmodified. Neither have we provided
% access to my code to another. The program we are submitting
% is our own original work.







    


end