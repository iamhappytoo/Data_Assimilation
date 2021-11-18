function [Z]=measurement_model(Y)
    H=zeros(5,185);
    H(1,8)=1;  %(19,3)
    H(2,70)=1;  %(15,8)
    H(3,94)=1; %(8,10)
    H(4,173)=1; %(12,16)
    H(5,182)=1;  %(14,17)
    
    % This example assumes spatial map of state #1 is measured; this 
    % function will vary depending on the particular measurement being 
    % assimilated, other examples could include 1) a spatial aggregation of
    % the state map to indicate a measurement coarser than the modeled
    % states, 2) states sampled at specific locations i.e. to indicate in
    % situ data being assimilated, 3) A function (i.e. radiative transfer
    % equation) representing that the measurement space is different than
    % the state space, etc.
    %% if true generation, do
    %Z=H * Y;
    %% if for Enkf, do 
    Z=H * Y(:,2);
return