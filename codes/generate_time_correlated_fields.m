function [field_time2]=generate_time_correlated_fields(field_time1,...
    autocorr,white_noise_field)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates time-varying fields that have a specified 
% de-correlation time scale.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS:
%
%   field_time1: a zero-mean random field at a particular time step
%
%   autcorr: auto-correlation (value between 0 and 1); for a value of 0, 
%   the new field at the next time step will be white noise; for a value 
%   of 1, the field at the next time step will be the same as the previous.
%
%   white_noise_field: a zero mean, unit-variance white noise field.
%
% OUTPUTS:
%
% field_time2: a zero-mean random field at the following time step
%

field_time2=autocorr*field_time1 + sqrt(1-autocorr^2)*white_noise_field;

return
