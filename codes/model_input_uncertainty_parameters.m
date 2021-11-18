function [IC_perturbations,params_perturbations,forcing_perturbations]=...
    model_input_uncertainty_parameters(uncertain_states,uncertain_params,...
    uncertain_forcings,restart_info,forcing_vector_length,seed1,seed2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function serves the role of introducing input uncertainties and 
% needs to be modified accordingly. 
% Based on the specified states, parameters, and forcings to be perturbed,
% the uncertainty statistics are specified here and the gaussian parameters
% are used to generate perturbations which are sent back and used to
% generate uncertain model inputs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grab size of domain for spatial field perturbations
[ny,nx]=size(restart_info.params.static_maps.mask);

%% Initial conditions
% diagnose number of uncertain states
nvar=length(uncertain_states);

% Specify IC uncertainty statistics; here using a multiplicative 
% (lognormal) model with an exponential correlation function for spatial 
% correlation.
% range for exponential spatial correlation function (in model pixels)
scale_factor_IC=5; 
% Specify the desired means and variances of LN variables
% unity mean
desired_mean = ones([1 nvar])'; 
% desired standard deviations (or coeffs. of variation if mean=1)
desired_stddev=[0.1 0.1]'; 
% Specify correlation matrix for log-normal variables
desired_corr_matrix = [1 0. ; 0. 1];
% Compute necessary Gaussian parameters
[mu_IC,covariance_matrix_IC]=compute_gaussian_parameters(desired_mean,...
    desired_stddev,desired_corr_matrix);

% Generate IC perturbations
% Assuming correlated spatial fields generate perturbations
[IC_perturbations]=tbcosim_uncond_func(nx,ny,2,scale_factor_IC,1,...
    covariance_matrix_IC,seed1);
% Need to add mean (since tbcosim generates zero-mean variables):
for ii=1:nvar
    IC_perturbations(:,:,ii)=IC_perturbations(:,:,ii)+mu_IC(ii);
end
% Transform variables to lognormal variates:
IC_perturbations=exp(IC_perturbations);

%% Time-invariant parameters
% Specify parameter uncertainty statistics (based on order shown above); 
% here using a multiplicative (lognormal) spatial fields with an 
% exponential correlation function for spatial correlation.

% range for exponential spatial correlation function (in model pixels)
scale_factor_params=5; 
% diagnose number of parameters
nvar=length(uncertain_params);

% Specify the desired means and variances of LN variables
% unity mean
desired_mean = ones([1 nvar])'; 
% desired standard deviations (or coeffs. of variation if mean=1) 
desired_stddev=[0.5 0.5]'; 
% Specify correlation matrix for log-normal variables
desired_corr_matrix = [1 1; 1 1];
% Compute necessary Gaussian parameters
[mu_params,covariance_matrix_params]=compute_gaussian_parameters(...
    desired_mean,desired_stddev,desired_corr_matrix);

% Generate parameter perturbations
% Assuming correlated spatial fields
[params_perturbations]=tbcosim_uncond_func(nx,ny,2,scale_factor_params,...
    1,covariance_matrix_params,seed2);
% Need to add mean (since tbcosim generates zero-mean variables):
for ii=1:size(params_perturbations,3)
    params_perturbations(:,:,ii)=params_perturbations(:,:,ii)+...
        mu_params(ii);
end
% Transform variables to lognormal variates:
params_perturbations=exp(params_perturbations);    
    
%% Time-varying forcings
% Specify parameter uncertainty statistics (based on order shown above); 
% here using a multiplicative (lognormal) model with time correlation only 
% (no spatial fields).
nvar=length(uncertain_forcings);
% Lag-1 Autocorrelation; Note: The "rho" (lag-1 autocorrelation) depends on
% the time step. Here specify the decorrelation time scale (where
% correlation falls off to 1/e) and then calculate the appropriate rho.
decorrelation_time_in_days=4;
% "rho" value:
autocorr_forcings=(1/exp(1))^(1/(decorrelation_time_in_days*24/...
    restart_info.control_params.dt));
% Specify the desired means and variances of LN variables
% unity mean
desired_mean = ones([1 nvar])'; 
% desired standard deviations (or coeffs. of variation if mean=1)
desired_stddev=[0.5]'; 
% Specify correlation matrix for log-normal variables
desired_corr_matrix = [1];
% Compute necessary Gaussian parameters
[mu_forcings,covariance_matrix_forcings]=compute_gaussian_parameters(...
    desired_mean,desired_stddev,desired_corr_matrix);

% Generate forcing perturbations
% if perfectly correlated in time (generate only one perturbation)
if (autocorr_forcings==1) 
    forcing_perturbations=exp(mvnrnd(mu_forcings,...
        covariance_matrix_forcings))';
% generate time correlated perturbations
else 
    forcing_perturbations=NaN(length(mu_forcings),forcing_vector_length);                
    % Generate time-correlated perturbations
    forcing_perturbations(:,1)=exp(mvnrnd(mu_forcings,...
        covariance_matrix_forcings))';
    for i_time=2:forcing_vector_length
        white_noise=mvnrnd(zeros(size(mu_forcings)),...
            covariance_matrix_forcings)';
        forcing_perturbations(:,i_time)=exp(...
            generate_time_correlated_fields(log(...
            forcing_perturbations(:,i_time-1)),autocorr_forcings,...
            white_noise));
    end 
end 

return