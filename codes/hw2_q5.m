% MOD_WET_wrapper.m (script)
% Written by Steve Margulis 
%
% This is a "wrapper" script that allows MOD-WET to be run for different
% inputs. This can be useful for running the model to perform sensitivity 
% tests or in an ensemble context for open loop simulations and/or adapted 
% for use with a data assimilation method.
%
% Note: The MOD-WET model can simply be run by executing the
% "MOD_WET_model_driver" function (without any input arguments). Here we 
% are just modifying inputs through the optional "restart_info" input
% variables.
%
% For simplicity all outputs and perturbed inputs are saved to file for
% later use.
%
% This code assumes that the nominal parameters, initial conditions, and
% meteorological forcing inputs are those specified in or linked to in the
% input files in the MOD-WET path pointed to below. For the specified path
% (i.e. "toolbox_path" below), these would include:
%
% initial conditions specified in:
% "toolbox_path"/chapter11/initialize_model.m
%
% parameters specified in:
% "toolbox_path"/chapter11/set_static_physical_parameters.m
%
% nominal meteorological input file specified in:
% "toolbox_path"/set_control_parameters.m
% 
% Based on these nominal inputs, they are perturbed as shown below before
% being input to a realization of MOD-WET. So for an ensemble of inputs,
% there will be an ensemble of outputs.
clear all;close all

%% Specify path to root directory for the MOD-WET functions
toolbox_path='C:/Users/BANZH/Downloads/2019spring/CEE steve/MOD_WET_2017a/MOD_WET_2017a/MOD_WET_src';
% Add path
addpath(genpath(toolbox_path))

% Root directory where to store outputs
output_root_dir='C:/Users/BANZH/Downloads/2019spring/CEE steve/output/ensemble/';

% Number of ensemble members to run (set to 1 if running only one
% simulation)
N_ensemble_members=200;

%% Call initialization codes to get nominal parameters for model simulation. 
% This first part is just a simple mechanism to load the nominal 
% parameters/initial from the existing MOD-WET parameter/initial condition 
% functions. The desired variable values should be set in those functions 
% and are loaded here. In other words, the nominal inputs should already be
% set in the normal MOD-WET inputs files. By running those functions you 
% will be storing all of the parameters into the "restart_info" structure
% array. The nominal values can then be perturbed and replaced as inputs 
% for a simulation (or ensemble of simulations) via the "restart_info" 
% variable that can be provided as an input to MOD-WET. Beyond
% parameter/forcing perturbations, this can also be used to set the
% start/end time of the simulation, etc. These will be useful in being able
% to use MOD-WET in data assimilation approaches.

% load nominal control parameters (i.e. start day, time step, number of 
% days, etc.)
[restart_info.control_params]=set_control_parameters;

% load nominal static physical parameters (soil parameters, TOPMODEL 
% parameters, etc.
[restart_info.params]=set_static_physical_parameters(restart_info.control_params);
% pre-allocation of variables
[restart_info.states,~,~,~]= ...
                pre_allocate_variables(restart_info.control_params,restart_info.params);
% load nominal initial conditions for state variables (i.e. "old_vars")
[~,~,restart_info.old_vars]= ...
        initialize_model(restart_info.control_params,restart_info.params,restart_info.states);
% set some time series variable and map matrix
tsSD=nan(200,721);
tsET=nan(200,720);
tsQ=nan(200,720);
mapSD=nan(26,19,30,200);
mapET=nan(26,19,30,200);
mapQ=nan(26,19,30,200);
%% Loop through an ensemble of simulations (with different inputs)        
for i_ensemble=1:N_ensemble_members

    %% Here the control parameters, physical parameters, initial conditions, and forcing can be perturbed:

    %% An example of changing the output file name (in order to store different files for different ensemble members
    % Replace output filename with one specific to this ensemble member
    restart_info.control_params.output_filename=[output_root_dir 'ensemble_member_' num2str(i_ensemble) '_output.mat'];
    
    %% An example of running a shorter simulation:
    % This might be useful for calling the model for different time
    % periods (i.e. when using a filter). These are just copied from the
    % set_control_parameters.m.
    % Replace the starting day of the simulation (0 is the beg. of the WY)
    restart_info.control_params.start_day=213; % starting day 
    % Replace the length of the simulation and related input parameters.
    restart_info.control_params.n_days=30; % number of days to run
    restart_info.control_params.ntime=1./restart_info.control_params.dt*24.*restart_info.control_params.n_days;
    % Length of time series to store (needs to be an integer)
    restart_info.control_params.nt=restart_info.control_params.ntime/restart_info.control_params.timeseries_frq2store; 
    % Set index to initialize forcings 
    restart_info.control_params.start_time=1+...     
        restart_info.control_params.start_day*24./restart_info.control_params.dt;

    %% An example of perturbing physical parameters:
    % In the example below, for simplicity it is assumed that the 
    % parameters are uncorrelated (and therefore can be generated 
    % independently) and perturbations are only applied uniformly to the 
    % spatial fields. In the more general case, a spatial field 
    % perturbation could be generated instead, correlations between 
    % parameters could be considered, etc.
    %
    % Perturb parameters with a multiplicative factor that is lognormally 
    % distributed. This can be done using mvnrnd with a specification of 
    % the mean and coefficient of variation (ratio of std. dev. to mean).
    % For spatial field perturbations, the TBCOSIM code should be used.
    nx=19; ny=26;
    % Covariance model type:
    model_type=2; % exponential
    % scale parameter
    scale_factor=3; % range for exponential dist.
    % b parameter:
    b=1;
    nvar=2;
    % Perturb the soil transmissivity:
    % Specify the desired means and variances of LN variables
    desired_mean = [1 1]; % desired mean
    desired_stddev=[1.0 0.15]; % desired standard deviations (or coeffs. of variation if mean=1) 
    % Compute log-normal parameters
    desired_variances = desired_stddev.^2; % variance  
    mu = log((desired_mean.^2) ./ sqrt(desired_variances + desired_mean.^2)); % parameter of the associated normal dist.
    sigma = sqrt(log(desired_variances ./ (desired_mean.^2)+1)); %  parameter of the associated normal dist.

    % Specify correlation matrix for log-normal variables
    corrmat_LN = [1 0;0 1];

    % Transform to covariance matrix for normal dist.
    for i=1:nvar
        for j=1:nvar
            covariance_matrix(i,j)=log( corrmat_LN(i,j) * sqrt( exp(sigma(i)^2)-1)*sqrt(exp(sigma(j)^2)-1) + 1 );        
        end
    end
    % Set seed (here it is set to a random number):
    seed=round(rand*1e7);
    % Call function:
    [SIM_OUTPUTS]=tbcosim_uncond_func(nx,ny,model_type,scale_factor,b,covariance_matrix,seed);
    % Need to add mean (since tbcosim generates zero-mean variables):
    for ii=1:nvar
        SIM_OUTPUTS(:,:,ii)=SIM_OUTPUTS(:,:,ii)+mu(ii);
    end
    % Transform variables to lognormal variates:
    multiplicative_factor=exp(SIM_OUTPUTS);
    % Perturb T0 field
    restart_info.params.T0=multiplicative_factor(:,:,1).*restart_info.params.T0;
    % Need to then re-set dependent parameters/maps
    restart_info.params.static_maps.T0=restart_info.params.T0.*restart_info.params.static_maps.maskNaN; % [-]
    % Perturb surface soil albedo field
    restart_info.params.albedo=multiplicative_factor(:,:,2).*restart_info.params.albedo;
    % Need to then re-set dependent parameters/maps
    restart_info.params.static_maps.albedo=restart_info.params.albedo.*restart_info.params.static_maps.maskNaN; % [-]
    
    %% An example of perturbing initial conditions:
    % Perturb rootzone soil moisture
    % Specify the desired means and variances of LN variables
    nvar=3;
    desired_mean = [1 1 1]; % desired mean
    desired_stddev=[0.2 0.1 0.005]; % desired standard deviations (or coeffs. of variation if mean=1) 
    % Compute log-normal parameters
    desired_variances = desired_stddev.^2; % variance  
    mu = log((desired_mean.^2) ./ sqrt(desired_variances + desired_mean.^2)); % parameter of the associated normal dist.
    sigma = sqrt(log(desired_variances ./ (desired_mean.^2)+1)); %  parameter of the associated normal dist.

    % Specify correlation matrix for log-normal variables
    corrmat_LN = [1 0 0;0 1 0;0 0 1];

    % Transform to covariance matrix for normal dist.
    for i=1:nvar
        for j=1:nvar
            covariance_matrix(i,j)=log( corrmat_LN(i,j) * sqrt( exp(sigma(i)^2)-1)*sqrt(exp(sigma(j)^2)-1) + 1 );        
        end
    end
    % Set seed (here it is set to a random number):
    seed=round(rand*1e7);
    % Call function:
    [SIM_OUTPUTS]=tbcosim_uncond_func(nx,ny,model_type,scale_factor,b,covariance_matrix,seed);
    % Need to add mean (since tbcosim generates zero-mean variables):
    for ii=1:nvar
        SIM_OUTPUTS(:,:,ii)=SIM_OUTPUTS(:,:,ii)+mu(ii);
    end
    % Transform variables to lognormal variates:
    multiplicative_factor=exp(SIM_OUTPUTS);
    % Perturb SD0 field
    restart_info.old_vars.SD0=multiplicative_factor(:,:,1).*restart_info.old_vars.SD0;
    % Pertub Srz0 field
    restart_info.old_vars.Srz0=multiplicative_factor(:,:,2).*restart_info.old_vars.Srz0;
    % Make sure variable is inbounds
    I=find(restart_info.old_vars.Srz0<restart_info.params.static_maps.Srzmin);
    restart_info.old_vars.Srz0(I)=restart_info.params.static_maps.Srzmin(I);
    I=find(restart_info.old_vars.Srz0>restart_info.params.static_maps.Srzmax);
    restart_info.old_vars.Srz0(I)=restart_info.params.static_maps.Srzmax(I);
    % Perturb Td0
    restart_info.old_vars.SD0=multiplicative_factor(:,:,3).*restart_info.old_vars.SD0;
    
    %% An example of perturbing meteorological forcing
    eval(['load ' '''' restart_info.control_params.met_data_filename ''''])
    
    % Perturb precip, Rs and Ta; this is just a simple temporally invariant
    % perturbation; could be temporally/spatially variable
    desired_mean = [1 1 1]; % unity mean
    desired_stddev=[0.5 0.1 0.005]; % desired standard deviations (or coeffs. of variation if mean=1) 
    % Compute log-normal parameters
    desired_variances = desired_stddev.^2; % variance  
    mu = log((desired_mean.^2) ./ sqrt(desired_variances + desired_mean.^2)); % parameter of the associated normal dist.
    sigma = sqrt(log(desired_variances ./ (desired_mean.^2)+1)); %  parameter of the associated normal dist.
    % Specify correlation matrix for log-normal variables
    corrmat_LN = [1 -0.1 -0.1;-0.1 1 0.3;-0.1 0.3 1];
    % Transform to covariance matrix for normal dist.
    for i=1:length(desired_mean)
        for j=1:length(desired_mean)
            covariance_matrix(i,j)=log( corrmat_LN(i,j) * sqrt( exp(sigma(i)^2)-1)*sqrt(exp(sigma(j)^2)-1) + 1 );        
        end
    end
    multiplicative_factor=exp(mvnrnd(mu,covariance_matrix));
    % Perturb precip. field
    PPT=multiplicative_factor(1)*PPT;
    % Perturb Rs field
    SW=multiplicative_factor(2)*SW;
    % Perturb Ta field
    Ta=multiplicative_factor(3)*Ta;
    restart_info.control_params.met_data_filename=[output_root_dir 'ensemble_member_' num2str(i_ensemble) '_met_forcing_input.mat'];
    eval(['save ' '''' restart_info.control_params.met_data_filename '''' ' time gage_elev PPT SW Ta qa U Psfc'])

    %% Call MOD-WET with perturbed inputs
    MOD_WET_model_driver(restart_info)

    %% Load outputs and plot some states/fluxes
    eval(['load ' '''' restart_info.control_params.output_filename ''''])
    tsSD(i_ensemble,:)=states.time_series.SD;
    tsSrz(i_ensemble,:)=states.time_series.Srz;
    tsTsurf(i_ensemble,:)=states.time_series.Tsurf;
    tsSWE(i_ensemble,:)=states.time_series.SWE;
    mapSWE(:,:,:,i_ensemble)=states.maps.SWE;
    tsLE(i_ensemble,:)=fluxes.time_series.LE;
    mapLE(:,:,:,i_ensemble)=fluxes.maps.LE;
    tsRn(i_ensemble,:)=fluxes.time_series.Rn;
    tsQ(i_ensemble,:)=fluxes.time_series.outlet_hydrograph;
end % end ensemble loop

DOY=1:1:201
figure(1);hold on
meanSD=mean(tsSD,1)
plot(DOY,mean(reshape(meanSD(1,2:4825),24,201)),'k','LineWidth',3)
for i=1:20
    figure(1);hold on
    plot(DOY,mean(reshape(tsSD(i,2:4825),24,201)),'c')
end
        xlabel('DOWY');ylabel('SD (m)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','south')
        pause

figure(2);hold on
meanSrz=mean(tsSrz,1)
plot(DOY,mean(reshape(meanSrz(1,2:4825),24,201)),'k','LineWidth',3)
for i=1:20
    figure(2);hold on
    plot(DOY,mean(reshape(tsSrz(i,2:4825),24,201)),'c')
end
        xlabel('DOWY');ylabel('S_{rz} (m)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','southeast')
        pause

figure(3);hold on
I=find(tsTsurf<-273);
tsTsurf(I)=273;
meanTsurf=mean(tsTsurf,1)
plot(DOY,mean(reshape(meanTsurf(1,2:4825),24,201)),'k','LineWidth',3)
for i=1:20
    figure(3);hold on
    plot(DOY,mean(reshape(tsTsurf(i,2:4825),24,201)),'c')
end
        xlabel('DOWY');ylabel('T_{surf} (K)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','southeast')
        pause

figure(4);hold on
meanSWE=mean(tsSWE,1)
plot(DOY,mean(reshape(meanSWE(1,2:4825),24,201)),'k','LineWidth',3)
for i=1:20
    figure(4);hold on
    plot(DOY,mean(reshape(tsSWE(i,2:4825),24,201)),'c')
end
        xlabel('DOWY');ylabel('SWE (m)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','northwest')
        pause
        
figure(5);hold on
meanLE=mean(tsLE,1)
plot(DOY,mean(reshape(meanLE(1,1:4824),24,201)),'k','LineWidth',3)
for i=1:20
    figure(5);hold on
    plot(DOY,mean(reshape(tsLE(i,1:4824),24,201)),'c')
end
        xlabel('DOWY');ylabel('LE (W/m^2)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','northeast')
        pause

figure(6);hold on
meanRn=mean(tsRn,1)
plot(DOY,mean(reshape(meanRn(1,1:4824),24,201)),'k','LineWidth',3)
for i=1:20
    figure(6);hold on
    plot(DOY,mean(reshape(tsRn(i,1:4824),24,201)),'c')
end
        xlabel('DOWY');ylabel('Rn (W/m^2)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','northeast')
        pause

figure(7);hold on
meanQ=mean(tsQ,1)
plot(DOY,mean(reshape(meanQ(1,1:4824),24,201)),'k','LineWidth',3)
for i=1:20
    figure(7);hold on
    plot(DOY,mean(reshape(tsQ(i,1:4824),24,201)),'c')
end
        xlabel('DOWY');ylabel('Q_{outlet} (m^3/s)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','northwest')
        pause

figure(8);hold on
stdSWE=std(mapSWE,0,4)
imagesc(stdSWE(:,:,200));axis image;colorbar;
    caxis([min(min(stdSWE(:,:,200))),max(max(stdSWE(:,:,200)))])
    title(['ensemble std of SWE at DOWY=200'])
    set(gca,'FontSize',18)  

figure(9);hold on
stdLE=std(mapLE,0,4)
imagesc(stdLE(:,:,200));axis image;colorbar;
    caxis([min(min(stdLE(:,:,200))),max(max(stdLE(:,:,200)))])
    title(['ensemble std of LE at DOWY=200'])
    set(gca,'FontSize',18)  
    
figure(10);hold on
meanSWE=mean(mapSWE,4)
imagesc(meanSWE(:,:,200));axis image;colorbar;
    caxis([min(min(meanSWE(:,:,200))),max(max(meanSWE(:,:,200)))])
    title(['ensemble mean of SWE at DOWY=200'])
    set(gca,'FontSize',18)  
    
figure(11);hold on
meanLE=mean(mapLE,4)
imagesc(meanLE(:,:,200));axis image;colorbar;
    caxis([min(min(meanLE(:,:,200))),max(max(meanLE(:,:,200)))])
    title(['ensemble mean of LE at DOWY=200'])
    set(gca,'FontSize',18)  