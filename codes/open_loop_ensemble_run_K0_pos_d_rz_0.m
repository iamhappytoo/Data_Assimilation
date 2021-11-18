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
output_root_dir='C:/Users/BANZH/Downloads/2019spring/CEE Steve/Project_outputs/open_loop_ensemble/';

% Number of ensemble members to run (set to 1 if running only one
% simulation)
N_ensemble_members=20;

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
%[~,~,restart_info.old_vars]= ...
%        initialize_model(restart_info.control_params,restart_info.params,restart_info.states);
% Do load nominal initial conditions in the loop, because the multiplier is
% not 1 mean, so it will have cumulative effects.
% set some time series variable and map matrix
tsSD=nan(20,745);
tsET=nan(20,744);
tsQ=nan(20,744);
tsSrz=nan(20,745);
mapSD=nan(26,19,31,20);
mapET=nan(26,19,31,20);
mapSrz=nan(26,19,31,20);
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
    restart_info.control_params.n_days=31; % number of days to run
    restart_info.control_params.ntime=1./restart_info.control_params.dt*24.*restart_info.control_params.n_days;
    % Length of time series to store (needs to be an integer)
    restart_info.control_params.nt=restart_info.control_params.ntime/restart_info.control_params.timeseries_frq2store; 
    % Set index to initialize forcings 
    restart_info.control_params.start_time=1+...     
        restart_info.control_params.start_day*24./restart_info.control_params.dt;

     %% An example of perturbing physical parameters K0 and d_rz:
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
    % Specify the desired means and SIGMA of multiplicative perturbations
    desired_mean = [0.2 0.1]; % desired mean
    desired_stddev= [0.1 0.01]; % desired standard deviations (or coeffs. of variation if mean=1) 
    % Generate random multipliers based on the given mean and stddev
    additive_factor=mvnrnd(desired_mean,desired_stddev);
    % Reload the restart info of params
    [restart_info.params]=set_static_physical_parameters(restart_info.control_params);
    % Perturb K0, d_rz 
    restart_info.params.static_maps.K0=restart_info.params.static_maps.K0*additive_factor(1)
    restart_info.params.d_rz=restart_info.params.d_rz+additive_factor(2)
    %restart_info.params.static_maps.K0=restart_info.params.static_maps.K0+2

    %% An example of perturbing initial conditions:
    % load nominal initial conditions for state variables (i.e. "old_vars")
    [~,~,restart_info.old_vars]= ...
            initialize_model(restart_info.control_params,restart_info.params,restart_info.states);
    % Perturb initial saturation deficit
    % Specify the desired means and SIGMA of multiplicative perturbations
    nvar=1;
    desired_mean = 1; % desired mean
    desired_stddev= 0.15; % desired standard deviations (or coeffs. of variation if mean=1) 
    % Generate random multipliers based on the given mean and stddev
    multiplicative_factor=mvnrnd(desired_mean,desired_stddev);
    % Perturb SD0 field
    restart_info.old_vars.SD0=multiplicative_factor*restart_info.old_vars.SD0;
 %% An example of perturbing meteorological forcing
    eval(['load ' '''' restart_info.control_params.met_data_filename ''''])
    
    % Perturb precip; this is just a simple temporally invariant
    % perturbation; could be temporally/spatially variable
    desired_mean = 1.2; % unity mean
    desired_stddev=0.1; % desired standard deviations (or coeffs. of variation if mean=1) 
    multiplicative_factor=mvnrnd(desired_mean,desired_stddev);
    PPT=multiplicative_factor(1)*PPT;
    restart_info.control_params.met_data_filename=[output_root_dir 'ensemble_member_' num2str(i_ensemble) '_met_forcing_input.mat'];
    eval(['save ' '''' restart_info.control_params.met_data_filename '''' ' time gage_elev PPT SW Ta qa U Psfc'])    %% Call MOD-WET with perturbed inputs
    %% Call MOD-WET with perturbed inputs
    MOD_WET_model_driver(restart_info)

    %% Load outputs and plot some states/fluxes
    eval(['load ' '''' restart_info.control_params.output_filename ''''])
    tsSD(i_ensemble,:)=states.time_series.SD;
    tsET(i_ensemble,:)=fluxes.time_series.ET;
    tsQ(i_ensemble,:)=fluxes.time_series.outlet_hydrograph;
    tsSrz(i_ensemble,:)=states.time_series.Srz;
    mapSD(:,:,:,i_ensemble)=states.maps.SD;
    mapET(:,:,:,i_ensemble)=fluxes.maps.ET;
    mapSrz(:,:,:,i_ensemble)=states.maps.Srz;
end % end ensemble loop
save('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble/tsQ.mat','tsQ')
save('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble/tsET.mat','tsET')
save('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble/tsSD.mat','tsSD')
save('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble/mapET.mat','mapET')
save('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble/mapSD.mat','mapSD')
save('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble/mapSrz.mat','mapSrz')
save('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble/tsSrz.mat','tsSrz')
DOY=213:1:243;
figure(1);hold on
pixSD=squeeze(mean(mapSD(14,17,:,:),4));
plot(DOY,pixSD,'k','LineWidth',3)
for i=1:20
    figure(1);hold on
    plot(DOY,squeeze(mapSD(14,17,:,i)),'c');
end
        xlabel('DOWY');ylabel('SD (m)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','south')
figure(2);hold on
pixSD=squeeze(mean(mapSD(19,3,:,:),4));
plot(DOY,pixSD,'k','LineWidth',3)
for i=1:20
    figure(2);hold on
    plot(DOY,squeeze(mapSD(19,3,:,i)),'c');
end
        xlabel('DOWY');ylabel('SD (m)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','south')
        pause

figure(3);hold on
pixSD=squeeze(mean(mapSD(8,10,:,:),4));
plot(DOY,pixSD,'k','LineWidth',3)
for i=1:20
    figure(3);hold on
    plot(DOY,squeeze(mapSD(8,10,:,i)),'c');
end
        xlabel('DOWY');ylabel('SD (m)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','south')
        pause

figure(4);hold on
pixSrz=squeeze(mean(mapSrz(14,17,:,:),4));
plot(DOY,pixSrz,'k','LineWidth',3)
for i=1:20
    figure(4);hold on
    plot(DOY,squeeze(mapSrz(14,17,:,i)),'c');
end
        xlabel('DOWY');ylabel('Srz (m)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','south')
        pause

figure(5);hold on
pixSrz=squeeze(mean(mapSrz(19,3,:,:),4));
plot(DOY,pixSrz,'k','LineWidth',3)
for i=1:20
    figure(5);hold on
    plot(DOY,squeeze(mapSrz(19,3,:,i)),'c');
end
        xlabel('DOWY');ylabel('Srz (m)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','south')
        pause

figure(6);hold on
pixSrz=squeeze(mean(mapSrz(8,10,:,:),4));
plot(DOY,pixSrz,'k','LineWidth',3)
for i=1:20
    figure(6);hold on
    plot(DOY,squeeze(mapSrz(8,10,:,i)),'c');
end
        xlabel('DOWY');ylabel('Srz (m)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','south')
        pause

figure(7);hold on
pixET=squeeze(mean(mapET(14,17,:,:),4));
plot(DOY,pixET,'k','LineWidth',3)
for i=1:20
    figure(7);hold on
    plot(DOY,squeeze(mapET(14,17,:,i)),'c');
end
        xlabel('DOWY');ylabel('ET (m/h)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','south')
        pause

figure(8);hold on
pixET=squeeze(mean(mapET(19,3,:,:),4));
plot(DOY,pixET,'k','LineWidth',3)
for i=1:20
    figure(8);hold on
    plot(DOY,squeeze(mapET(19,3,:,i)),'c');
end
        xlabel('DOWY');ylabel('ET (m/h)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','south')
        pause

figure(9);hold on
pixET=squeeze(mean(mapET(8,10,:,:),4));
plot(DOY,pixET,'k','LineWidth',3)
for i=1:20
    figure(9);hold on
    plot(DOY,squeeze(mapET(8,10,:,i)),'c');
end
        xlabel('DOWY');ylabel('ET (m/h)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','south')
        pause








































figure(4);hold on
meanET=mean(tsET,1);
plot(DOY,mean(reshape(meanET(1,1:744),24,31)),'k','LineWidth',3)
for i=1:20
    figure(4);hold on
    plot(DOY,mean(reshape(tsET(i,1:744),24,31)),'c')
end
        xlabel('DOWY');ylabel('ET (m/h)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','southeast')
        pause

figure(11);hold on
meanQ=mean(tsQ,1);
plot(DOY,mean(reshape(meanQ(1,1:744),24,31)),'k','LineWidth',3)
for i=1:20
    figure(11);hold on
    plot(DOY,mean(reshape(tsQ(i,1:744),24,31)),'c')
end
        xlabel('DOWY');ylabel('outlet hydrograph (m^3/s)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','southeast')
        pause

figure(4);hold on
stdSD=std(mapSD,0,4);
imagesc(stdSD(:,:,31));axis image;colorbar;
    caxis([min(min(stdSD(:,:,31))),max(max(stdSD(:,:,31)))])
    title(['ensemble std of SD at DOWY=243'])
    set(gca,'FontSize',18)  

figure(5);hold on
stdET=std(mapET,0,4);
imagesc(stdET(:,:,31));axis image;colorbar;
    caxis([min(min(stdET(:,:,31))),max(max(stdET(:,:,31)))])
    title(['ensemble std of LE at DOWY=243'])
    set(gca,'FontSize',18)  