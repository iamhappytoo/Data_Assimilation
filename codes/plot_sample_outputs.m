%% Script to plot sample outputs from MOD-WET DA examples
% This is setup specifically for the default code case where Srz is being 
% assimilated (pixel-wise). Need to generalize code if want to run 
% automatically for other scenarios.

clear all;close all
rng('default');

% Set case (either 'true_gen' or 'EnKF') to generate some sample plots:
%run_case='true_gen'
run_case='EnKF'

% If run_case=EnKF this is the number of ensemble members that were used
N_ensemble_members=6;

% Directory where output files are located
root_dir='C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/true_output_daily_0.1';
%root_dir='/Users/margulis/Box Sync/CEE251D/MOD_WET_outputs/low_elev_basin/DA_outputs/without_localization/';

switch run_case
   
    case 'true_gen'
        % load true states
        load([root_dir 'true_output.mat'])
        % load measurements
        load([root_dir 'measurements.mat'])
        
        % Generate measurement index to match states
        meas_index=meas_day_vector*24*restart_info.control_params.map_frq2store*restart_info.control_params.dt;
        % Store measurements for ease of plotting
        zmeas_tmp=NaN(size(params.static_maps.mask,1),size(params.static_maps.mask,2),length(meas_index));
        for i=1:length(meas_index)
            tmp=NaN(size(params.static_maps.mask));
            tmp(params.static_maps.mask==1)=ZMEAS(:,i);
            zmeas_tmp(:,:,i)=tmp;
        end
        
        %% Plot basin DEM for reference
        figure(1)
        imagesc(params.static_maps.elev.*params.static_maps.maskNaN)
        axis image
        colorbar
        title('DEM of watershed')
        set(gca,'FontSize',14)
        
        %% Plot select states maps and measurements maps
        for i=1:length(meas_index)
            figure
            subplot(2,3,5)
            imagesc(zmeas_tmp(:,:,i)); axis image; colorbar
            title(['Measurement on model time index: ' num2str(meas_index(i))])
            set(gca,'FontSize',14)
            subplot(2,3,1)
            if i==1
                plot_index=mean([0;meas_index(i)]);
            else
                plot_index=mean([meas_index(i-1);meas_index(i)]);
            end
            imagesc(states.maps.Srz(:,:,plot_index)); axis image; colorbar
            title(['True Srz on model time index: ' num2str(plot_index)])
            set(gca,'FontSize',14)
            subplot(2,3,2)
            plot_index=meas_index(i);
            imagesc(states.maps.Srz(:,:,plot_index)); axis image; colorbar
            title(['True Srz on model time index: ' num2str(plot_index)])
            set(gca,'FontSize',14)
            subplot(2,3,3)
            if i==length(meas_index)
                plot_index=mean([meas_index(i);size(states.maps.Srz,3)]);
            else
                plot_index=mean([meas_index(i+1);meas_index(i)]);
            end
            imagesc(states.maps.Srz(:,:,plot_index)); axis image; colorbar
            title(['True Srz on model time index: ' num2str(plot_index)])
            set(gca,'FontSize',14)
        end
        
        %% Plot individual pixels and their measurements
        % Pick three random pixels
        I=find(params.static_maps.mask==1);
        tmp=I(randperm(length(I)));
        pixel_indices=tmp(1:3);
        for ipix=1:length(pixel_indices)
           figure
           [row,col] = ind2sub(size(params.static_maps.mask),pixel_indices(ipix));
           plot(squeeze(states.maps.Srz(row,col,:)));hold on
           plot(meas_index,squeeze(zmeas_tmp(row,col,:)),'o','MarkerSize',10)
           xlabel('Model time step');ylabel('Srz (m)')
           title(['True time series and measurements for pixel: (' num2str(row) ',' num2str(col) ')'])
           set(gca,'FontSize',14);grid
        end
        
    case 'EnKF'
        
        % load true states
        load([root_dir 'true_output.mat'])
        % load measurements
        load([root_dir 'measurements.mat'])
        % Generate measurement index to match states
        meas_index=meas_day_vector*24*restart_info.control_params.map_frq2store*restart_info.control_params.dt;
        % Store measurements for ease of plotting
        zmeas_tmp=NaN(size(params.static_maps.mask,1),size(params.static_maps.mask,2),length(meas_index));
        for i=1:length(meas_index)
            tmp=NaN(size(params.static_maps.mask));
            tmp(params.static_maps.mask==1)=ZMEAS(:,i);
            zmeas_tmp(:,:,i)=tmp;
        end
        
        % Compute ensemble mean/std. dev. fields
        for irep=1:N_ensemble_members
            load([root_dir 'ensemble_member_' num2str(irep) '_merged_output.mat'],'STATE_maps')
            clear tmp
            if irep==1
                % initializes
                Srz_mean=STATE_maps.Srz;
                Srz2=(STATE_maps.Srz).^2;
                Srz_max=STATE_maps.Srz;
                Srz_min=STATE_maps.Srz;
            elseif irep<N_ensemble_members
                % augments sum across ensemble
                Srz_mean=Srz_mean+STATE_maps.Srz;
                Srz2=Srz2+(STATE_maps.Srz).^2;
                tmp(:,:,:,1)=Srz_max;tmp(:,:,:,2)=STATE_maps.Srz;
                Srz_max=squeeze(max(tmp,[],4));
                tmp(:,:,:,1)=Srz_min;tmp(:,:,:,2)=STATE_maps.Srz;
                Srz_min=squeeze(min(tmp,[],4));
            elseif irep==N_ensemble_members
                % Computes mean on last replicate
                Srz_mean=(Srz_mean+STATE_maps.Srz)/N_ensemble_members;
                Srz2=(Srz2+(STATE_maps.Srz).^2)/N_ensemble_members - ...
                    Srz_mean.^2;
                Srz_stddev=sqrt(Srz2);
                tmp(:,:,:,1)=Srz_max;tmp(:,:,:,2)=STATE_maps.Srz;
                Srz_max=squeeze(max(tmp,[],4));
                tmp(:,:,:,1)=Srz_min;tmp(:,:,:,2)=STATE_maps.Srz;
                Srz_min=squeeze(min(tmp,[],4));
            end      
        end
        
        %% Plot some sample time series with full ensemble statistics
        % Pick three random pixels
        I=find(params.static_maps.mask==1);
        tmp=I(randperm(length(I)));
        pixel_indices=tmp(1:3);
        for ipix=1:length(pixel_indices)
           figure
           [row,col] = ind2sub(size(params.static_maps.mask),pixel_indices(ipix));
           plot(squeeze(states.maps.Srz(row,col,:)),'Color',[0 0.45 0.74],'LineWidth',3);hold on
           plot(meas_index,squeeze(zmeas_tmp(row,col,:)),'ko','LineWidth',3,'MarkerSize',15)
           plot(squeeze(Srz_mean(row,col,:)),'Color',[0.85 0.33 0.1],'LineWidth',2);
           plot(squeeze(Srz_mean(row,col,:)+2*Srz_stddev(row,col,:)),'--','Color',[0.85 0.33 0.1],'LineWidth',1);
           plot(squeeze(Srz_mean(row,col,:)-2*Srz_stddev(row,col,:)),'--','Color',[0.85 0.33 0.1],'LineWidth',1);
           %plot(squeeze(Srz_max(row,col,:)),'--','Color',[0.85 0.33 0.1],'LineWidth',1);
           %plot(squeeze(Srz_min(row,col,:)),'--','Color',[0.85 0.33 0.1],'LineWidth',1);
           xlabel('Model time step');ylabel('Srz (m)')
           title(['Time series and measurements for pixel: (' num2str(row) ',' num2str(col) ')'])
           set(gca,'FontSize',14);grid
           legend('True','Meas.','EnKF Ens. Mean','+/- 2 Std. Dev.')
           %legend('True','Meas.','EnKF Ens. Mean','EnKF Ens. Max.','EnKF Ens. Min.')
        end
        
        
        %% Plot ensemble mean maps
        for i=1:length(meas_index)
            figure
            subplot(3,3,5)
            imagesc(zmeas_tmp(:,:,i)); axis image; colorbar
            title(['Measurement on model time index: ' num2str(meas_index(i))])
            set(gca,'FontSize',14)
            subplot(3,3,1)
            if i==1
                plot_index=mean([0;meas_index(i)]);
            else
                plot_index=mean([meas_index(i-1);meas_index(i)]);
            end
            imagesc(states.maps.Srz(:,:,plot_index)); axis image; colorbar
            title(['True Srz on model time index: ' num2str(plot_index)])
            set(gca,'FontSize',14)
            subplot(3,3,7)
            imagesc(Srz_mean(:,:,plot_index)); axis image; colorbar
            title(['EnKF Ens. Mean Srz on model time index: ' num2str(plot_index)])
            set(gca,'FontSize',14)
            
            subplot(3,3,2)
            plot_index=meas_index(i);
            imagesc(states.maps.Srz(:,:,plot_index)); axis image; colorbar
            title(['True Srz on model time index: ' num2str(plot_index)])
            set(gca,'FontSize',14)
            subplot(3,3,8)
            imagesc(Srz_mean(:,:,plot_index)); axis image; colorbar
            title(['EnKF Ens. Srz on model time index: ' num2str(plot_index)])
            set(gca,'FontSize',14)
            
            subplot(3,3,3)
            if i==length(meas_index)
                plot_index=mean([meas_index(i);size(states.maps.Srz,3)]);
            else
                plot_index=mean([meas_index(i+1);meas_index(i)]);
            end
            imagesc(states.maps.Srz(:,:,plot_index)); axis image; colorbar
            title(['Srz on model time index: ' num2str(plot_index)])
            set(gca,'FontSize',14)
            subplot(3,3,9)
            imagesc(Srz_mean(:,:,plot_index)); axis image; colorbar
            title(['EnKF Ens. Srz on model time index: ' num2str(plot_index)])
            set(gca,'FontSize',14)
        end
end