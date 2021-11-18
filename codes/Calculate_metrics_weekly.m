%% Load output files for calculating metrics
%% Load open loop output
openmapET=load('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble1year/mapET.mat');
openmapSD=load('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble1year/mapSD.mat');
openmapSrz=load('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble1year/mapSrz.mat');
opentsET=load('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble1year/tsET.mat');
opentsSD=load('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble1year/tsSD.mat');
opentsSrz=load('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble1year/tsSrz.mat');
opentsQ=load('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble1year/tsQ.mat');
%% weekly
% weekly 0.002 Enkf output
enkf1=load('C:/Users/BANZH/Downloads/2019spring/CEE Steve/Project_outputs/long_enkf_output_weekly_0.002/merged_enkf_outputs.mat');
truth1=load('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/long_true_output_weekly_0.002/true_output.mat');
% weekly 0.02 Enkf output
enkf2=load('C:/Users/BANZH/Downloads/2019spring/CEE Steve/Project_outputs/long_enkf_output_weekly_0.02/merged_enkf_outputs.mat');
truth2=load('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/long_true_output_weekly_0.02/true_output.mat');
% weekly 0.05 Enkf output
enkf3=load('C:/Users/BANZH/Downloads/2019spring/CEE Steve/Project_outputs/long_enkf_output_weekly_0.05/merged_EnKF_outputs.mat');
truth3=load('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/long_true_output_weekly_0.05/true_output.mat');
% weekly 0.1 Enkf output
enkf4=load('C:/Users/BANZH/Downloads/2019spring/CEE Steve/Project_outputs/long_enkf_output_weekly_0.1/merged_enkf_outputs.mat');
truth4=load('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/long_true_output_weekly_0.1/true_output.mat');
%% Calculate RMSE prior using truthgen when V=0.002
sumseSD=zeros(1,5832);
sumseSrz=zeros(1,5832);
sumseET=zeros(1,5832);
sumseQ=zeros(1,5832);
for i=1:20
    sumseSD=sumseSD+(truth1.states.time_series.SD(2:5833)-opentsSD.tsSD(i,2:5833)).^2;
    sumseSrz=sumseSrz+(truth1.states.time_series.Srz(2:5833)-opentsSrz.tsSrz(i,2:5833)).^2;
    sumseET=sumseET+(truth1.fluxes.time_series.ET-opentsET.tsET(i,:)).^2;
    sumseQ=sumseQ+(truth1.fluxes.time_series.outlet_hydrograph-opentsQ.tsQ(i,:)).^2;
end
rmseSD=sqrt(sumseSD/20);
rmseSrz=sqrt(sumseSrz/20);
rmseET=sqrt(sumseET/20);
rmseQ=sqrt(sumseQ/20);
DOY=1:1:5832;
%% Calculate RMSE posterior for measV=0.002
sumseSD1=zeros(1,5832);
sumseSrz1=zeros(1,5832);
sumseET1=zeros(1,5832);
sumseQ1=zeros(1,5832);
pixel_indices=find(restart_info.params.static_maps.mask==1);
psd=permute(enkf1.STATE_maps.SD,[3,4,1,2]);
psrz=permute(enkf1.STATE_maps.Srz,[3,4,1,2]);
pet=permute(enkf1.FLUX_maps.ET,[3,4,1,2]);
for i=1:20
    mpsd=mean(psd(:,i,pixel_indices),3);
    mpsrz=mean(psrz(:,i,pixel_indices),3);
    mpet=mean(pet(:,i,pixel_indices),3);
    mpq=enkf1.FLUX_time_series.outlet_hydrograph(:,i);
    sumseSD1=sumseSD1+(truth1.states.time_series.SD(2:5833)-permute(mpsd,[2,1])).^2;
    sumseSrz1=sumseSrz1+(truth1.states.time_series.Srz(2:5833)-permute(mpsrz,[2,1])).^2;
    sumseET1=sumseET1+(truth1.fluxes.time_series.ET-permute(mpet,[2,1])).^2;
    sumseQ1=sumseQ1+(truth1.fluxes.time_series.outlet_hydrograph-permute(mpq,[2,1])).^2;
end
rmseSD1=sqrt(sumseSD1/20);
rmseSrz1=sqrt(sumseSrz1/20);
rmseET1=sqrt(sumseET1/20);
rmseQ1=sqrt(sumseQ1/20);

%% RMSE plot when V=0.002
figure(1);
subplot(4,1,1);
plot(DOY,squeeze(rmseSD),'--r','LineWidth',1);hold on
plot(DOY,squeeze(rmseSD1),'--b','LineWidth',1);
ylabel('RMSE SD');
set(gca,'FontSize',8);grid
legend('prior','post.');
subplot(4,1,2);
plot(DOY,squeeze(rmseSrz),'--r','LineWidth',1);hold on
plot(DOY,squeeze(rmseSrz1),'--b','LineWidth',1);
ylabel('RMSE Srz');
set(gca,'FontSize',8);grid
legend('prior','post.');
subplot(4,1,3);
plot(DOY,squeeze(rmseET),'--r','LineWidth',1);hold on
plot(DOY,squeeze(rmseET1),'--b','LineWidth',1);
ylabel('RMSE ET');
set(gca,'FontSize',8);grid
legend('prior','post.');
subplot(4,1,4);
plot(DOY,squeeze(rmseQ),'--r','LineWidth',1);hold on
plot(DOY,squeeze(rmseQ1),'--b','LineWidth',1);
xlabel('Model time step');ylabel('prior.RMSE Q');
set(gca,'FontSize',8);grid
legend('prior','post.');
mkdir('C:\Users\BANZH\Downloads\2019spring\CEE Steve\Project_outputs\Calculate_metrics_long\Figures_weekly');
cd('C:\Users\BANZH\Downloads\2019spring\CEE Steve\Project_outputs\Calculate_metrics_long\Figures_weekly');
saveas(figure(1),[pwd '/RMSE0.002.png']);        
close(figure(1));

%
%% Calculate RMSE prior using truthgen when V=0.02
sumseSD=zeros(1,5832);
sumseSrz=zeros(1,5832);
sumseET=zeros(1,5832);
sumseQ=zeros(1,5832);
for i=1:20
    sumseSD=sumseSD+(truth2.states.time_series.SD(2:5833)-opentsSD.tsSD(i,2:5833)).^2;
    sumseSrz=sumseSrz+(truth2.states.time_series.Srz(2:5833)-opentsSrz.tsSrz(i,2:5833)).^2;
    sumseET=sumseET+(truth2.fluxes.time_series.ET-opentsET.tsET(i,:)).^2;
    sumseQ=sumseQ+(truth2.fluxes.time_series.outlet_hydrograph-opentsQ.tsQ(i,:)).^2;
end
rmseSD=sqrt(sumseSD/20);
rmseSrz=sqrt(sumseSrz/20);
rmseET=sqrt(sumseET/20);
rmseQ=sqrt(sumseQ/20);
%% Calculate RMSE posterior for measV=0.02
sumseSD1=zeros(1,5832);
sumseSrz1=zeros(1,5832);
sumseET1=zeros(1,5832);
sumseQ1=zeros(1,5832);
pixel_indices=find(restart_info.params.static_maps.mask==1);
psd=permute(enkf2.STATE_maps.SD,[3,4,1,2]);
psrz=permute(enkf2.STATE_maps.Srz,[3,4,1,2]);
pet=permute(enkf2.FLUX_maps.ET,[3,4,1,2]);
for i=1:20
    mpsd=mean(psd(:,i,pixel_indices),3);
    mpsrz=mean(psrz(:,i,pixel_indices),3);
    mpet=mean(pet(:,i,pixel_indices),3);
    mpq=enkf2.FLUX_time_series.outlet_hydrograph(:,i);
    sumseSD1=sumseSD1+(truth2.states.time_series.SD(2:5833)-permute(mpsd,[2,1])).^2;
    sumseSrz1=sumseSrz1+(truth2.states.time_series.Srz(2:5833)-permute(mpsrz,[2,1])).^2;
    sumseET1=sumseET1+(truth2.fluxes.time_series.ET-permute(mpet,[2,1])).^2;
    sumseQ1=sumseQ1+(truth2.fluxes.time_series.outlet_hydrograph-permute(mpq,[2,1])).^2;
end
rmseSD1=sqrt(sumseSD1/20);
rmseSrz1=sqrt(sumseSrz1/20);
rmseET1=sqrt(sumseET1/20);
rmseQ1=sqrt(sumseQ1/20);
%% RMSE plot when V=0.02, since the truth is also changed, so we cannot plot the same prior, so figure2.
figure(2);
subplot(4,1,1);
plot(DOY,squeeze(rmseSD),'--r','LineWidth',1);hold on
plot(DOY,squeeze(rmseSD1),'--b','LineWidth',1);
ylabel('RMSE SD');
set(gca,'FontSize',8);grid
legend('prior','post.');
subplot(4,1,2);
plot(DOY,squeeze(rmseSrz),'--r','LineWidth',1);hold on
plot(DOY,squeeze(rmseSrz1),'--b','LineWidth',1);
ylabel('RMSE Srz');
set(gca,'FontSize',8);grid
legend('prior','post.');
subplot(4,1,3);
plot(DOY,squeeze(rmseET),'--r','LineWidth',1);hold on
plot(DOY,squeeze(rmseET1),'--b','LineWidth',1);
ylabel('RMSE ET');
set(gca,'FontSize',8);grid
legend('prior','post.');
subplot(4,1,4);
plot(DOY,squeeze(rmseQ),'--r','LineWidth',1);hold on
plot(DOY,squeeze(rmseQ1),'--b','LineWidth',1);
xlabel('Model time step');ylabel('prior.RMSE Q');
set(gca,'FontSize',8);grid
legend('prior','post.');
cd('C:\Users\BANZH\Downloads\2019spring\CEE Steve\Project_outputs\Calculate_metrics_long\Figures_weekly');
saveas(figure(2),[pwd '/RMSE0.02.png']);        
close(figure(2));

%
%% Calculate RMSE prior using truthgen when V=0.05
sumseSD=zeros(1,5832);
sumseSrz=zeros(1,5832);
sumseET=zeros(1,5832);
sumseQ=zeros(1,5832);
for i=1:20
    sumseSD=sumseSD+(truth3.states.time_series.SD(2:5833)-opentsSD.tsSD(i,2:5833)).^2;
    sumseSrz=sumseSrz+(truth3.states.time_series.Srz(2:5833)-opentsSrz.tsSrz(i,2:5833)).^2;
    sumseET=sumseET+(truth3.fluxes.time_series.ET-opentsET.tsET(i,:)).^2;
    sumseQ=sumseQ+(truth3.fluxes.time_series.outlet_hydrograph-opentsQ.tsQ(i,:)).^2;
end
rmseSD=sqrt(sumseSD/20);
rmseSrz=sqrt(sumseSrz/20);
rmseET=sqrt(sumseET/20);
rmseQ=sqrt(sumseQ/20);
%% Calculate RMSE posterior for measV=0.05
sumseSD1=zeros(1,5832);
sumseSrz1=zeros(1,5832);
sumseET1=zeros(1,5832);
sumseQ1=zeros(1,5832);
pixel_indices=find(restart_info.params.static_maps.mask==1);
psd=permute(enkf3.STATE_maps.SD,[3,4,1,2]);
psrz=permute(enkf3.STATE_maps.Srz,[3,4,1,2]);
pet=permute(enkf3.FLUX_maps.ET,[3,4,1,2]);
for i=1:20
    mpsd=mean(psd(:,i,pixel_indices),3);
    mpsrz=mean(psrz(:,i,pixel_indices),3);
    mpet=mean(pet(:,i,pixel_indices),3);
    mpq=enkf3.FLUX_time_series.outlet_hydrograph(:,i);
    sumseSD1=sumseSD1+(truth3.states.time_series.SD(2:5833)-permute(mpsd,[2,1])).^2;
    sumseSrz1=sumseSrz1+(truth3.states.time_series.Srz(2:5833)-permute(mpsrz,[2,1])).^2;
    sumseET1=sumseET1+(truth3.fluxes.time_series.ET-permute(mpet,[2,1])).^2;
    sumseQ1=sumseQ1+(truth3.fluxes.time_series.outlet_hydrograph-permute(mpq,[2,1])).^2;
end
rmseSD1=sqrt(sumseSD1/20);
rmseSrz1=sqrt(sumseSrz1/20);
rmseET1=sqrt(sumseET1/20);
rmseQ1=sqrt(sumseQ1/20);
%% RMSE plot when V=0.05, since the truth is also changed, so we cannot plot the same prior, so figure2.
figure(3);
subplot(4,1,1);
plot(DOY,squeeze(rmseSD),'--r','LineWidth',1);hold on
plot(DOY,squeeze(rmseSD1),'--b','LineWidth',1);
ylabel('RMSE SD');
set(gca,'FontSize',8);grid
legend('prior','post.');
subplot(4,1,2);
plot(DOY,squeeze(rmseSrz),'--r','LineWidth',1);hold on
plot(DOY,squeeze(rmseSrz1),'--b','LineWidth',1);
ylabel('RMSE Srz');
set(gca,'FontSize',8);grid
legend('prior','post.');
subplot(4,1,3);
plot(DOY,squeeze(rmseET),'--r','LineWidth',1);hold on
plot(DOY,squeeze(rmseET1),'--b','LineWidth',1);
ylabel('RMSE ET');
set(gca,'FontSize',8);grid
legend('prior','post.');
subplot(4,1,4);
plot(DOY,squeeze(rmseQ),'--r','LineWidth',1);hold on
plot(DOY,squeeze(rmseQ1),'--b','LineWidth',1);
xlabel('Model time step');ylabel('prior.RMSE Q');
set(gca,'FontSize',8);grid
legend('prior','post.');
cd('C:\Users\BANZH\Downloads\2019spring\CEE Steve\Project_outputs\Calculate_metrics_long\Figures_weekly');
saveas(figure(3),[pwd '/RMSE0.05.png']);        
close(figure(3));

%
%% Calculate RMSE prior using truthgen when V=0.1
sumseSD=zeros(1,5832);
sumseSrz=zeros(1,5832);
sumseET=zeros(1,5832);
sumseQ=zeros(1,5832);
for i=1:20
    sumseSD=sumseSD+(truth4.states.time_series.SD(2:5833)-opentsSD.tsSD(i,2:5833)).^2;
    sumseSrz=sumseSrz+(truth4.states.time_series.Srz(2:5833)-opentsSrz.tsSrz(i,2:5833)).^2;
    sumseET=sumseET+(truth4.fluxes.time_series.ET-opentsET.tsET(i,:)).^2;
    sumseQ=sumseQ+(truth4.fluxes.time_series.outlet_hydrograph-opentsQ.tsQ(i,:)).^2;
end
rmseSD=sqrt(sumseSD/20);
rmseSrz=sqrt(sumseSrz/20);
rmseET=sqrt(sumseET/20);
rmseQ=sqrt(sumseQ/20);
%% Calculate RMSE posterior for measV=0.05
sumseSD1=zeros(1,5832);
sumseSrz1=zeros(1,5832);
sumseET1=zeros(1,5832);
sumseQ1=zeros(1,5832);
pixel_indices=find(restart_info.params.static_maps.mask==1);
psd=permute(enkf4.STATE_maps.SD,[3,4,1,2]);
psrz=permute(enkf4.STATE_maps.Srz,[3,4,1,2]);
pet=permute(enkf4.FLUX_maps.ET,[3,4,1,2]);
for i=1:20
    mpsd=mean(psd(:,i,pixel_indices),3);
    mpsrz=mean(psrz(:,i,pixel_indices),3);
    mpet=mean(pet(:,i,pixel_indices),3);
    mpq=enkf4.FLUX_time_series.outlet_hydrograph(:,i);
    sumseSD1=sumseSD1+(truth4.states.time_series.SD(2:5833)-permute(mpsd,[2,1])).^2;
    sumseSrz1=sumseSrz1+(truth4.states.time_series.Srz(2:5833)-permute(mpsrz,[2,1])).^2;
    sumseET1=sumseET1+(truth4.fluxes.time_series.ET-permute(mpet,[2,1])).^2;
    sumseQ1=sumseQ1+(truth4.fluxes.time_series.outlet_hydrograph-permute(mpq,[2,1])).^2;
end
rmseSD1=sqrt(sumseSD1/20);
rmseSrz1=sqrt(sumseSrz1/20);
rmseET1=sqrt(sumseET1/20);
rmseQ1=sqrt(sumseQ1/20);
%% RMSE plot when V=0.05, since the truth is also changed, so we cannot plot the same prior, so figure2.
figure(4);
subplot(4,1,1);
plot(DOY,squeeze(rmseSD),'--r','LineWidth',1);hold on
plot(DOY,squeeze(rmseSD1),'--b','LineWidth',1);
ylabel('RMSE SD');
set(gca,'FontSize',8);grid
legend('prior','post.');
subplot(4,1,2);
plot(DOY,squeeze(rmseSrz),'--r','LineWidth',1);hold on
plot(DOY,squeeze(rmseSrz1),'--b','LineWidth',1);
ylabel('RMSE Srz');
set(gca,'FontSize',8);grid
legend('prior','post.');
subplot(4,1,3);
plot(DOY,squeeze(rmseET),'--r','LineWidth',1);hold on
plot(DOY,squeeze(rmseET1),'--b','LineWidth',1);
ylabel('RMSE ET');
set(gca,'FontSize',8);grid
legend('prior','post.');
subplot(4,1,4);
plot(DOY,squeeze(rmseQ),'--r','LineWidth',1);hold on
plot(DOY,squeeze(rmseQ1),'--b','LineWidth',1);
xlabel('Model time step');ylabel('prior.RMSE Q');
set(gca,'FontSize',8);grid
legend('prior','post.');
cd('C:\Users\BANZH\Downloads\2019spring\CEE Steve\Project_outputs\Calculate_metrics_long\Figures_weekly');
saveas(figure(4),[pwd '/RMSE0.1.png']);        
close(figure(4));
%% Calculate Variance prior and posterior and plot
meanSD=mean(opentsSD.tsSD,1);
meanSrz=mean(opentsSrz.tsSrz,1);
meanET=mean(opentsET.tsET,1);
meanQ=mean(opentsQ.tsQ,1);
sumvSD=zeros(5,5832);
sumvSrz=zeros(5,5832);
sumvET=zeros(5,5832);
sumvQ=zeros(5,5832);
for i=1:4
    command=['psd' num2str(i) '=permute(enkf' num2str(i) '.STATE_maps.SD, [4,3,1,2]);'];
    eval(command)
    command=['psd' num2str(i) '=psd' num2str(i) '(:,:,pixel_indices);'];
    eval(command)
    command=['mpsd' num2str(i) '=mean(psd' num2str(i) ',3);'];
    eval(command)
    command=['psrz' num2str(i) '=permute(enkf' num2str(i) '.STATE_maps.Srz, [4,3,1,2]);'];
    eval(command)
    command=['psrz' num2str(i) '=psrz' num2str(i) '(:,:,pixel_indices);']
    eval(command)
    command=['mpsrz' num2str(i) '=mean(psrz' num2str(i) ',3);'];
    eval(command)
    command=['pet' num2str(i) '=permute(enkf' num2str(i) '.FLUX_maps.ET, [4,3,1,2]);'];
    eval(command)
    command=['pet' num2str(i) '=pet' num2str(i) '(:,:,pixel_indices);']
    eval(command)
    command=['mpet' num2str(i) '=mean(pet' num2str(i) ',3);'];
    eval(command)
    command=['mpq' num2str(i) '=permute(enkf' num2str(i) '.FLUX_time_series.outlet_hydrograph, [2,1]);'];
    eval(command)
end
for i=1:20
    sumvSD(1,:)=sumvSD(1,:)+(opentsSD.tsSD(i,2:5833)-meanSD(2:5833)).^2;
    sumvSrz(1,:)=sumvSrz(1,:)+(opentsSrz.tsSrz(i,2:5833)-meanSrz(2:5833)).^2;
    sumvET(1,:)=sumvET(1,:)+(opentsET.tsET(i,:)-meanET).^2;
    sumvQ(1,:)=sumvQ(1,:)+(opentsQ.tsQ(i,:)-meanQ).^2;
    for j=2:5
        %sumvSD(j,:)=sumvSD(j,:)+(mpsd1(i,:)-mean(mpsd1,1)).^2;
        command=['sumvSD(' num2str(j) ',:)=sumvSD(' num2str(j) ',:)+(mpsd' num2str(j-1) '(' num2str(i) ',:)-mean(mpsd' num2str(j-1) ',1)).^2;'];
        eval(command)
        command=['sumvSrz(' num2str(j) ',:)=sumvSrz(' num2str(j) ',:)+(mpsrz' num2str(j-1) '(' num2str(i) ',:)-mean(mpsrz' num2str(j-1) ',1)).^2;'];
        eval(command)
        command=['sumvET(' num2str(j) ',:)=sumvET(' num2str(j) ',:)+(mpet' num2str(j-1) '(' num2str(i) ',:)-mean(mpet' num2str(j-1) ',1)).^2;'];
        eval(command)
        command=['sumvQ(' num2str(j) ',:)=sumvQ(' num2str(j) ',:)+(mpq' num2str(j-1) '(' num2str(i) ',:)-mean(mpq' num2str(j-1) ',1)).^2;'];
        eval(command)
    end
end
varSD=sumvSD/20;
varSrz=sumvSrz/20;
varET=sumvET/20;
varQ=sumvQ/20;

figure(5);
subplot(2,2,1);
plot(DOY,squeeze(varSD(1,:)),'--r','LineWidth',1);hold on
plot(DOY,squeeze(varSD(2,:)),'--b','LineWidth',1);hold on
plot(DOY,squeeze(varSD(3,:)),'--g','LineWidth',1);hold on
plot(DOY,squeeze(varSD(4,:)),'--k','LineWidth',1);hold on
plot(DOY,squeeze(varSD(5,:)),'--m','LineWidth',1);hold on
ylabel('Variance SD');
set(gca,'FontSize',8);grid
legend('prior','post 0.002','post 0.02','post 0.05','post 0.1');
%legend('post 0.002','post 0.02','post 0.05','post 0.1');
subplot(2,2,2);
plot(DOY,squeeze(varSrz(1,:)),'--r','LineWidth',1);hold on
plot(DOY,squeeze(varSrz(2,:)),'--b','LineWidth',1);hold on
plot(DOY,squeeze(varSrz(3,:)),'--g','LineWidth',1);hold on
plot(DOY,squeeze(varSrz(4,:)),'--k','LineWidth',1);hold on
plot(DOY,squeeze(varSrz(5,:)),'--m','LineWidth',1);hold on
ylabel('Variance Srz');
set(gca,'FontSize',8);grid
legend('prior','post 0.002','post 0.02','post 0.05','post 0.1');
subplot(2,2,3);
plot(DOY,squeeze(varET(1,:)),'--r','LineWidth',1);hold on
plot(DOY,squeeze(varET(2,:)),'--b','LineWidth',1);hold on
plot(DOY,squeeze(varET(3,:)),'--g','LineWidth',1);hold on
plot(DOY,squeeze(varET(4,:)),'--k','LineWidth',1);hold on
plot(DOY,squeeze(varET(5,:)),'--m','LineWidth',1);hold on
ylabel('Variance ET');
set(gca,'FontSize',8);grid
legend('prior','post 0.002','post 0.02','post 0.05','post 0.1');
%legend('post 0.002','post 0.02','post 0.05','post 0.1');
subplot(2,2,4);
plot(DOY,squeeze(varQ(1,:)),'--r','LineWidth',1);hold on
plot(DOY,squeeze(varQ(2,:)),'--b','LineWidth',1);hold on
plot(DOY,squeeze(varQ(3,:)),'--g','LineWidth',1);hold on
plot(DOY,squeeze(varQ(4,:)),'--k','LineWidth',1);hold on
plot(DOY,squeeze(varQ(5,:)),'--m','LineWidth',1);hold on
ylabel('Variance Q');
set(gca,'FontSize',8);grid
legend('prior','post 0.002','post 0.02','post 0.05','post 0.1');
%legend('post 0.002','post 0.02','post 0.05','post 0.1');
cd('C:\Users\BANZH\Downloads\2019spring\CEE Steve\Project_outputs\Calculate_metrics_long\Figures_weekly');
saveas(figure(5),[pwd '/Variance.png']);        
close(figure(5));
%% Calculate Bias prior and posterior and plot
bsd=zeros(5,5832);
bsrz=zeros(5,5832);
bet=zeros(5,5832);
bq=zeros(5,5832);
mosd=mean(opentsSD.tsSD,1);
mosrz=mean(opentsSrz.tsSrz,1);
moet=mean(opentsET.tsET,1);
moq=mean(opentsQ.tsQ,1);
bsd(1,:)=mosd(2:5833)-truth1.states.time_series.SD(2:5833);
bsrz(1,:)=mosrz(2:5833)-truth1.states.time_series.Srz(2:5833);
bet(1,:)=moet-truth1.fluxes.time_series.ET;
bq(1,:)=moq-truth1.fluxes.time_series.outlet_hydrograph;
for i=1:4
    command=['bsd(' num2str(i+1) ',:)=mean(mpsd' num2str(i) ',1)-truth' num2str(i) '.states.time_series.SD(2:5833);'];
    eval(command)
    command=['bsrz(' num2str(i+1) ',:)=mean(mpsrz' num2str(i) ',1)-truth' num2str(i) '.states.time_series.Srz(2:5833);'];
    eval(command)
    command=['bet(' num2str(i+1) ',:)=mean(mpet' num2str(i) ',1)-truth' num2str(i) '.fluxes.time_series.ET;'];
    eval(command)
    command=['bq(' num2str(i+1) ',:)=mean(mpq' num2str(i) ',1)-truth' num2str(i) '.fluxes.time_series.outlet_hydrograph;'];
    eval(command)
end
figure(6);
subplot(2,2,1);
plot(DOY,squeeze(bsd(1,:)),'--r','LineWidth',1);hold on
plot(DOY,squeeze(bsd(2,:)),'--b','LineWidth',1);hold on
plot(DOY,squeeze(bsd(3,:)),'--g','LineWidth',1);hold on
plot(DOY,squeeze(bsd(4,:)),'--k','LineWidth',1);hold on
plot(DOY,squeeze(bsd(5,:)),'--m','LineWidth',1);hold on
ylabel('Bias SD');
set(gca,'FontSize',11);grid
legend('prior','post 0.002','post 0.02','post 0.05','post 0.1');
%legend('post 0.002','post 0.02','post 0.05','post 0.1');
subplot(2,2,2);
plot(DOY,squeeze(bsrz(1,:)),'--r','LineWidth',1);hold on
plot(DOY,squeeze(bsrz(2,:)),'--b','LineWidth',1);hold on
plot(DOY,squeeze(bsrz(3,:)),'--g','LineWidth',1);hold on
plot(DOY,squeeze(bsrz(4,:)),'--k','LineWidth',1);hold on
plot(DOY,squeeze(bsrz(5,:)),'--m','LineWidth',1);hold on
ylabel('Bias Srz');
set(gca,'FontSize',11);grid
legend('prior','post 0.002','post 0.02','post 0.05','post 0.1');
subplot(2,2,3);
plot(DOY,squeeze(bet(1,:)),'--r','LineWidth',1);hold on
plot(DOY,squeeze(bet(2,:)),'--b','LineWidth',1);hold on
plot(DOY,squeeze(bet(3,:)),'--g','LineWidth',1);hold on
plot(DOY,squeeze(bet(4,:)),'--k','LineWidth',1);hold on
plot(DOY,squeeze(bet(5,:)),'--m','LineWidth',1);hold on
ylabel('Bias ET');
set(gca,'FontSize',11);grid
legend('prior','post 0.002','post 0.02','post 0.05','post 0.1');
%legend('post 0.002','post 0.02','post 0.05','post 0.1');
subplot(2,2,4);
plot(DOY,squeeze(bq(1,:)),'--r','LineWidth',1);hold on
plot(DOY,squeeze(bq(2,:)),'--b','LineWidth',1);hold on
plot(DOY,squeeze(bq(3,:)),'--g','LineWidth',1);hold on
plot(DOY,squeeze(bq(4,:)),'--k','LineWidth',1);hold on
plot(DOY,squeeze(bq(5,:)),'--m','LineWidth',1);hold on
ylabel('Bias Q');
set(gca,'FontSize',11);grid
legend('prior','post 0.002','post 0.02','post 0.05','post 0.1');
%legend('post 0.002','post 0.02','post 0.05','post 0.1');
cd('C:\Users\BANZH\Downloads\2019spring\CEE Steve\Project_outputs\Calculate_metrics_long\Figures_weekly');
saveas(figure(6),[pwd '/Bias.png']);        
close(figure(6));
%% End
