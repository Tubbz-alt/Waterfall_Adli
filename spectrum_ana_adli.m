clear all;
set(gcf, 'Color', 'w');
set(0, 'defaultaxesfontsize', 20);
addpath('/Users/eadli/Dropbox/SLAC/E200/E200_waterfall/');
addpath('/Users/eadli/Dropbox/SLAC/E200/framework/E200_data');
addpath('/Users/eadli/Dropbox/SLAC/E200/E200_cher');

select_roi = 0;

header = '/Volumes/PWFA_4big';
%header = ''; % use no header for MCC operation
nas  ='/nas/nas-li20-pm00/';
expt = 'E200';
year = '/2014/';
%day  = '20140428/';
%dataset = '12722';
%day  = '20140425/';
%dataset = '12664';

% 5e16/cm3
day  = '20140424/';  dataset = '12615';  % Witness left edge at 435 pixels on SYAG.
day  = '20140424/';  dataset = '12614';  % Witness left edge at 550 pixels on SYAG.
day  = '20140424/';  dataset = '12616';  % Witness left edge at 645 pixels on SYAG.
day  = '20140424/'; dataset = '12582'; % Oven IN. Laser 1/2 Hz. Notch out. Left jaw in. Right jaw out. QS scan.
day  = '20140424/'; dataset = '12581'; % Oven IN. Laser 1/2 Hz. Notch and Jaw in. QS scan.


% 8e16/cm3
day  = '20140425/'; dataset = '12661'; % Oven in. Laser 1 Hz. Nominal notch and jaw positions. Phase ramp -17.6. QS scan
day  = '20140425/'; dataset = '12664'; % Oven in. Laser 1/2 Hz. Notch out. Phase ramp -17.6. QS scan.

% 3e16/cm3
day  = '20140427/'; dataset = '12700'; % QS Scan, Oven In, Laser On 0.5 Hz, Witness slit at pix 460.
day  = '20140427/'; dataset = '12692'; % Oven IN. Laser 1/2 Hz. Notch and Jaw in. QS scan.
day  = '20140427/'; dataset = '12693'; % Oven IN. Laser 1/2 Hz. Notch OUT. QS scan.
day  = '20140427/'; dataset = '12695'; % QS Scan, Oven In, Laser On 0.5 Hz, Witness only
day  = '20140427/'; dataset = '12698'; % Oven IN. Laser 1/2 Hz. Notch and Jaw in. QS scan.


if( strcmp(dataset, '12700') )
  Eaxismin = 20;
  Eaxismax = 40;
  caxismin = 50;
  caxismax = 500;
else
  Eaxismin = 20;
  Eaxismax = 31;
  caxismin = 50;
  caxismax = 200;
end% if


data_path = [nas expt year day expt '_' dataset '/' expt '_' dataset '.mat'];

% load data
load([header data_path]);

%IP2B = data.raw.images.IP2B;
%IP2B_bg = load([header IP2B.background_dat{1}]);

SYAG = data.raw.images.SYAG;
SYAG_bg = load([header SYAG.background_dat{1}]);

if( isfield( data.raw.images, 'B6') )
  have_B6_laser_data = 1;
else
  have_B6_laser_data = 0;
end% if

% TEMP, do not care about B6
have_B6_laser_data = 0;


if( have_B6_laser_data )
  B6 = data.raw.images.B6;
  B6_bg = load([header B6.background_dat{1}]);
end% if
  
CMOS_NEAR = data.raw.images.CMOS_NEAR;
CMOS_NEAR_bg = load([header CMOS_NEAR.background_dat{1}]);
%CN_bg_img = fliplr(CMOS_bg.img);

CMOS_FAR = data.raw.images.CMOS_FAR;
CMOS_FAR_bg = load([header CMOS_FAR.background_dat{1}]);
Eaxis_CF = E200_cher_get_E_axis('20131116', 'CMOS', 0, 1:2559, 0, 20.35);
esub = fliplr(Eaxis_CF);



EPICS_UID = data.raw.scalars.PATT_SYS1_1_PULSEID.UID;
[dummy,EPICS_CMOS_NEAR,CMOS_NEAR_index] = intersect(EPICS_UID,CMOS_NEAR.UID);
[dummy,EPICS_CMOS_FAR,CMOS_FAR_index] = intersect(EPICS_UID,CMOS_FAR.UID);
[dummy,EPICS_SYAG,SYAG_index] = intersect(EPICS_UID,SYAG.UID);
if( have_B6_laser_data )
  [dummy,EPICS_B6,B6_index] = intersect(EPICS_UID,B6.UID);
  n_b6 = numel(B6_index);
end% if

% EA: cut on pyro
do_cut_pyro = 0;
pyro_struct = data.raw.scalars.BLEN_LI20_3014_BRAW;
[psort, psort_idx] = sort(pyro_struct.dat);
if(do_cut_pyro)
  disp('EA: warning, cutting on pyro');
  figure(1);
  hh = plot(1:numel(psort), psort, 'b');
  set(hh, 'LineWidth', 3);
  hold on;
  xlabel('shot #');
  ylabel('pyro 3014');
  title(['Dataset ' dataset '     '], 'fontsize',16);
%   for n=2:N,
%    psort(round(end*(n-1)/N))
%   end% if
  N = 3;
  n = 3;
  n_min =round(numel(psort)*(n-1)/N+1);
  n_max =round( min(numel(psort)*(n)/N+1, numel(psort)) );
  pyro_min = psort(n_min);
  pyro_max = psort(n_max);
  pyro_idx = (pyro_struct.dat > pyro_min) &  (pyro_struct.dat < pyro_max);
  [dummy,pyro_idx2,CMOS_FAR_index] = intersect(pyro_struct.UID(pyro_idx),CMOS_FAR.UID);
  hh = plot(n_min:n_max, psort(n_min:n_max), 'r');
  set(hh, 'LineWidth', 3);
  hold off;
else
  figure(1);
  hh = plot(1:numel(psort), psort, 'r');
  set(hh, 'LineWidth', 3);
  xlabel('shot #');
  ylabel('pyro 3014');
  title(['Dataset ' dataset '     '], 'fontsize',16);
end% if

n_CMOS_NEAR = numel(CMOS_NEAR_index);
n_CMOS_FAR = numel(CMOS_FAR_index);
n_syag = numel(SYAG_index);




toro_2452_tmit = data.raw.scalars.GADC0_LI20_EX01_AI_CH0_.dat;
toro_3163_tmit = data.raw.scalars.GADC0_LI20_EX01_AI_CH2_.dat;
toro_3255_tmit = data.raw.scalars.GADC0_LI20_EX01_AI_CH3_.dat;
bpms_2445_x    = data.raw.scalars.BPMS_LI20_2445_X.dat;
bpms_2445_y    = data.raw.scalars.BPMS_LI20_2445_Y.dat;
bpms_2445_tmit = data.raw.scalars.BPMS_LI20_2445_TMIT.dat;
bpms_3156_x    = data.raw.scalars.BPMS_LI20_3156_X.dat;
bpms_3156_y    = data.raw.scalars.BPMS_LI20_3156_Y.dat;
bpms_3156_tmit = data.raw.scalars.BPMS_LI20_3156_TMIT.dat;
bpms_3265_x    = data.raw.scalars.BPMS_LI20_3265_X.dat;
bpms_3265_y    = data.raw.scalars.BPMS_LI20_3265_Y.dat;
bpms_3265_tmit = data.raw.scalars.BPMS_LI20_3265_TMIT.dat;
bpms_3315_x    = data.raw.scalars.BPMS_LI20_3315_X.dat;
bpms_3315_y    = data.raw.scalars.BPMS_LI20_3315_Y.dat;
bpms_3315_tmit = data.raw.scalars.BPMS_LI20_3315_TMIT.dat;
pyro           = data.raw.scalars.BLEN_LI20_3014_BRAW.dat;
laser          = data.raw.scalars.PMTR_LA20_10_PWR.dat;
laser_on = laser > 5;
laser_off = laser < 5;

if isfield(data.raw.metadata,'n_steps')
    is_scan = 1;
    
    n_step         = data.raw.metadata.n_steps;
    step_num       = data.raw.scalars.step_num.dat;
    step_val       = data.raw.scalars.step_value.dat;
    
    step_num_cmos_near = step_num(CMOS_NEAR_index);
    step_num_cmos_far = step_num(CMOS_FAR_index);
    step_num_syag = step_num(SYAG_index);
    
    step_val_cmos_near = step_val(CMOS_NEAR_index);
    step_val_cmos_far = step_val(CMOS_FAR_index);
    step_val_syag = step_val(SYAG_index);
    
else
    is_scan = 0;
    
    n_step = 1;
    step_num       = data.raw.scalars.step_num.dat;
    step_val       = data.raw.scalars.step_num.dat;

    step_val_cmos_near = step_val(CMOS_NEAR_index);
    step_val_cmos_far = step_val(CMOS_FAR_index);
    step_val_syag = step_val(SYAG_index);
    
    step_num_cmos_near = step_num(CMOS_NEAR_index);
    step_num_cmos_far = step_num(CMOS_FAR_index);
    step_num_syag = step_num(SYAG_index);
    
end

% syag_roi.top = 175;
% syag_roi.bottom = 225;
% syag_roi.left = 300;
% syag_roi.right = 1100;
% syag_roi.rot = 2;
% SYAG_ANA = basic_image_ana(SYAG,0,syag_roi,header);
% wit = sum(SYAG_ANA.x_profs(95:310,:));
% drive = sum(SYAG_ANA.x_profs(365:755,:));
% SYAG_ANA.sum = sum(SYAG_ANA.x_profs);
% SYAG_ANA.witness_charge = wit;%./SYAG_ANA.sum;
% SYAG_ANA.drive_charge = drive;%./SYAG_ANA.sum;

% roi.top = 200;
% roi.bottom = 730;
% roi.left = 350;
% roi.right = 850;
% roi.rot = 0;
%CMOS_NEAR_ANA = basic_image_ana(CMOS_NEAR,0,roi,header);

%%

if( strcmp(dataset, '12616')  ||  strcmp(dataset, '12615') ||  strcmp(dataset, '12614') || ~select_roi  )
     disp('EA: warning, ROI bybassed for this dataset');
     my_roi = [366 85; 571 979];
else
if select_roi || ~exist('my_roi')
    image = imread([header CMOS_FAR.dat{40}]);
    image = rot90(image,3);
    figure(1); imagesc(image); caxis([0 1000]);
    my_roi = ginput(2);
end
end% if

cmos_roi.top = round(my_roi(1,2));
cmos_roi.bottom = round(my_roi(2,2));
cmos_roi.left = round(my_roi(1,1));
cmos_roi.right = round(my_roi(2,1));

esub = esub(cmos_roi.top:cmos_roi.bottom);
cmos_roi.rot = 3;
use_bg = 1;
CMOS_FAR_ANA = basic_image_ana2(CMOS_FAR,use_bg,cmos_roi,header, CMOS_FAR_index);

qs_x = 1:n_CMOS_FAR;
qs_l = 1:n_CMOS_FAR/2;
%cmos_specs = CMOS_FAR_ANA.y_profs(:,CMOS_FAR_index);
cmos_specs = CMOS_FAR_ANA.y_profs;




%%
if( have_B6_laser_data )
  b6_roi.top = 1;
  b6_roi.bottom = 700;
  b6_roi.left = 1;
  b6_roi.right = 1200;
  b6_roi.rot = 0;
  B6_ANA = basic_image_ana(B6,0,b6_roi,header);

%%

b6_laser_on = B6_ANA.sum > 1E8;
b6_laser_off = B6_ANA.sum < 1E8;
% below: make intersect of UIDs
%cmos_on = cmos_specs(:,b6_laser_on);
%cmos_off = cmos_specs(:,b6_laser_off);
% temp hack
cmos_on = cmos_specs;
cmos_off = cmos_specs;
qs_x = 1:n_CMOS_FAR;
qs_l = 1:n_CMOS_FAR;

end% if( have_B6_laser_data )


%
%
% plotting
%
%

%
figure(2); pcolor(qs_x,esub,cmos_specs); shading flat; box off; colorbar;
if is_scan
    xlabel('Imaging Energy Relative to 20.35 GeV [GeV]','fontsize',16);
    steps = unique(step_num_cmos_far);
    n_steps = numel(steps);
    for n=1:n_steps,
      steps_val(n) = data.raw.metadata.param.dat{1}.PV_scan_list(n);
      steps_idx(n) = min(find(step_num_cmos_far == steps(n)));
    end% for
    set(gca, 'XTick', [0 steps_idx]);
    set(gca, 'XTickLabel', [0 steps_val]);
    for n=2:n_steps,
      hh = line([steps_idx(n) steps_idx(n)+eps],[esub(1) esub(end)],'color','k','linestyle','--');
      set(hh, 'LineWidth', 2.5);
    end
end

cmap  = custom_cmap();
colormap(cmap.wbgyr);

ylim([Eaxismin Eaxismax]);
caxis([caxismin caxismax]);

if is_scan 
    xlabel('Imaging Energy Relative to 20.35 GeV [GeV]','fontsize',16);
else
    xlabel('Shot Number','fontsize',16);
end

ylabel('Energy [GeV]','fontsize',16);
title(['Dataset ' dataset '     '], 'fontsize',16);
set(gca,'fontsize',16);

if( have_B6_laser_data )

figure(6); pcolor(qs_l,esub,cmos_on); shading flat; box off;
if is_scan
    xlabel('Imaging Energy Relative to 20.35 GeV [GeV]','fontsize',16);
    set(gca, 'XTick', (data.raw.metadata.param.dat{1}.n_shot/4):(data.raw.metadata.param.dat{1}.n_shot/2):(n_CMOS_FAR-data.raw.metadata.param.dat{1}.n_shot/4));
    set(gca, 'XTickLabel', data.raw.metadata.param.dat{1}.PV_scan_list);
    for i = 1:(n_step-1)
        line([i*(data.raw.metadata.param.dat{1}.n_shot/2) i*(data.raw.metadata.param.dat{1}.n_shot/2)],[esub(1) esub(end)],'color','k','linestyle','--');
    end
end

cmap  = custom_cmap();
colormap(cmap.wbgyr);
caxis([caxismin caxismax]);

if is_scan 
    xlabel('Imaging Energy Relative to 20.35 GeV [GeV]','fontsize',16);
else
    xlabel('Shot Number','fontsize',16);
end

ylabel('Energy [GeV]','fontsize',16);
title(['Dataset ' dataset '. Plasma in. Laser on shots.'],'fontsize',16);
set(gca,'fontsize',16);

figure(7); pcolor(qs_l,esub,cmos_off); shading flat; box off;
if is_scan
    xlabel('Imaging Energy Relative to 20.35 GeV [GeV]','fontsize',16);
    set(gca, 'XTick', (data.raw.metadata.param.dat{1}.n_shot/4):(data.raw.metadata.param.dat{1}.n_shot/2):(n_CMOS_FAR-data.raw.metadata.param.dat{1}.n_shot/4));
    set(gca, 'XTickLabel', data.raw.metadata.param.dat{1}.PV_scan_list);
    for i = 1:(n_step-1)
        line([i*(data.raw.metadata.param.dat{1}.n_shot/2) i*(data.raw.metadata.param.dat{1}.n_shot/2)],[esub(1) esub(end)],'color','k','linestyle','--');
    end
end

cmap  = custom_cmap();
colormap(cmap.wbgyr);
caxis([caxismin caxismax]);

if is_scan 
    xlabel('Imaging Energy Relative to 20.35 GeV [GeV]','fontsize',16);
else
    xlabel('Shot Number','fontsize',16);
end

ylabel('Energy [GeV]','fontsize',16);
title(['Dataset ' dataset '. Plasma in. Laser off shots.'],'fontsize',16);
set(gca,'fontsize',16);

end% if( have_B6_laser_data )


%
% AVG CHARGE VS QS
%

%%
if is_scan

figure(3);
specs = CMOS_FAR_ANA.y_profs(:,:);
if( have_B6_laser_data )
  on_specs = specs(:,b6_laser_on);
  off_specs = specs(:,b6_laser_on);
end% if( have_B6_laser_data )
mean_specs = zeros(numel(esub),n_step+1);
%mean_on = zeros(numel(esub),n_step+1);
%mean_off = zeros(numel(esub),n_step+1);
%dns_on = zeros(1,numel(CMOS_FAR_index)/2);
%dns_off = zeros(1,numel(CMOS_FAR_index)/2);
stds = zeros(1,n_step);
%stds_on = zeros(1,n_step);
%stds_off = zeros(1,n_step);
means = zeros(1,n_step);
%means_on = zeros(1,n_step);
%means_off = zeros(1,n_step);
for i = 1:n_step
    
    step_specs = specs(:,step_num_cmos_far==i);
if( have_B6_laser_data )
    on = b6_laser_on(step_num_cmos_far==i);
    off = b6_laser_on(step_num_cmos_far==i);
end% if

    dGeV = 0.2;
    hi = 20.35 + data.raw.metadata.param.dat{1}.PV_scan_list(i) + dGeV;
    lo = 20.35 + data.raw.metadata.param.dat{1}.PV_scan_list(i) - dGeV;
    range = esub > lo & esub < hi;
        
    dns(step_num_cmos_far==i) = mean(step_specs(range,:));
    stds(i) = std(dns(step_num_cmos_far==i));
    means(i) = mean(dns(step_num_cmos_far==i));
    mean_specs(:,i) = mean(step_specs,2);
    
    
end



pcolor(1:(n_step+1),esub,mean_specs); shading flat; box off; colorbar;
set(gca, 'XTick', (1:n_step)+0.5);
set(gca, 'XTickLabel', data.raw.metadata.param.dat{1}.PV_scan_list);
xlabel('Imaging Energy Relative to 20.35 GeV [GeV]','fontsize',16);
ylabel('Energy [GeV]','fontsize',16);
title(['Dataset ' dataset '     '], 'fontsize',16);

colormap(cmap.wbgyr);
caxis([caxismin caxismax]);
ylim([Eaxismin Eaxismax]);

end% if scan

if is_scan

figure(4);
plot(step_val_cmos_far,dns,'s');
xlabel('Imaging Energy Relative to 20.35 GeV [GeV]','fontsize',16);
ylabel('dN/dE [Counts/GeV]','fontsize',16);
ylim([0 caxismax*2]);
title(['Dataset ' dataset '     '], 'fontsize',16);
set(gca, 'XTick', steps_val);

end% if scan

if is_scan

figure(5);
errorbar(data.raw.metadata.param.dat{1}.PV_scan_list,means,stds,'s');
%axis([-1 data.raw.metadata.param.dat{1}.PV_scan_list(end)+1 0 500]);
xlabel('Imaging Energy Relative to 20.35 GeV [GeV]','fontsize',16);
ylabel('dN/dE [Counts/GeV]','fontsize',16);
ylim([0 caxismax]);
title(['Dataset ' dataset '     '], 'fontsize',16);
set(gca, 'XTick', steps_val);

end% if scan
