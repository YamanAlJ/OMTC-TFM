
%%  Section 1: Indices Filtering (May be necessary
clc;
clear; close all
set(0,'DefaultFigureWindowStyle','docked')

[text_only, numbers, text_cats] = xlsread('Exp_summary_vCompFinal.xlsx');

% Identify Cells 
RV1_indices = find(ismember(text_cats(:, 1), '22RV1')==1); 
PC3_indices = find(ismember(text_cats(:, 1), 'PC3')==1); 
DU145_indices = find(ismember(text_cats(:,1), 'DU145')==1);
LNCaP_indices = find(ismember(text_cats(:,1), 'LNCaP')==1); 

% Find indices for cell_type 
ex_soft_indices = find(ismember(text_cats(:, 2), '1 kPa')==1); 
soft_indices = find(ismember(text_cats(:, 2), '3 kPa')==1); 
stiff_indices = find(ismember(text_cats(:, 2), '12 kPa')==1); 
ex_stiff_indices = find(ismember(text_cats(:, 2), '50 kPa')==1); 


% Find indices for incubation subsets 
inc24hr_indices = find(ismember(text_cats(:, 4), '24hrs')==1);
inc72hr_indices = find(ismember(text_cats(:, 4), '72hrs')==1);
inc48hr_indices = find(ismember(text_cats(:, 4), '48hrs')==1);

% Find indices for TREATMENT subsets
% No_treatment_indices = find(ismember(text_cats(:, 6), 'No'));
% FBSmediaonly_indices = find(ismember(text_cats(:, 3), 'FBS') == 1);
% 
% CSSmediaonly_indices = find(ismember(text_cats(:,3), 'CSS') == 1); 
% CDXTreatment_indices = find(ismember(text_cats(:, 'fill'), 'CDX')); 
% FBSTreatment_indices = find(ismember(text_cats(:, 'fill'), 'FBS')); 
% 
% CellType and Stiffness, FBS Only 
% RV1_3kPa_ind = intersect(RV1_indices, soft_indices);
% RV1_12kPa_ind = intersect(RV1_indices, stiff_indices); 
% RV1_1kPa_ind = intersect(RV1_indices, ex_soft_indices);
% RV1_50kPa_stiff = intersect(RV1_indices, ex_stiff_indices);
% 
% Cell/Stiffness selected, 72 hours
% RV1_3kPa_72_ind = intersect(RV1_soft_ind, inc72hr_indices); 
% RV1_12kPa_72_ind = intersect(RV1_stiff_ind, inc72hr_indices); 
% 
% Cell/Stiffness/FBS selected, 72 hours
% RV1_3kPa_72_FBS_ind = intersect(FBSmediaonly_indices, RV1_3kPa_72_ind ); 
% RV1_12kPa_72_FBS_ind = intersect(FBSmediaonly_indices, RV1_12kPa_72_ind); 
% 
% % As above 
% PC3_3kPa_ind = intersect(PC3_indices, soft_indices);
% PC3_12kPa_ind = intersect(PC3_indices, stiff_indices); 
% PC3_1kPa_ind = intersect(PC3_indices, ex_soft_indices);
% PC3_50kPa_ind = intersect(PC3_indices, ex_stiff_indices); 
% 
% 
% PC3_3kPa_72_ind = intersect(PC3_3kPa_ind, inc72hr_indices); 
% PC3_12kPa_72_ind = intersect(PC3_12kPa_ind, inc72hr_indices); 
% 
% 


%% Want to do a bulk comparison of pre-treatment of different 

PC3_3kPa_ind = intersect(PC3_indices, soft_indices);
PC3_12kPa_ind = intersect(PC3_indices, stiff_indices); 
PC3_1kPa_ind = intersect(PC3_indices, ex_soft_indices);
PC3_50kPa_ind = intersect(PC3_indices, ex_stiff_indices);


% PC3_3kPa_72_ind = intersect(PC3_3kPa_ind, inc72hr_indices); 
% PC3_12kPa_72_ind = intersect(PC3_12kPa_ind, inc72hr_indices); 

% PC3_3kPa_72_FBS_ind = intersect(PC3_3kPa_72_ind, FBSmediaonly_indices); 
% PC3_3kPa_72_CSS_ind = intersect(PC3_3kPa_ind, CSSmediaonly_indices);

% PC3_12kPa_72_FBS_ind = intersect(PC3_12kPa_ind, FBSmediaonly_indices); 
% PC3_12kPa_72_CSS_ind = intersect(PC3_12kPa_ind, CSSmediaonly_indices);


RV1_3kPa_ind = intersect(RV1_indices, soft_indices);
RV1_12kPa_ind = intersect(RV1_indices, stiff_indices); 
RV1_1kPa_ind = intersect(RV1_indices, ex_soft_indices);
RV1_50kPa_ind = intersect(RV1_indices, ex_stiff_indices); 
% 
% RV1_3kPa_72_ind = intersect(RV1_3kPa_ind, inc72hr_indices); 
% RV1_12kPa_72_ind = intersect(RV1_12kPa_ind, inc72hr_indices); 
% 
% RV1_3kPa_72_FBS_ind = intersect(FBSmediaonly_indices, RV1_3kPa_72_ind); 
% RV1_12kPa_72_FBS_ind = intersect(FBSmediaonly_indices, RV1_12kPa_72_ind); 

% RV1_3kPa_72_CSS_ind = intersect(RV1_3kPa_ind, CSSmediaonly_indices);
% RV1_12kPa_72_CSS_ind = intersect(RV1_12kPa_ind, CSSmediaonly_indices);

DU145_3kPa_ind = intersect(DU145_indices, soft_indices);
DU145_12kPa_ind = intersect(DU145_indices, stiff_indices); 
DU145_1kPa_ind = intersect(DU145_indices, ex_soft_indices);
DU145_50kPa_ind = intersect(DU145_indices, ex_stiff_indices); 

% DU145_3kPa_72_ind = intersect(DU145_3kPa_ind, inc72hr_indices); 
% DU145_12kPa_72_ind = intersect(DU145_12kPa_ind, inc72hr_indices); 
% 
% DU145_3kPa_72_FBS_ind = intersect(FBSmediaonly_indices, DU145_3kPa_72_ind);
% DU145_12kPa_72_FBS_ind= intersect(FBSmediaonly_indices, DU145_12kPa_72_ind);

% DU145_3kPa_72_CSS_ind = intersect(DU145_3kPa_ind, CSSmediaonly_indices);
% DU145_12kPa_72_CSS_ind = intersect(DU145_12kPa_ind, CSSmediaonly_indices);

LNCaP_3kPa_ind = intersect(LNCaP_indices, soft_indices);
LNCaP_12kPa_ind = intersect(LNCaP_indices, stiff_indices); 
LNCaP_1kPa_ind = intersect(LNCaP_indices, ex_soft_indices);
LNCaP_50kPa_ind = intersect(LNCaP_indices, ex_stiff_indices); 

% LNCaP_3kPa_72_ind = intersect(LNCaP_3kPa_ind, inc72hr_indices); 
% LNCaP_12kPa_72_ind = intersect(LNCaP_12kPa_ind, inc72hr_indices); 
% 
% LNCaP_3kPa_72_FBS_ind = intersect(FBSmediaonly_indices, LNCaP_3kPa_72_ind);
% LNCaP_12kPa_72_FBS_ind = intersect(FBSmediaonly_indices, LNCaP_12kPa_72_ind);
% 
% LNCaP_3kPa_72_CSS_ind = intersect(LNCaP_3kPa_ind, CSSmediaonly_indices);
% LNCaP_12kPa_72_CSS_ind = intersect(LNCaP_12kPa_ind, CSSmediaonly_indices);
% 

%% SECTION 2: 22RV1 vs. PC3 vs. DU145 vs. LNCaP as a Function of   Stiffness 
f_stiff_rmst = figure('Name', 'RMST vs. Time, 22RV1 & PC3 12kPa');
ylabel('RMST (Pa)');
xlabel('Time (Minutes)'); 

f_stiff_SE = figure('Name', 'SE vs. Time 22RV1 vs. PC3 12kPa');
ylabel('Strain  Energy (pJ)');
xlabel('Time (Minutes)'); 
  
% Plot RMST and SE of 12 kPa 22RV1

% Identify Relevant Exps
to_open = text_cats(RV1_12kPa_ind, 9);

% Identify tie interval
time_interval = cell2mat(text_cats(RV1_12kPa_ind,10))/60;
frame_add = cell2mat(text_cats(RV1_12kPa_ind,6));

% Identify bad positions
series_omit_array = text_cats(RV1_12kPa_ind,7);

% Run through each position and plot 
q=1;

for i = 1:length(to_open) 
    % Open file and extract values rxc = pos x time
  data_temp = open([to_open{i}]);
  rmst = data_temp.rmst';
  strain_energy = data_temp.strain_energy'*10^12; % Adjust to pJ
  average_displacement = data_temp.average_displacement';
  
  % Identify positions to omit and do so
  if series_omit_array{i} ~=0 
    try 
    omit_indices = str2double(strsplit(series_omit_array{i}, ','));
    catch
    omit_indices = series_omit_array{i}; 
    end 
    rmst(omit_indices,:) = []; 
    strain_energy(omit_indices, :) = [];
    average_displacement(omit_indices,:) = [];
  end     

  % Identify the frame at which the treatment was applied 
  [pos_num, time_num] = size(rmst);
  if frame_add(i) ~= 0
      time_num = frame_add(i) - 1;
  end 
  t = ([1:time_num]-1)*time_interval(i); 

  for j = 1:pos_num  
  figure(f_stiff_rmst)  
  hold on
  p1 = plot(t, rmst(j,1:length(t)), 'r'); 
  time_aveRMST_12RV1(q) = mean(rmst(j, 1:length(t)));
% text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
  figure(f_stiff_SE) 
  hold on
  q1 =  plot(t, strain_energy(j,1:length(t)), 'r'); 
  time_aveSE_12RV1(q) = mean(strain_energy(j, 1:length(t)));
  time_aveAD_12RV1(q) = mean(average_displacement(j, 1:length(t)));
  
  q = q+1;
  end 
end 

% Plot RMST and SE of 12 kPa PC3
% Identify Relevant Exps
to_open = text_cats(PC3_12kPa_ind , 9);

% Identify time interval
time_interval = cell2mat(text_cats(PC3_12kPa_ind ,10))/60;

% Identify when the frame for adding 
frame_add = cell2mat(text_cats(PC3_12kPa_ind ,6));

% Identify bad positions
series_omit_array = text_cats(PC3_12kPa_ind ,7);

% Run through each position and plot 
q=1;
for i = 1:length(to_open) 
    % Open file and extract values rxc = pos x time
  data_temp = open(['Mats_So_Far/' to_open{i}]);
  rmst = data_temp.rmst';
  strain_energy = data_temp.strain_energy'*10^12;
  average_displacement = data_temp.average_displacement';
  
  % Identify positions to omit and do so
  if series_omit_array{i} ~=0 
    try 
    omit_indices = str2double(strsplit(series_omit_array{i}, ','));
    catch
    omit_indices = series_omit_array{i}; 
    end 
    rmst(omit_indices,:) = []; 
    strain_energy(omit_indices, :) = [];
    average_displacement(omit_indices,:) = [];
  end     



  [pos_num, time_num] = size(rmst);
  if frame_add(i) ~= 0
      time_num = frame_add(i) - 1;
  end 
  t = ([1:time_num]-1)*time_interval(i); 
  
  for j = 1:pos_num  
      figure(f_stiff_rmst)  
      hold on
      p2 = plot(t, rmst(j,1:length(t)), 'b');
      time_aveRMST_12PC3(q) = mean(rmst(j, 1:length(t)));
% text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
      figure(f_stiff_SE)
      hold on
      q2 = plot(t, strain_energy(j,1:length(t)), 'b');
      time_aveSE_12PC3(q) = mean(strain_energy(j, 1:length(t)));
        time_aveAD_12PC3(q) = mean(average_displacement(j, 1:length(t)));
  q = q+1;
  end 
end 
   
% Plot RMST and SE of 12 kPa DU145 in FBS (non-CSS) 
% Identify Relevant Exps
to_open = text_cats(DU145_12kPa_ind, 9);

% Identify time interval
time_interval = cell2mat(text_cats(DU145_12kPa_ind ,10))/60;

% Identify when the frame for adding 
frame_add = cell2mat(text_cats(DU145_12kPa_ind ,6));

% Identify bad positions
series_omit_array = text_cats(DU145_12kPa_ind ,7);

% Run through each position and plot 
q=1;
for i = 1:length(to_open) 
    % Open file and extract values rxc = pos x time
  data_temp = open(['Mats_So_Far/' to_open{i}]);
  rmst = data_temp.rmst';
  strain_energy = data_temp.strain_energy'*10^12;
  average_displacement = data_temp.average_displacement';
  
  % Identify positions to omit and do so
  if series_omit_array{i} ~=0 
    try 
    omit_indices = str2double(strsplit(series_omit_array{i}, ','));
    catch
    omit_indices = series_omit_array{i}; 
    end 
    rmst(omit_indices,:) = []; 
    strain_energy(omit_indices, :) = [];
    average_displacement(omit_indices,:) = [];
  end     


  [pos_num, time_num] = size(rmst);
  if frame_add(i) ~= 0
      time_num = frame_add(i) - 1;
  end 
  t = ([1:time_num]-1)*time_interval(i); 
  
  for j = 1:pos_num  
      figure(f_stiff_rmst)  
      hold on
      p3 = plot(t, rmst(j,1:length(t)), 'g');
      time_aveRMST_12DU145(q) = mean(rmst(j, 1:length(t)));
% text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
      figure(f_stiff_SE)
      hold on
      q3 = plot(t, strain_energy(j,1:length(t)), 'g');
      time_aveSE_12DU145(q) = mean(strain_energy(j, 1:length(t)));
      time_aveAD_12DU145(q) = mean(strain_energy(j, 1:length(t)));
  q = q+1;
  end 
  
end 

  
% Plot RMST and SE of 12 kPa LNCaP in FBS (non-CSS) 
% Identify Relevant Exps
to_open = text_cats(LNCaP_12kPa_ind, 9);

% Identify time interval
time_interval = cell2mat(text_cats(LNCaP_12kPa_ind ,10))/60;

% Identify when the frame for adding 
frame_add = cell2mat(text_cats(LNCaP_12kPa_ind ,6));

% Identify bad positions
series_omit_array = text_cats(LNCaP_12kPa_ind ,7);

% Run through each position and plot 
q=1;
for i = 1:length(to_open) 
    % Open file and extract values rxc = pos x time
  data_temp = open(['Mats_So_Far/' to_open{i}]);
  rmst = data_temp.rmst';
  strain_energy = data_temp.strain_energy'*10^12;
   average_displacement = data_temp.average_displacement';
  
  % Identify positions to omit and do so
  if series_omit_array{i} ~=0 
    try 
    omit_indices = str2double(strsplit(series_omit_array{i}, ','));
    catch
    omit_indices = series_omit_array{i}; 
    end 
    rmst(omit_indices,:) = []; 
    strain_energy(omit_indices, :) = [];
    average_displacement(omit_indices,:) = [];
  end     

  [pos_num, time_num] = size(rmst);
  if frame_add(i) ~= 0
      time_num = frame_add(i) - 1;
  end 
  t = ([1:time_num]-1)*time_interval(i); 
  
  for j = 1:pos_num  
      figure(f_stiff_rmst)  
      hold on
      p4 = plot(t, rmst(j,1:length(t)), 'm');
      time_aveRMST_12LNCaP(q) = mean(rmst(j, 1:length(t)));
% text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
      figure(f_stiff_SE)
      hold on
      q4 = plot(t, strain_energy(j,1:length(t)), 'm');
      time_aveSE_12LNCaP(q) = mean(strain_energy(j, 1:length(t)));
        time_aveAD_12LNCaP(q) = mean(average_displacement(j, 1:length(t)));
  q = q+1;
  end 
  
end 


hold off 
figure(f_stiff_rmst)
hold on
legend([p1,p4,p3,p2], '22RV1, 12 kPa','LNCaP, 12 kPa', 'DU145, 12 kPa' , 'PC3, 12 kPa')
saveas(gcf, 's12kPa_rmst.tiff');

figure(f_stiff_SE)
hold on
legend([q1,q4,q3,q2],'22RV1, 12 kPa','LNCaP, 12 kPa', 'DU145, 12 kPa' , 'PC3, 12 kPa')
saveas(gcf, 'stiff_SE.tiff');


%% RMST and Strain Energy Boxplots STIFF
% f_box_stiff_rmst = figure('Name', 'RMST boxplot 22RV1 vs. PC3 12kPa');
% ylabel('RMST (Pa)');
% 
% f_box_stiff_SE = figure('Name', 'SE boxplot 22RV1 vs. PC3 12kPa');
% ylabel('Strain  Energy (pJ)');
% 
% figure(f_box_stiff_rmst)
% hold on 
% boxplot([time_aveRMST_12RV1,time_aveRMST_12LNCaP ,time_aveRMST_12DU145, time_aveRMST_12PC3], [1*ones((length(time_aveRMST_12RV1)),1); 2*ones((length(time_aveRMST_12LNCaP)),1); ...
%           3*ones((length(time_aveRMST_12DU145)),1); 4*ones((length(time_aveRMST_12PC3)),1)], 'Positions',[1*ones((length(time_aveRMST_12RV1)),1); 2*ones((length(time_aveRMST_12LNCaP)),1); ...
%           3*ones((length(time_aveRMST_12DU145)),1); 4*ones((length(time_aveRMST_12PC3)),1)],  'Labels', {'22RV1 12 kPa', 'LNCaP', 'DU145', 'PC3 12 kPa'});
%       
% RMST_scatter_ind12RV1 = 1*ones((length(time_aveRMST_12RV1)),1);
% f1 = scatter(RMST_scatter_ind12RV1, time_aveRMST_12RV1, 'r');
% 
% RMST_scatter_ind12DU145 = 2*ones((length(time_aveRMST_12LNCaP)),1);
% f1 = scatter(RMST_scatter_ind12DU145, time_aveRMST_12LNCaP, 'r');
% 
% RMST_scatter_ind12DU145 = 3*ones((length(time_aveRMST_12DU145)),1);
% f1 = scatter(RMST_scatter_ind12DU145, time_aveRMST_12DU145, 'r');
% 
% RMST_scatter_ind12PC3 = 4*ones((length(time_aveRMST_12PC3)),1);
% f1 = scatter(RMST_scatter_ind12PC3, time_aveRMST_12PC3, 'b');
%       
% % saveas('f_box_stiff_rmst.tiff');
% 
% figure(f_box_stiff_SE)
% hold on 
% boxplot([time_aveSE_12RV1,time_aveSE_12LNCaP ,time_aveSE_12DU145, time_aveSE_12PC3], [1*ones((length(time_aveSE_12RV1)),1); 2*ones((length(time_aveSE_12LNCaP)),1); ...
%           3*ones((length(time_aveSE_12DU145)),1); 4*ones((length(time_aveSE_12PC3)),1)], 'Positions',[1*ones((length(time_aveSE_12RV1)),1); 2*ones((length(time_aveSE_12LNCaP)),1); ...
%           3*ones((length(time_aveSE_12DU145)),1); 4*ones((length(time_aveSE_12PC3)),1)],  'Labels', {'22RV1 12 kPa', 'LNCaP', 'DU145', 'PC3 12 kPa'});
%       
% SE_scatter_ind12RV1 = 1*ones((length(time_aveSE_12RV1)),1);
% f1 = scatter(SE_scatter_ind12RV1, time_aveSE_12RV1, 'r');
% 
% SE_scatter_ind12LNCaP = 2*ones((length(time_aveSE_12LNCaP)),1);
% f1 = scatter(SE_scatter_ind12LNCaP, time_aveSE_12LNCaP, 'r');
% 
% SE_scatter_ind12DU145 = 3*ones((length(time_aveSE_12DU145)),1);
% f1 = scatter(SE_scatter_ind12DU145, time_aveSE_12DU145, 'r');
% 
% SE_scatter_ind12PC3 = 4*ones((length(time_aveRMST_12PC3)),1);
% f1 = scatter(SE_scatter_ind12PC3, time_aveSE_12PC3, 'b');

% saveas(gcf, 'f_box_stiff_SE.tiff');
%% 3kPa TIME PLOTS 
f_soft_rmst = figure('Name', 'RMST vs. Time, 22RV1 & PC3 3kPa');
ylabel('RMST (Pa)');
xlabel('Time (Minutes)'); 

f_soft_SE = figure('Name', 'SE vs. Time 22RV1 vs. PC3 3kPa');
ylabel('Strain  Energy (pJ)');
xlabel('Time (Minutes)'); 
  
% Plot RMST and SE of 3 kPa 22RV1

% Identify Relevant Exps
to_open = text_cats(RV1_3kPa_ind, 9);

% Identify tie interval
time_interval = cell2mat(text_cats(RV1_3kPa_ind,10))/60;
frame_add = cell2mat(text_cats(RV1_3kPa_ind,6));

% Identify bad positions
series_omit_array = text_cats(RV1_3kPa_ind,7);

% Run through each position and plot 
q=1;

for i = 1:length(to_open) 
    % Open file and extract values rxc = pos x time
  data_temp = open(['Mats_So_Far/' to_open{i}]);
  rmst = data_temp.rmst';
  strain_energy = data_temp.strain_energy'*10^12; % Adjust to pJ
   average_displacement = data_temp.average_displacement';
  
  % Identify positions to omit and do so
  if series_omit_array{i} ~=0 
    try 
    omit_indices = str2double(strsplit(series_omit_array{i}, ','));
    catch
    omit_indices = series_omit_array{i}; 
    end 
    rmst(omit_indices,:) = []; 
    strain_energy(omit_indices, :) = [];
    average_displacement(omit_indices,:) = [];
  end     


  % Identify the frame at which the treatment was applied 
  [pos_num, time_num] = size(rmst);
  if frame_add(i) ~= 0
      time_num = frame_add(i) - 1;
  end 
  t = ([1:time_num]-1)*time_interval(i); 

  for j = 1:pos_num  
  figure(f_soft_rmst)  
  hold on
  p1 = plot(t, rmst(j,1:length(t)), 'r'); 
  time_aveRMST_3RV1(q) = mean(rmst(j, 1:length(t)));
% text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
  figure(f_soft_SE) 
  hold on
  q1 =  plot(t, strain_energy(j,1:length(t)), 'r'); 
  time_aveSE_3RV1(q) = mean(strain_energy(j, 1:length(t)));
    time_aveAD_3RV1(q) = mean(average_displacement(j, 1:length(t)));
  q = q+1;
  end 
end 

% Plot RMST and SE of 3 kPa PC3
% Identify Relevant Exps
to_open = text_cats(PC3_3kPa_ind , 9);

% Identify time interval
time_interval = cell2mat(text_cats(PC3_3kPa_ind ,10))/60;

% Identify when the frame for adding 
frame_add = cell2mat(text_cats(PC3_3kPa_ind ,6));

% Identify bad positions
series_omit_array = text_cats(PC3_3kPa_ind ,7);

% Run through each position and plot 
q=1;
for i = 1:length(to_open) 
    % Open file and extract values rxc = pos x time
  data_temp = open(['Mats_So_Far/' to_open{i}]);
  rmst = data_temp.rmst';
  strain_energy = data_temp.strain_energy'*10^12;
   average_displacement = data_temp.average_displacement';
  
  % Identify positions to omit and do so
  if series_omit_array{i} ~=0 
    try 
    omit_indices = str2double(strsplit(series_omit_array{i}, ','));
    catch
    omit_indices = series_omit_array{i}; 
    end 
    rmst(omit_indices,:) = []; 
    strain_energy(omit_indices, :) = [];
    average_displacement(omit_indices,:) = [];
  end     

  [pos_num, time_num] = size(rmst);
  if frame_add(i) ~= 0
      time_num = frame_add(i) - 1;
  end 
  t = ([1:time_num]-1)*time_interval(i); 
  
  for j = 1:pos_num  
      figure(f_soft_rmst)  
      hold on
      p2 = plot(t, rmst(j,1:length(t)), 'b');
      time_aveRMST_3PC3(q) = mean(rmst(j, 1:length(t)));
% text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
      figure(f_soft_SE)
      hold on
      q2 = plot(t, strain_energy(j,1:length(t)), 'b');
      time_aveSE_3PC3(q) = mean(strain_energy(j, 1:length(t)));
        time_aveAD_3PC3(q) = mean(average_displacement(j, 1:length(t)));
  q = q+1;
  end 
end 
   
% Plot RMST and SE of 3 kPa DU145 in FBS (non-CSS) 
% Identify Relevant Exps
to_open = text_cats(DU145_3kPa_ind, 9);

% Identify time interval
time_interval = cell2mat(text_cats(DU145_3kPa_ind ,10))/60;

% Identify when the frame for adding 
frame_add = cell2mat(text_cats(DU145_3kPa_ind ,6));

% Identify bad positions
series_omit_array = text_cats(DU145_3kPa_ind ,7);

% Run through each position and plot 
q=1;
for i = 1:length(to_open) 
    % Open file and extract values rxc = pos x time
  data_temp = open(['Mats_So_Far/' to_open{i}]);
  rmst = data_temp.rmst';
  strain_energy = data_temp.strain_energy'*10^12;
   average_displacement = data_temp.average_displacement';
  
  % Identify positions to omit and do so
  if series_omit_array{i} ~=0 
    try 
    omit_indices = str2double(strsplit(series_omit_array{i}, ','));
    catch
    omit_indices = series_omit_array{i}; 
    end 
    rmst(omit_indices,:) = []; 
    strain_energy(omit_indices, :) = [];
    average_displacement(omit_indices,:) = [];
  end     


  [pos_num, time_num] = size(rmst);
  if frame_add(i) ~= 0
      time_num = frame_add(i) - 1;
  end 
  t = ([1:time_num]-1)*time_interval(i); 
  
  for j = 1:pos_num  
      figure(f_soft_rmst)  
      hold on
      p3 = plot(t, rmst(j,1:length(t)), 'g');
      time_aveRMST_3DU145(q) = mean(rmst(j, 1:length(t)));
% text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
      figure(f_soft_SE)
      hold on
      q3 = plot(t, strain_energy(j,1:length(t)), 'g');
      time_aveSE_3DU145(q) = mean(strain_energy(j, 1:length(t)));
      time_aveAD_3DU145(q) = mean(average_displacement(j, 1:length(t)));
  q = q+1;
  end 
  
end 

  
% Plot RMST and SE of 3 kPa DU145 in FBS (non-CSS) 
% Identify Relevant Exps
to_open = text_cats(LNCaP_3kPa_ind, 9);

% Identify time interval
time_interval = cell2mat(text_cats(LNCaP_3kPa_ind ,10))/60;

% Identify when the frame for adding 
frame_add = cell2mat(text_cats(LNCaP_3kPa_ind ,6));

% Identify bad positions
series_omit_array = text_cats(LNCaP_3kPa_ind ,7);

% Run through each position and plot 
q=1;
for i = 1:length(to_open) 
    % Open file and extract values rxc = pos x time
  data_temp = open(['Mats_So_Far/' to_open{i}]);
  rmst = data_temp.rmst';
  strain_energy = data_temp.strain_energy'*10^12;
   average_displacement = data_temp.average_displacement';
  
  % Identify positions to omit and do so
  if series_omit_array{i} ~=0 
    try 
    omit_indices = str2double(strsplit(series_omit_array{i}, ','));
    catch
    omit_indices = series_omit_array{i}; 
    end 
    rmst(omit_indices,:) = []; 
    strain_energy(omit_indices, :) = [];
    average_displacement(omit_indices,:) = [];
  end     

  [pos_num, time_num] = size(rmst);
  if frame_add(i) ~= 0
      time_num = frame_add(i) - 1;
  end 
  t = ([1:time_num]-1)*time_interval(i); 
  
  for j = 1:pos_num  
      figure(f_soft_rmst)  
      hold on
      p4 = plot(t, rmst(j,1:length(t)), 'm');
      time_aveRMST_3LNCaP(q) = mean(rmst(j, 1:length(t)));
% text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
      figure(f_soft_SE)
      hold on
      q4 = plot(t, strain_energy(j,1:length(t)), 'm');575
      time_aveSE_3LNCaP(q) = mean(strain_energy(j, 1:length(t)));
       time_aveAD_3LNCaP(q) = mean(average_displacement(j, 1:length(t)));
  q = q+1;
  end 
  
end 


hold off 
figure(f_soft_rmst)
hold on
legend([p1,p4,p3,p2],'22RV1, 3 kPa','LNCaP, 3 kPa', 'DU145, 3 kPa' , 'PC3, 3 kPa')
saveas(gcf, 'soft3_rmst.tiff');

figure(f_soft_SE)
hold on
legend([q1,q4,q3,q2],'22RV1, 3 kPa','LNCaP, 3 kPa', 'DU145, 3 kPa' , 'PC3, 3 kPa')
saveas(gcf, 'soft3_SE.tiff');

%% 1 kPa 

f_soft_rmst = figure('Name', 'RMST vs. Time, 22RV1 & PC3 1kPa');
ylabel('RMST (Pa)');
xlabel('Time (Minutes)'); 

f_soft_SE = figure('Name', 'SE vs. Time 22RV1 vs. PC3 1kPa');
ylabel('Strain  Energy (pJ)');
xlabel('Time (Minutes)'); 
  
% Plot RMST and SE of 3 kPa 22RV1

% Identify Relevant Exps
to_open = text_cats(RV1_1kPa_ind, 9);

% Identify tie interval
time_interval = cell2mat(text_cats(RV1_1kPa_ind,10))/60;
frame_add = cell2mat(text_cats(RV1_1kPa_ind,6));

% Identify bad positions
series_omit_array = text_cats(RV1_1kPa_ind,7);

% Run through each position and plot 
q=1;

for i = 1:length(to_open) 
    % Open file and extract values rxc = pos x time
  data_temp = open(['Mats_So_Far/' to_open{i}]);
  rmst = data_temp.rmst';
  strain_energy = data_temp.strain_energy'*10^12; % Adjust to pJ
   average_displacement = data_temp.average_displacement';
  
  % Identify positions to omit and do so
  if series_omit_array{i} ~=0 
    try 
    omit_indices = str2double(strsplit(series_omit_array{i}, ','));
    catch
    omit_indices = series_omit_array{i}; 
    end 
    rmst(omit_indices,:) = []; 
    strain_energy(omit_indices, :) = [];
    average_displacement(omit_indices,:) = [];
  end     

  % Identify the frame at which the treatment was applied 
  [pos_num, time_num] = size(rmst);
  if frame_add(i) ~= 0
      time_num = frame_add(i) - 1;
  end 
  t = ([1:time_num]-1)*time_interval(i); 

  for j = 1:pos_num  
  figure(f_soft_rmst)  
  hold on
  p1 = plot(t, rmst(j,1:length(t)), 'r'); 
  time_aveRMST_1RV1(q) = mean(rmst(j, 1:length(t)));
% text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
  figure(f_soft_SE) 
  hold on
  q1 =  plot(t, strain_energy(j,1:length(t)), 'r'); 
  time_aveSE_1RV1(q) = mean(strain_energy(j, 1:length(t)));
    time_aveAD_1RV1(q) = mean(average_displacement(j, 1:length(t)));
  q = q+1;
  end 
end 

% Plot RMST and SE of 3 kPa PC3
% Identify Relevant Exps
to_open = text_cats(PC3_1kPa_ind , 9);

% Identify time interval
time_interval = cell2mat(text_cats(PC3_1kPa_ind ,10))/60;

% Identify when the frame for adding 
frame_add = cell2mat(text_cats(PC3_1kPa_ind ,6));

% Identify bad positions
series_omit_array = text_cats(PC3_1kPa_ind ,7);

% Run through each position and plot 
q=1;
for i = 1:length(to_open) 
    % Open file and extract values rxc = pos x time
  data_temp = open(['Mats_So_Far/' to_open{i}]);
  rmst = data_temp.rmst';
  strain_energy = data_temp.strain_energy'*10^12;
   average_displacement = data_temp.average_displacement';
  
  % Identify positions to omit and do so
  if series_omit_array{i} ~=0 
    try 
    omit_indices = str2double(strsplit(series_omit_array{i}, ','));
    catch
    omit_indices = series_omit_array{i}; 
    end 
    rmst(omit_indices,:) = []; 
    strain_energy(omit_indices, :) = [];
    average_displacement(omit_indices,:) = [];
  end     

  [pos_num, time_num] = size(rmst);
  if frame_add(i) ~= 0
      time_num = frame_add(i) - 1;
  end 
  t = ([1:time_num]-1)*time_interval(i); 
  
  for j = 1:pos_num  
      figure(f_soft_rmst)  
      hold on
      p2 = plot(t, rmst(j,1:length(t)), 'b');
      time_aveRMST_1PC3(q) = mean(rmst(j, 1:length(t)));
% text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
      figure(f_soft_SE)
      hold on
      q2 = plot(t, strain_energy(j,1:length(t)), 'b');
      time_aveSE_1PC3(q) = mean(strain_energy(j, 1:length(t)));
        time_aveAD_1PC3(q) = mean(average_displacement(j, 1:length(t)));
  q = q+1;
  end 
end 
   
% Plot RMST and SE of 3 kPa DU145 in FBS (non-CSS) 
% Identify Relevant Exps
to_open = text_cats(DU145_1kPa_ind, 9);

% Identify time interval
time_interval = cell2mat(text_cats(DU145_1kPa_ind ,10))/60;

% Identify when the frame for adding 
frame_add = cell2mat(text_cats(DU145_1kPa_ind ,6));

% Identify bad positions
series_omit_array = text_cats(DU145_1kPa_ind ,7);

% Run through each position and plot 
q=1;
for i = 1:length(to_open) 
    % Open file and extract values rxc = pos x time
  data_temp = open(['Mats_So_Far/' to_open{i}]);
  rmst = data_temp.rmst';
  strain_energy = data_temp.strain_energy'*10^12;
   average_displacement = data_temp.average_displacement';
  
  % Identify positions to omit and do so
  if series_omit_array{i} ~=0 
    try 
    omit_indices = str2double(strsplit(series_omit_array{i}, ','));
    catch
    omit_indices = series_omit_array{i}; 
    end 
    rmst(omit_indices,:) = []; 
    strain_energy(omit_indices, :) = [];
    average_displacement(omit_indices,:) = [];
  end     


  [pos_num, time_num] = size(rmst);
  if frame_add(i) ~= 0
      time_num = frame_add(i) - 1;
  end 
  t = ([1:time_num]-1)*time_interval(i); 
  
  for j = 1:pos_num  
      figure(f_soft_rmst)  
      hold on
      p3 = plot(t, rmst(j,1:length(t)), 'g');
      time_aveRMST_1DU145(q) = mean(rmst(j, 1:length(t)));
% text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
      figure(f_soft_SE)
      hold on
      q3 = plot(t, strain_energy(j,1:length(t)), 'g');
      time_aveSE_1DU145(q) = mean(strain_energy(j, 1:length(t)));
       time_aveAD_1DU145(q) = mean(average_displacement(j, 1:length(t)));
  q = q+1;
  end 
  
end 

  
% Plot RMST and SE of 3 kPa DU145 in FBS (non-CSS) 
% Identify Relevant Exps
to_open = text_cats(LNCaP_1kPa_ind, 9);

% Identify time interval
time_interval = cell2mat(text_cats(LNCaP_1kPa_ind ,10))/60;

% Identify when the frame for adding 
frame_add = cell2mat(text_cats(LNCaP_1kPa_ind ,6));

% Identify bad positions
series_omit_array = text_cats(LNCaP_1kPa_ind ,7);

% Run through each position and plot 
q=1;
for i = 1:length(to_open) 
    % Open file and extract values rxc = pos x time
  data_temp = open(['Mats_So_Far/' to_open{i}]);
  rmst = data_temp.rmst';
  strain_energy = data_temp.strain_energy'*10^12;
   average_displacement = data_temp.average_displacement';
  
  % Identify positions to omit and do so
  if series_omit_array{i} ~=0 
    try 
    omit_indices = str2double(strsplit(series_omit_array{i}, ','));
    catch
    omit_indices = series_omit_array{i}; 
    end 
    rmst(omit_indices,:) = []; 
    strain_energy(omit_indices, :) = [];
    average_displacement(omit_indices,:) = [];
  end     

  [pos_num, time_num] = size(rmst);
  if frame_add(i) ~= 0
      time_num = frame_add(i) - 1;
  end 
  t = ([1:time_num]-1)*time_interval(i); 
  
  for j = 1:pos_num  
      figure(f_soft_rmst)  
      hold on
      p4 = plot(t, rmst(j,1:length(t)), 'm');
      time_aveRMST_1LNCaP(q) = mean(rmst(j, 1:length(t)));
% text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
      figure(f_soft_SE)
      hold on
      q4 = plot(t, strain_energy(j,1:length(t)), 'm');
      time_aveSE_1LNCaP(q) = mean(strain_energy(j, 1:length(t)));
        time_aveAD_1LNCaP(q) = mean(average_displacement(j, 1:length(t)));
  q = q+1;
  end 
  
end 


hold off 
figure(f_soft_rmst)
hold on
legend([p1,p4,p3,p2],'22RV1, 1 kPa','LNCaP, 1 kPa', 'DU145, 1 kPa' , 'PC3, 1 kPa')
saveas(gcf, 'soft1_rmst.tiff');

figure(f_soft_SE)
hold on
legend([q1,q4,q3,q2],'22RV1, 1 kPa','LNCaP, 1 kPa', 'DU145, 1 kPa' , 'PC3, 1 kPa')
saveas(gcf, 'soft1_SE.tiff');


%% 50 kPa

f_soft_rmst = figure('Name', 'RMST vs. Time, 22RV1 & PC3 50kPa');
ylabel('RMST (Pa)');
xlabel('Time (Minutes)'); 

f_soft_SE = figure('Name', 'SE vs. Time 22RV1 vs. PC3 50kPa');
ylabel('Strain  Energy (pJ)');
xlabel('Time (Minutes)'); 
  
% Plot RMST and SE of 3 kPa 22RV1

% Identify Relevant Exps
to_open = text_cats(RV1_50kPa_ind, 9);

% Identify tie interval
time_interval = cell2mat(text_cats(RV1_50kPa_ind,10))/60;
frame_add = cell2mat(text_cats(RV1_50kPa_ind,6));

% Identify bad positions
series_omit_array = text_cats(RV1_50kPa_ind,7);

% Run through each position and plot 
q=1;

for i = 1:length(to_open) 
    % Open file and extract values rxc = pos x time
  data_temp = open(['Mats_So_Far/' to_open{i}]);
  rmst = data_temp.rmst';
  strain_energy = data_temp.strain_energy'*10^12; % Adjust to pJ
   average_displacement = data_temp.average_displacement';
  
  % Identify positions to omit and do so
  if series_omit_array{i} ~=0 
    try 
    omit_indices = str2double(strsplit(series_omit_array{i}, ','));
    catch
    omit_indices = series_omit_array{i}; 
    end 
    rmst(omit_indices,:) = []; 
    strain_energy(omit_indices, :) = [];
    average_displacement(omit_indices,:) = [];
  end     


  % Identify the frame at which the treatment was applied 
  [pos_num, time_num] = size(rmst);
  if frame_add(i) ~= 0
      time_num = frame_add(i) - 1;
  end 
  t = ([1:time_num]-1)*time_interval(i); 

  for j = 1:pos_num  
  figure(f_soft_rmst)  
  hold on
  p1 = plot(t, rmst(j,1:length(t)), 'r'); 
  time_aveRMST_50RV1(q) = mean(rmst(j, 1:length(t)));
% text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
  figure(f_soft_SE) 
  hold on
  q1 =  plot(t, strain_energy(j,1:length(t)), 'r'); 
  time_aveSE_50RV1(q) = mean(strain_energy(j, 1:length(t)));
    time_aveAD_50RV1(q) = mean(average_displacement(j, 1:length(t)));
  q = q+1;
  end 
end 

% Plot RMST and SE of 3 kPa PC3
% Identify Relevant Exps
to_open = text_cats(PC3_50kPa_ind , 9);

% Identify time interval
time_interval = cell2mat(text_cats(PC3_50kPa_ind ,10))/60;

% Identify when the frame for adding 
frame_add = cell2mat(text_cats(PC3_50kPa_ind ,6));

% Identify bad positions
series_omit_array = text_cats(PC3_50kPa_ind ,7);

% Run through each position and plot 
q=1;
for i = 1:length(to_open) 
    % Open file and extract values rxc = pos x time
  data_temp = open(['Mats_So_Far/' to_open{i}]);
  rmst = data_temp.rmst';
  strain_energy = data_temp.strain_energy'*10^12;
   average_displacement = data_temp.average_displacement';
  
  % Identify positions to omit and do so
  if series_omit_array{i} ~=0 
    try 
    omit_indices = str2double(strsplit(series_omit_array{i}, ','));
    catch
    omit_indices = series_omit_array{i}; 
    end 
    rmst(omit_indices,:) = []; 
    strain_energy(omit_indices, :) = [];
    average_displacement(omit_indices,:) = [];
  end     


  [pos_num, time_num] = size(rmst);
  if frame_add(i) ~= 0
      time_num = frame_add(i) - 1;
  end 
  t = ([1:time_num]-1)*time_interval(i); 
  
  for j = 1:pos_num  
      figure(f_soft_rmst)  
      hold on
      p2 = plot(t, rmst(j,1:length(t)), 'b');
      time_aveRMST_50PC3(q) = mean(rmst(j, 1:length(t)));
% text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
      figure(f_soft_SE)
      hold on
      q2 = plot(t, strain_energy(j,1:length(t)), 'b');
      time_aveSE_50PC3(q) = mean(strain_energy(j, 1:length(t)));
      time_aveAD_50PC3(q) = mean(average_displacement(j, 1:length(t)));
  q = q+1;
  end 
end 
   
% Plot RMST and SE of 3 kPa DU145 in FBS (non-CSS) 
% Identify Relevant Exps
to_open = text_cats(DU145_50kPa_ind, 9);

% Identify time interval
time_interval = cell2mat(text_cats(DU145_50kPa_ind ,10))/60;

% Identify when the frame for adding 
frame_add = cell2mat(text_cats(DU145_50kPa_ind ,6));

% Identify bad positions
series_omit_array = text_cats(DU145_50kPa_ind ,7);

% Run through each position and plot 
q=1;
for i = 1:length(to_open) 
    % Open file and extract values rxc = pos x time
  data_temp = open(['Mats_So_Far/' to_open{i}]);
  rmst = data_temp.rmst';
  strain_energy = data_temp.strain_energy'*10^12;
   average_displacement = data_temp.average_displacement';
  
  % Identify positions to omit and do so
  if series_omit_array{i} ~=0 
    try 
    omit_indices = str2double(strsplit(series_omit_array{i}, ','));
    catch
    omit_indices = series_omit_array{i}; 
    end 
    rmst(omit_indices,:) = []; 
    strain_energy(omit_indices, :) = [];
    average_displacement(omit_indices,:) = [];
  end     


  [pos_num, time_num] = size(rmst);
  if frame_add(i) ~= 0
      time_num = frame_add(i) - 1;
  end 
  t = ([1:time_num]-1)*time_interval(i); 
  
  for j = 1:pos_num  
      figure(f_soft_rmst)  
      hold on
      p3 = plot(t, rmst(j,1:length(t)), 'g');
      time_aveRMST_50DU145(q) = mean(rmst(j, 1:length(t)));
% text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
      figure(f_soft_SE)
      hold on
      q3 = plot(t, strain_energy(j,1:length(t)), 'g');
      time_aveSE_50DU145(q) = mean(strain_energy(j, 1:length(t)));
        time_aveAD_50DU145(q) = mean(average_displacement(j, 1:length(t)));
  q = q+1;
  end 
  
end 

  
% Plot RMST and SE of 3 kPa DU145 in FBS (non-CSS) 
% Identify Relevant Exps
to_open = text_cats(LNCaP_50kPa_ind, 9);

% Identify time interval
time_interval = cell2mat(text_cats(LNCaP_50kPa_ind ,10))/60;

% Identify when the frame for adding 
frame_add = cell2mat(text_cats(LNCaP_50kPa_ind ,6));

% Identify bad positions
series_omit_array = text_cats(LNCaP_50kPa_ind ,7);

% Run through each position and plot 
q=1;
for i = 1:length(to_open) 
    % Open file and extract values rxc = pos x time
  data_temp = open(['Mats_So_Far/' to_open{i}]);
  rmst = data_temp.rmst';
  strain_energy = data_temp.strain_energy'*10^12;
   average_displacement = data_temp.average_displacement';
  
  % Identify positions to omit and do so
  if series_omit_array{i} ~=0 
    try 
    omit_indices = str2double(strsplit(series_omit_array{i}, ','));
    catch
    omit_indices = series_omit_array{i}; 
    end 
    rmst(omit_indices,:) = []; 
    strain_energy(omit_indices, :) = [];
    average_displacement(omit_indices,:) = [];
  end     

  [pos_num, time_num] = size(rmst);
  if frame_add(i) ~= 0
      time_num = frame_add(i) - 1;
  end 
  t = ([1:time_num]-1)*time_interval(i); 
  
  for j = 1:pos_num  
      figure(f_soft_rmst)  
      hold on
      p4 = plot(t, rmst(j,1:length(t)), 'm');
      time_aveRMST_50LNCaP(q) = mean(rmst(j, 1:length(t)));
% text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
      figure(f_soft_SE)
      hold on
      q4 = plot(t, strain_energy(j,1:length(t)), 'm');
      time_aveSE_50LNCaP(q) = mean(strain_energy(j, 1:length(t)));
      time_aveAD_50LNCaP(q) = mean(average_displacement(j, 1:length(t)));
  q = q+1;
  end 
  
end 


hold off 
figure(f_soft_rmst)
hold on
legend([p1,p4,p3,p2],'22RV1, 50 kPa','LNCaP, 50 kPa', 'DU145, 50 kPa' , 'PC3, 50 kPa')
saveas(gcf, 'soft1_rmst.tiff');

figure(f_soft_SE)
hold on
legend([q1,q4,q3,q2],'22RV1, 50 kPa','LNCaP, 50 kPa', 'DU145, 50 kPa' , 'PC3, 50 kPa')
saveas(gcf, 'soft1_SE.tiff');


%% RMST and Strain Energy Boxplots Soft
% f_box_soft_rmst = figure('Name', 'RMST boxplot 22RV1 vs. PC3 3kPa');
% ylabel('RMST (Pa)');
% 
% f_box_soft_SE = figure('Name', 'SE boxplot 22RV1 vs. PC3 3kPa');
% ylabel('Strain  Energy (pJ)');
% 
% figure(f_box_soft_rmst)
% hold on 
% boxplot([time_aveRMST_3RV1,time_aveRMST_3LNCaP ,time_aveRMST_3DU145, time_aveRMST_3PC3], [1*ones((length(time_aveRMST_3RV1)),1); 2*ones((length(time_aveRMST_3LNCaP)),1); ...
%           3*ones((length(time_aveRMST_3DU145)),1); 4*ones((length(time_aveRMST_3PC3)),1)], 'Positions',[1*ones((length(time_aveRMST_3RV1)),1); 2*ones((length(time_aveRMST_3LNCaP)),1); ...
%           3*ones((length(time_aveRMST_3DU145)),1); 4*ones((length(time_aveRMST_3PC3)),1)],  'Labels', {'22RV1 3 kPa', 'LNCaP', 'DU145', 'PC3 3 kPa'});
%       
% RMST_scatter_ind3RV1 = 1*ones((length(time_aveRMST_3RV1)),1);
% f1 = scatter(RMST_scatter_ind3RV1, time_aveRMST_3RV1, 'r');
% 
% RMST_scatter_ind3LNCaP = 2*ones((length(time_aveRMST_3LNCaP)),1);
% f1 = scatter(RMST_scatter_ind3LNCaP, time_aveRMST_3LNCaP, 'r');
% 
% RMST_scatter_ind3DU145 = 3*ones((length(time_aveRMST_3DU145)),1);
% f1 = scatter(RMST_scatter_ind3DU145, time_aveRMST_3DU145, 'r');
% 
% RMST_scatter_ind3PC3 = 4*ones((length(time_aveRMST_3PC3)),1);
% f1 = scatter(RMST_scatter_ind3PC3, time_aveRMST_3PC3, 'b');
%       
% % saveas('f_box_soft_rmst.tiff');
% 
% figure(f_box_soft_SE)
% hold on 
% boxplot([time_aveSE_3RV1,time_aveSE_3LNCaP ,time_aveSE_3DU145, time_aveSE_3PC3], [1*ones((length(time_aveSE_3RV1)),1); 2*ones((length(time_aveSE_3LNCaP)),1); ...
%           3*ones((length(time_aveSE_3DU145)),1); 4*ones((length(time_aveSE_3PC3)),1)], 'Positions',[1*ones((length(time_aveSE_3RV1)),1); 2*ones((length(time_aveSE_3LNCaP)),1); ...
%           3*ones((length(time_aveSE_3DU145)),1); 4*ones((length(time_aveSE_3PC3)),1)],  'Labels', {'22RV1 3 kPa', 'LNCaP', 'DU145', 'PC3 3 kPa'});
%       
% SE_scatter_ind3RV1 = 1*ones((length(time_aveSE_3RV1)),1);
% f1 = scatter(SE_scatter_ind3RV1, time_aveSE_3RV1, 'r');
% 
% SE_scatter_ind3LNCaP = 2*ones((length(time_aveSE_3LNCaP)),1);
% f1 = scatter(SE_scatter_ind3LNCaP, time_aveSE_3LNCaP, 'r');
% 
% SE_scatter_ind3DU145 = 3*ones((length(time_aveSE_3DU145)),1);
% f1 = scatter(SE_scatter_ind3DU145, time_aveSE_3DU145, 'r');
% 
% SE_scatter_ind3PC3 = 4*ones((length(time_aveRMST_3PC3)),1);
% f1 = scatter(SE_scatter_ind3PC3, time_aveSE_3PC3, 'b');

% saveas(gcf, 'f_box_soft_SE.tiff');


% %%
% % RMST boxplots 
% 
% f_box_soft_rmst = figure('Name', 'RMST boxplot, 22RV1 vs. PC3 3 kPa');
% ylabel('RMST (Pa)');
% % xlabel('Time (Minutes)'); 
% 
% f_box_soft_SE = figure('Name', 'SE boxplot, 22RV1 vs. PC3 3 kPa');
% ylabel('Strain  Energy (pJ)');
% % xlabel('Time (Minutes)'); 
% 
% 
% figure(f_box_soft_rmst)
% hold on 
% boxplot([time_aveRMST_3RV1, time_aveRMST_3PC3], [1*ones((length(time_aveRMST_3RV1)),1);...
%           2*ones((length(time_aveRMST_3PC3)),1)], 'Positions',[1*ones((length(time_aveRMST_3RV1)),1);...
%           2*ones((length(time_aveRMST_3PC3)),1)],  'Labels', {'22RV1, 3 kPa', 'PC3, 3 kPa'});
%       
% RMST_scatter_ind3RV1 = 1*ones((length(time_aveRMST_3RV1)),1);
% f1 = scatter(RMST_scatter_ind3RV1, time_aveRMST_3RV1, 'r');
%   
% RMST_scatter_ind3PC3 = 2*ones((length(time_aveRMST_3PC3)),1);
% f1 = scatter(RMST_scatter_ind3PC3, time_aveRMST_3PC3, 'b');
% 
% saveas(gcf, 'f_box_soft_rmst.tiff')
%  
% figure(f_box_soft_SE)
% hold on 
% boxplot([time_aveSE_3RV1, time_aveSE_3PC3], [1*ones((length(time_aveSE_3RV1)),1);...
%           2*ones((length(time_aveSE_3PC3)),1)], 'Positions',[1*ones((length(time_aveSE_3RV1)),1);...
%           2*ones((length(time_aveSE_3PC3)),1)],  'Labels', {'22RV1, 3 kPa', 'PC3, 3 kPa'});
%       
% SE_scatter_ind3RV1 = 1*ones((length(time_aveSE_3RV1)),1);
% f1 = scatter(SE_scatter_ind3RV1, time_aveSE_3RV1, 'r');
%   
% SE_scatter_ind3PC3 = 2*ones((length(time_aveSE_3PC3)),1);
% f1 = scatter(SE_scatter_ind3PC3, time_aveSE_3PC3, 'b');
% 
% saveas(gcf, 'f_box_soft_SE.tiff')
 
%%  Stiffness Sensitivity Summary 
f_box_rmst_softness_compare = figure('Name', 'FBS_PC3/22RV1 RMST Stiffness Compare');
hold on
ylabel('RMST (Pa)');

boxplot([time_aveRMST_1RV1, time_aveRMST_3RV1, time_aveRMST_12RV1,  ...
         time_aveRMST_1LNCaP, time_aveRMST_3LNCaP, time_aveRMST_12LNCaP,...
         time_aveRMST_1DU145, time_aveRMST_3DU145, time_aveRMST_12DU145,...
         time_aveRMST_1PC3, time_aveRMST_3PC3, time_aveRMST_12PC3],...
         [0.4*ones((length(time_aveRMST_1RV1)),1); 0.8*ones((length(time_aveRMST_3RV1)),1);... 
          1.2*ones((length(time_aveRMST_12RV1)),1);  ...
          2.4*ones((length(time_aveRMST_1LNCaP)),1); 2.8*ones((length(time_aveRMST_3LNCaP)),1); ...
          3.2*ones((length(time_aveRMST_12LNCaP)),1);  ...
          4.4*ones((length(time_aveRMST_1DU145)),1); 4.8*ones((length(time_aveRMST_3DU145)),1);...
          5.2*ones((length(time_aveRMST_12DU145)),1); ...
          6.4*ones((length( time_aveRMST_1PC3)),1); 6.8*ones((length( time_aveRMST_3PC3)),1);......
          7.2*ones((length(time_aveRMST_12PC3)),1)],...
          'Positions', ...
          [0.4*ones((length(time_aveRMST_1RV1)),1); 0.8*ones((length(time_aveRMST_3RV1)),1);... 
          1.2*ones((length(time_aveRMST_12RV1)),1);  ...
          2.4*ones((length(time_aveRMST_1LNCaP)),1); 2.8*ones((length(time_aveRMST_3LNCaP)),1); ...
          3.2*ones((length(time_aveRMST_12LNCaP)),1); ...
          4.4*ones((length(time_aveRMST_1DU145)),1); 4.8*ones((length(time_aveRMST_3DU145)),1);...
          5.2*ones((length(time_aveRMST_12DU145)),1);...
          6.4*ones((length( time_aveRMST_1PC3)),1); 6.8*ones((length( time_aveRMST_3PC3)),1);...
          7.2*ones((length(time_aveRMST_12PC3)),1)], ...
          'Labels', ...
          {'22RV1, 1 kPa', '22RV1, 3 kPa', '22RV1, 12 kPa',...
          'LNCaP, 1 kPa', 'LNCaP, 3 kPa', 'LNCaP, 12 kPa',...
          'DU145, 1 kPa', 'DU145, 3 kPa', 'DU145, 12 kPa',...
          'PC3, 1 kPa', 'PC3, 3 kPa', 'PC3, 12 kPa'}, 'widths', 0.3)
set(gca, 'FontSize', 15);
    
delete(findobj(gca,'type','text'))  % wipe out the boxplot labels that are problematic
set(gca,'xtick',0.8:2:6.8,'xticklabel',{ '22RV1' 'LNCaP' 'DU145' 'PC3'})  % and write the labels wanted
% xlim([0.5,4.5])
RMST_scatter_ind1RV1 = 0.4*ones((length(time_aveRMST_1RV1)),1);
f1 = scatter(RMST_scatter_ind1RV1, time_aveRMST_1RV1, 'ro');
f2 = scatter(0.4, mean(time_aveRMST_1RV1), 'ko', 'LineWidth', 2);

RMST_scatter_ind3RV1 = 0.8*ones((length(time_aveRMST_3RV1)),1);
f1 = scatter(RMST_scatter_ind3RV1, time_aveRMST_3RV1, 'rd');
f2 = scatter(0.8, mean(time_aveRMST_3RV1), 'kd', 'LineWidth', 2);

RMST_scatter_ind12RV1 = 1.2*ones((length(time_aveRMST_12RV1)),1);
f1 = scatter(RMST_scatter_ind12RV1, time_aveRMST_12RV1, 'r*');
f2 = scatter(1.2, mean(time_aveRMST_12RV1), 'k*', 'LineWidth', 1);
% 
% RMST_scatter_ind50RV1 =1.6*ones((length(time_aveRMST_50RV1)),1);
% f1 = scatter(RMST_scatter_ind50RV1, time_aveRMST_50RV1, 'r^');
% f2 = scatter(1.6, mean(time_aveRMST_50RV1), 'k^', 'LineWidth', 2);
%%%%%%%%%%%%%%%%%%%

RMST_scatter_ind1LNCaP = 2.4*ones((length(time_aveRMST_1LNCaP)),1);
f1 = scatter(RMST_scatter_ind1LNCaP, time_aveRMST_1LNCaP, 'mo');
f2 = scatter(2.4, mean(time_aveRMST_1LNCaP), 'kd', 'LineWidth', 1);


RMST_scatter_ind3LNCaP = 2.8*ones((length(time_aveRMST_3LNCaP)),1);
f1 = scatter(RMST_scatter_ind3LNCaP, time_aveRMST_3LNCaP, 'md');
f2 = scatter(2.8, mean(time_aveRMST_3LNCaP), 'kd', 'LineWidth', 1);


RMST_scatter_ind12LNCaP = 3.2*ones((length(time_aveRMST_12LNCaP)),1);
f1 = scatter(RMST_scatter_ind12LNCaP, time_aveRMST_12LNCaP, 'm*');
f2 = scatter(3.2, mean(time_aveRMST_12LNCaP), 'k*', 'LineWidth', 1);
% 
% RMST_scatter_ind50LNCaP = 3.6*ones((length(time_aveRMST_50LNCaP)),1);
% f1 = scatter(RMST_scatter_ind50LNCaP, time_aveRMST_50LNCaP, 'm^');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RMST_scatter_ind1DU145 = 4.4*ones((length(time_aveRMST_1DU145)),1);
f1 = scatter(RMST_scatter_ind1DU145, time_aveRMST_1DU145, 'go');
f2 = scatter(4.4, mean(time_aveRMST_1DU145), 'ko', 'LineWidth', 1);

RMST_scatter_ind3DU145 = 4.8*ones((length(time_aveRMST_3DU145)),1);
f1 = scatter(RMST_scatter_ind3DU145, time_aveRMST_3DU145, 'gd');
f2 = scatter(4.8, mean(time_aveRMST_3DU145), 'kd', 'LineWidth', 1);

RMST_scatter_ind12DU145 = 5.2*ones((length(time_aveRMST_12DU145)),1);
f1 = scatter(RMST_scatter_ind12DU145, time_aveRMST_12DU145, 'g*');
f2 = scatter(5.2, mean(time_aveRMST_12DU145), 'k*', 'LineWidth', 1);
% 

% RMST_scatter_ind50DU145 = 5.6*ones((length(time_aveRMST_50DU145)),1);
% f1 = scatter(RMST_scatter_ind50DU145, time_aveRMST_50DU145, 'g^');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RMST_scatter_ind1PC3 = 6.4*ones((length(time_aveRMST_1PC3)),1);
f1 = scatter(RMST_scatter_ind1PC3, time_aveRMST_1PC3, 'bo');
f2 = scatter(6.4, mean(time_aveRMST_1PC3), 'ko', 'LineWidth', 1);


RMST_scatter_ind3PC3 = 6.8*ones((length(time_aveRMST_3PC3)),1);
f1 = scatter(RMST_scatter_ind3PC3, time_aveRMST_3PC3, 'bd');
f2 = scatter(6.8, mean(time_aveRMST_3PC3), 'kd', 'LineWidth', 1);


RMST_scatter_ind12PC3 = 7.2*ones((length(time_aveRMST_12PC3)),1);
f1 = scatter(RMST_scatter_ind12PC3, time_aveRMST_12PC3, 'b*');
f2 = scatter(7.2, mean(time_aveRMST_12PC3), 'k*', 'LineWidth', 1);

% 
% RMST_scatter_ind50PC3 = 7.6*ones((length(time_aveRMST_50PC3)),1);
% f1 = scatter(RMST_scatter_ind50PC3, time_aveRMST_50PC3, 'b^');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlim([-0.2, 7.8]);


mlabel1 = scatter(10,1, 'ko'); 
mlabel2 = scatter(10,1, 'kd'); 
mlabel3 = scatter(10,1, 'k*'); 
mlabel4 = scatter(10,1, 'k^'); 


legend([mlabel1, mlabel2, mlabel3, mlabel4], '1 kPa','3 kPa', '12 kPa', '50 kPa', 'Location', 'northwest')
set(gca, 'FontSize', 15, 'LineWidth', 1);
% ylim([0,16])
ylim([0,400])


saveas(gcf, 'RMST_Summary', 'tiff')

%% Strain Energy 

f_box_rmst_softness_compare = figure('Name', 'FBS_PC3/22RV1 SE Stiffness Compare');
hold on
ylabel('Strain Energy (pJ)');
boxplot([time_aveSE_1RV1, time_aveSE_3RV1, time_aveSE_12RV1,  ...
         time_aveSE_1LNCaP, time_aveSE_3LNCaP, time_aveSE_12LNCaP,...
         time_aveSE_1DU145, time_aveSE_3DU145, time_aveSE_12DU145,...
         time_aveSE_1PC3, time_aveSE_3PC3, time_aveSE_12PC3],...
         [0.4*ones((length(time_aveSE_1RV1)),1); 0.8*ones((length(time_aveSE_3RV1)),1);... 
          1.2*ones((length(time_aveSE_12RV1)),1);  ...
          2.4*ones((length(time_aveSE_1LNCaP)),1); 2.8*ones((length(time_aveSE_3LNCaP)),1); ...
          3.2*ones((length(time_aveSE_12LNCaP)),1);  ...
          4.4*ones((length(time_aveSE_1DU145)),1); 4.8*ones((length(time_aveSE_3DU145)),1);...
          5.2*ones((length(time_aveSE_12DU145)),1); ...
          6.4*ones((length( time_aveSE_1PC3)),1); 6.8*ones((length( time_aveSE_3PC3)),1);......
          7.2*ones((length(time_aveSE_12PC3)),1)],...
          'Positions', ...
          [0.4*ones((length(time_aveSE_1RV1)),1); 0.8*ones((length(time_aveSE_3RV1)),1);... 
          1.2*ones((length(time_aveSE_12RV1)),1);  ...
          2.4*ones((length(time_aveSE_1LNCaP)),1); 2.8*ones((length(time_aveSE_3LNCaP)),1); ...
          3.2*ones((length(time_aveSE_12LNCaP)),1); ...
          4.4*ones((length(time_aveSE_1DU145)),1); 4.8*ones((length(time_aveSE_3DU145)),1);...
          5.2*ones((length(time_aveSE_12DU145)),1);...
          6.4*ones((length( time_aveSE_1PC3)),1); 6.8*ones((length( time_aveSE_3PC3)),1);...
          7.2*ones((length(time_aveSE_12PC3)),1)], ...
          'Labels', ...
          {'22RV1, 1 kPa', '22RV1, 3 kPa', '22RV1, 12 kPa',...
          'LNCaP, 1 kPa', 'LNCaP, 3 kPa', 'LNCaP, 12 kPa',...
          'DU145, 1 kPa', 'DU145, 3 kPa', 'DU145, 12 kPa',...
          'PC3, 1 kPa', 'PC3, 3 kPa', 'PC3, 12 kPa'}, 'widths', 0.3)
set(gca, 'FontSize', 15);
    
delete(findobj(gca,'type','text'))  % wipe out the boxplot labels that are problematic
set(gca,'xtick',0.8:2:6.8,'xticklabel',{ '22RV1' 'LNCaP' 'DU145' 'PC3'})  % and write the labels wanted
% xlim([0.5,4.5])
SE_scatter_ind1RV1 = 0.4*ones((length(time_aveSE_1RV1)),1);
f1 = scatter(SE_scatter_ind1RV1, time_aveSE_1RV1, 'ro');
f2 = scatter(0.4, mean(time_aveSE_1RV1), 'ko', 'LineWidth', 2);

SE_scatter_ind3RV1 = 0.8*ones((length(time_aveSE_3RV1)),1);
f1 = scatter(SE_scatter_ind3RV1, time_aveSE_3RV1, 'rd');
f2 = scatter(0.8, mean(time_aveSE_3RV1), 'kd', 'LineWidth', 2);

SE_scatter_ind12RV1 = 1.2*ones((length(time_aveSE_12RV1)),1);
f1 = scatter(SE_scatter_ind12RV1, time_aveSE_12RV1, 'r*');
f2 = scatter(1.2, mean(time_aveSE_12RV1), 'k*', 'LineWidth', 1);
% 
% SE_scatter_ind50RV1 =1.6*ones((length(time_aveSE_50RV1)),1);
% f1 = scatter(SE_scatter_ind50RV1, time_aveSE_50RV1, 'r^');
% f2 = scatter(1.6, mean(time_aveSE_50RV1), 'k^', 'LineWidth', 2);
%%%%%%%%%%%%%%%%%%%

SE_scatter_ind1LNCaP = 2.4*ones((length(time_aveSE_1LNCaP)),1);
f1 = scatter(SE_scatter_ind1LNCaP, time_aveSE_1LNCaP, 'mo');
f2 = scatter(2.4, mean(time_aveSE_1LNCaP), 'kd', 'LineWidth', 1);


SE_scatter_ind3LNCaP = 2.8*ones((length(time_aveSE_3LNCaP)),1);
f1 = scatter(SE_scatter_ind3LNCaP, time_aveSE_3LNCaP, 'md');
f2 = scatter(2.8, mean(time_aveSE_3LNCaP), 'kd', 'LineWidth', 1);


SE_scatter_ind12LNCaP = 3.2*ones((length(time_aveSE_12LNCaP)),1);
f1 = scatter(SE_scatter_ind12LNCaP, time_aveSE_12LNCaP, 'm*');
f2 = scatter(3.2, mean(time_aveSE_12LNCaP), 'k*', 'LineWidth', 1);
% 
% SE_scatter_ind50LNCaP = 3.6*ones((length(time_aveSE_50LNCaP)),1);
% f1 = scatter(SE_scatter_ind50LNCaP, time_aveSE_50LNCaP, 'm^');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SE_scatter_ind1DU145 = 4.4*ones((length(time_aveSE_1DU145)),1);
f1 = scatter(SE_scatter_ind1DU145, time_aveSE_1DU145, 'go');
f2 = scatter(4.4, mean(time_aveSE_1DU145), 'ko', 'LineWidth', 1);

SE_scatter_ind3DU145 = 4.8*ones((length(time_aveSE_3DU145)),1);
f1 = scatter(SE_scatter_ind3DU145, time_aveSE_3DU145, 'gd');
f2 = scatter(4.8, mean(time_aveSE_3DU145), 'kd', 'LineWidth', 1);

SE_scatter_ind12DU145 = 5.2*ones((length(time_aveSE_12DU145)),1);
f1 = scatter(SE_scatter_ind12DU145, time_aveSE_12DU145, 'g*');
f2 = scatter(5.2, mean(time_aveSE_12DU145), 'k*', 'LineWidth', 1);
% 

% SE_scatter_ind50DU145 = 5.6*ones((length(time_aveSE_50DU145)),1);
% f1 = scatter(SE_scatter_ind50DU145, time_aveSE_50DU145, 'g^');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SE_scatter_ind1PC3 = 6.4*ones((length(time_aveSE_1PC3)),1);
f1 = scatter(SE_scatter_ind1PC3, time_aveSE_1PC3, 'bo');
f2 = scatter(6.4, mean(time_aveSE_1PC3), 'ko', 'LineWidth', 1);


SE_scatter_ind3PC3 = 6.8*ones((length(time_aveSE_3PC3)),1);
f1 = scatter(SE_scatter_ind3PC3, time_aveSE_3PC3, 'bd');
f2 = scatter(6.8, mean(time_aveSE_3PC3), 'kd', 'LineWidth', 1);


SE_scatter_ind12PC3 = 7.2*ones((length(time_aveSE_12PC3)),1);
f1 = scatter(SE_scatter_ind12PC3, time_aveSE_12PC3, 'b*');
f2 = scatter(7.2, mean(time_aveSE_12PC3), 'k*', 'LineWidth', 1);

% 
% RMST_scatter_ind50PC3 = 7.6*ones((length(time_aveRMST_50PC3)),1);
% f1 = scatter(RMST_scatter_ind50PC3, time_aveRMST_50PC3, 'b^');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlim([-0.2, 7.8]);


mlabel1 = scatter(10,1, 'ko'); 
mlabel2 = scatter(10,1, 'kd'); 
mlabel3 = scatter(10,1, 'k*'); 
mlabel4 = scatter(10,1, 'k^'); 


legend([mlabel1, mlabel2, mlabel3, mlabel4], '1 kPa','3 kPa', '12 kPa', '50 kPa', 'Location', 'northwest')
set(gca, 'FontSize', 15, 'LineWidth', 1);
% ylim([0,16])
ylim([0,10])


saveas(gcf, 'SE_Summary', 'tiff')

%%

f_box_rmst_softness_compare = figure('Name', 'FBS_PC3/22RV1 Displacement Compare');
hold on
ylabel('Average Displacement (\mum)');

boxplot([time_aveAD_1RV1, time_aveAD_3RV1, time_aveAD_12RV1,  ...
         time_aveAD_1LNCaP, time_aveAD_3LNCaP, time_aveAD_12LNCaP,...
         time_aveAD_1DU145, time_aveAD_3DU145, time_aveAD_12DU145,...
         time_aveAD_1PC3, time_aveAD_3PC3, time_aveAD_12PC3],...
         [0.4*ones((length(time_aveAD_1RV1)),1); 0.8*ones((length(time_aveAD_3RV1)),1);... 
          1.2*ones((length(time_aveAD_12RV1)),1);  ...
          2.4*ones((length(time_aveAD_1LNCaP)),1); 2.8*ones((length(time_aveAD_3LNCaP)),1); ...
          3.2*ones((length(time_aveAD_12LNCaP)),1);  ...
          4.4*ones((length(time_aveAD_1DU145)),1); 4.8*ones((length(time_aveAD_3DU145)),1);...
          5.2*ones((length(time_aveAD_12DU145)),1); ...
          6.4*ones((length( time_aveAD_1PC3)),1); 6.8*ones((length( time_aveAD_3PC3)),1);......
          7.2*ones((length(time_aveAD_12PC3)),1)],...
          'Positions', ...
          [0.4*ones((length(time_aveAD_1RV1)),1); 0.8*ones((length(time_aveAD_3RV1)),1);... 
          1.2*ones((length(time_aveAD_12RV1)),1);  ...
          2.4*ones((length(time_aveAD_1LNCaP)),1); 2.8*ones((length(time_aveAD_3LNCaP)),1); ...
          3.2*ones((length(time_aveAD_12LNCaP)),1); ...
          4.4*ones((length(time_aveAD_1DU145)),1); 4.8*ones((length(time_aveAD_3DU145)),1);...
          5.2*ones((length(time_aveAD_12DU145)),1);...
          6.4*ones((length( time_aveAD_1PC3)),1); 6.8*ones((length( time_aveAD_3PC3)),1);...
          7.2*ones((length(time_aveAD_12PC3)),1)], ...
          'Labels', ...
          {'22RV1, 1 kPa', '22RV1, 3 kPa', '22RV1, 12 kPa',...
          'LNCaP, 1 kPa', 'LNCaP, 3 kPa', 'LNCaP, 12 kPa',...
          'DU145, 1 kPa', 'DU145, 3 kPa', 'DU145, 12 kPa',...
          'PC3, 1 kPa', 'PC3, 3 kPa', 'PC3, 12 kPa'}, 'widths', 0.3)
set(gca, 'FontSize', 15);
    
delete(findobj(gca,'type','text'))  % wipe out the boxplot labels that are problematic
set(gca,'xtick',0.8:2:6.8,'xticklabel',{ '22RV1' 'LNCaP' 'DU145' 'PC3'})  % and write the labels wanted
% xlim([0.5,4.5])
AD_scatter_ind1RV1 = 0.4*ones((length(time_aveAD_1RV1)),1);
f1 = scatter(AD_scatter_ind1RV1, time_aveAD_1RV1, 'ro');
f2 = scatter(0.4, mean(time_aveAD_1RV1), 'ko', 'LineWidth', 2);

AD_scatter_ind3RV1 = 0.8*ones((length(time_aveAD_3RV1)),1);
f1 = scatter(AD_scatter_ind3RV1, time_aveAD_3RV1, 'rd');
f2 = scatter(0.8, mean(time_aveAD_3RV1), 'kd', 'LineWidth', 2);

AD_scatter_ind12RV1 = 1.2*ones((length(time_aveAD_12RV1)),1);
f1 = scatter(AD_scatter_ind12RV1, time_aveAD_12RV1, 'r*');
f2 = scatter(1.2, mean(time_aveAD_12RV1), 'k*', 'LineWidth', 1);
% 
% AD_scatter_ind50RV1 =1.6*ones((length(time_aveAD_50RV1)),1);
% f1 = scatter(AD_scatter_ind50RV1, time_aveAD_50RV1, 'r^');
% f2 = scatter(1.6, mean(time_aveAD_50RV1), 'k^', 'LineWidth', 2);
%%%%%%%%%%%%%%%%%%%

AD_scatter_ind1LNCaP = 2.4*ones((length(time_aveAD_1LNCaP)),1);
f1 = scatter(AD_scatter_ind1LNCaP, time_aveAD_1LNCaP, 'mo');
f2 = scatter(2.4, mean(time_aveAD_1LNCaP), 'kd', 'LineWidth', 1);


AD_scatter_ind3LNCaP = 2.8*ones((length(time_aveAD_3LNCaP)),1);
f1 = scatter(AD_scatter_ind3LNCaP, time_aveAD_3LNCaP, 'md');
f2 = scatter(2.8, mean(time_aveAD_3LNCaP), 'kd', 'LineWidth', 1);


AD_scatter_ind12LNCaP = 3.2*ones((length(time_aveAD_12LNCaP)),1);
f1 = scatter(AD_scatter_ind12LNCaP, time_aveAD_12LNCaP, 'm*');
f2 = scatter(3.2, mean(time_aveAD_12LNCaP), 'k*', 'LineWidth', 1);
% 
% AD_scatter_ind50LNCaP = 3.6*ones((length(time_aveAD_50LNCaP)),1);
% f1 = scatter(AD_scatter_ind50LNCaP, time_aveAD_50LNCaP, 'm^');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AD_scatter_ind1DU145 = 4.4*ones((length(time_aveAD_1DU145)),1);
f1 = scatter(AD_scatter_ind1DU145, time_aveAD_1DU145, 'go');
f2 = scatter(4.4, mean(time_aveAD_1DU145), 'ko', 'LineWidth', 1);

AD_scatter_ind3DU145 = 4.8*ones((length(time_aveAD_3DU145)),1);
f1 = scatter(AD_scatter_ind3DU145, time_aveAD_3DU145, 'gd');
f2 = scatter(4.8, mean(time_aveAD_3DU145), 'kd', 'LineWidth', 1);

AD_scatter_ind12DU145 = 5.2*ones((length(time_aveAD_12DU145)),1);
f1 = scatter(AD_scatter_ind12DU145, time_aveAD_12DU145, 'g*');
f2 = scatter(5.2, mean(time_aveAD_12DU145), 'k*', 'LineWidth', 1);
% 

% AD_scatter_ind50DU145 = 5.6*ones((length(time_aveAD_50DU145)),1);
% f1 = scatter(AD_scatter_ind50DU145, time_aveAD_50DU145, 'g^');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AD_scatter_ind1PC3 = 6.4*ones((length(time_aveAD_1PC3)),1);
f1 = scatter(AD_scatter_ind1PC3, time_aveAD_1PC3, 'bo');
f2 = scatter(6.4, mean(time_aveAD_1PC3), 'ko', 'LineWidth', 1);


AD_scatter_ind3PC3 = 6.8*ones((length(time_aveAD_3PC3)),1);
f1 = scatter(AD_scatter_ind3PC3, time_aveAD_3PC3, 'bd');
f2 = scatter(6.8, mean(time_aveAD_3PC3), 'kd', 'LineWidth', 1);


AD_scatter_ind12PC3 = 7.2*ones((length(time_aveAD_12PC3)),1);
f1 = scatter(AD_scatter_ind12PC3, time_aveAD_12PC3, 'b*');
f2 = scatter(7.2, mean(time_aveAD_12PC3), 'k*', 'LineWidth', 1);

% 
% RMST_scatter_ind50PC3 = 7.6*ones((length(time_aveRMST_50PC3)),1);
% f1 = scatter(RMST_scatter_ind50PC3, time_aveRMST_50PC3, 'b^');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlim([-0.2, 7.8]);


mlabel1 = scatter(10,1, 'ko'); 
mlabel2 = scatter(10,1, 'kd'); 
mlabel3 = scatter(10,1, 'k*'); 
% mlabel4 = scatter(10,1, 'k^'); 


legend([mlabel1, mlabel2, mlabel3], '1 kPa','3 kPa', '12 kPa', 'Location', 'northwest') % '50 kPa',
set(gca, 'FontSize', 15, 'LineWidth', 1);

mlabel1 = scatter(10,1, 'ko'); 
mlabel2 = scatter(10,1, 'kd'); 
mlabel3 = scatter(10,1, 'k*'); 
% mlabel4 = scatter(10,1, 'k^'); 


legend([mlabel1, mlabel2, mlabel3], '1 kPa','3 kPa', '12 kPa', 'Location', 'northeast')
set(gca, 'FontSize', 15, 'LineWidth', 1);
% ylim([0,16])
ylim([0,8])
saveas(gcf, 'Average Displacement_Summary', 'tiff')

%%


f_box_rmst_softness_compare = figure('Name', 'Cluster SE');
hold on
ylabel('Strain Energy (pJ)');

boxplot([time_aveSE_1RV1,  time_aveSE_1LNCaP, time_aveSE_1DU145,time_aveSE_1PC3, ...
         time_aveSE_3RV1, time_aveSE_3LNCaP, time_aveSE_3DU145,  time_aveSE_3PC3,...
         time_aveSE_12RV1, time_aveSE_12LNCaP, time_aveSE_12DU145,  time_aveSE_12PC3],... time_aveSE_50RV1, time_aveSE_50LNCaP, time_aveSE_50DU145, time_aveSE_50PC3
         [0.4*ones((length(time_aveSE_1RV1)),1); 0.8*ones((length(time_aveSE_1LNCaP)),1);... 
          1.2*ones((length(time_aveSE_1DU145)),1); 1.6*ones((length( time_aveSE_1PC3)),1); ...
          2.4*ones((length(time_aveSE_3RV1)),1); 2.8*ones((length( time_aveSE_3LNCaP)),1); ...
          3.2*ones((length(time_aveSE_3DU145)),1); 3.6*ones((length( time_aveSE_3PC3)),1); ...
          4.4*ones((length( time_aveSE_12RV1)),1); 4.8*ones((length( time_aveSE_12LNCaP)),1);...
          5.2*ones((length(time_aveSE_12DU145)),1); 5.6*ones((length(time_aveSE_12PC3)),1);...
%           6.4*ones((length( time_aveSE_50RV1)),1); 6.8*ones((length(time_aveSE_50LNCaP)),1);...
%           7.2*ones((length(time_aveSE_50DU145)),1); 7.6*ones((length(time_aveSE_50PC3)),1)
            ],...
          'Positions', ...
           [0.4*ones((length(time_aveSE_1RV1)),1); 0.8*ones((length(time_aveSE_1LNCaP)),1);... 
          1.2*ones((length(time_aveSE_1DU145)),1); 1.6*ones((length(time_aveSE_1PC3)),1); ...
          2.4*ones((length(time_aveSE_3RV1)),1); 2.8*ones((length(time_aveSE_3LNCaP)),1); ...
          3.2*ones((length(time_aveSE_3DU145)),1); 3.6*ones((length(time_aveSE_3PC3)),1); ...
          4.4*ones((length( time_aveSE_12RV1)),1); 4.8*ones((length(time_aveSE_12LNCaP)),1);...
          5.2*ones((length(time_aveSE_12DU145)),1); 5.6*ones((length(time_aveSE_12PC3)),1);...
%           6.4*ones((length(time_aveSE_50RV1)),1); 6.8*ones((length(time_aveSE_50LNCaP)),1);......
%           7.2*ones((length(time_aveSE_50DU145)),1); 7.6*ones((length(time_aveSE_50PC3)),1)
           ],...
          'Labels', ...
          {'22RV1, 1 kPa', '22RV1, 3 kPa', '22RV1, 12 kPa'... , '22RV1, 50 kPa',
          'LNCaP, 1 kPa', 'LNCaP, 3 kPa', 'LNCaP, 12 kPa', ...
          'DU145, 1 kPa', 'DU145, 3 kPa', 'DU145, 12 kPa',...
          'PC3, 1 kPa', 'PC3, 3 kPa', 'PC3, 12 kPa'}, 'widths', 0.3);

set(gca, 'FontSize', 15);
    
delete(findobj(gca,'type','text'))  % wipe out the boxplot labels that are problematic
set(gca,'xtick',1:2:7,'xticklabel',{ '1 kPa',  '3 kPa' '12 kPa'})  % and write the labels wanted
% xlim([0.5,4.5])
SE_scatter_ind1RV1 = 0.4*ones((length(time_aveSE_1RV1)),1);
f1 = scatter(SE_scatter_ind1RV1, time_aveSE_1RV1, 'ro');
f2 = scatter(0.4, mean(time_aveSE_1RV1), 'ko', 'LineWidth', 1);

SE_scatter_ind1LNCaP = 0.8*ones((length(time_aveSE_1LNCaP)),1);
f1 = scatter(SE_scatter_ind1LNCaP, time_aveSE_1LNCaP, 'mo');
f2 = scatter(0.8, mean(time_aveSE_1LNCaP), 'ko', 'LineWidth', 1);


SE_scatter_ind1DU145 = 1.2*ones((length(time_aveSE_1DU145)),1);
f1 = scatter(SE_scatter_ind1DU145, time_aveSE_1DU145, 'go');
f2 = scatter(1.2, mean(time_aveSE_1DU145), 'ko', 'LineWidth', 1);

SE_scatter_ind1PC3 = 1.6*ones((length(time_aveSE_1PC3)),1);
f1 = scatter(SE_scatter_ind1PC3, time_aveSE_1PC3, 'bo');
f2 = scatter(1.6, mean(time_aveSE_1PC3), 'ko', 'LineWidth', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SE_scatter_ind3RV1 = 2.4*ones((length(time_aveSE_3RV1)),1);
f1 = scatter(SE_scatter_ind3RV1, time_aveSE_3RV1, 'rd');
f2 = scatter(2.4, mean(time_aveSE_3RV1), 'kd', 'LineWidth', 1);

SE_scatter_ind3LNCaP = 2.8*ones((length(time_aveSE_3LNCaP)),1);
f1 = scatter(SE_scatter_ind3LNCaP, time_aveSE_3LNCaP, 'md');
f2 = scatter(2.8, mean(time_aveSE_3LNCaP), 'kd', 'LineWidth', 1);

SE_scatter_ind3DU145 = 3.2*ones((length(time_aveSE_3DU145)),1);
f1 = scatter(SE_scatter_ind3DU145, time_aveSE_3DU145, 'gd');
f2 = scatter(3.2, mean(time_aveSE_3DU145), 'kd', 'LineWidth', 1);

SE_scatter_ind3PC3 = 3.6*ones((length(time_aveSE_3PC3)),1);
f1 = scatter(SE_scatter_ind3PC3, time_aveSE_3PC3, 'bd');
f2 = scatter(3.6, mean(time_aveSE_3PC3), 'kd', 'LineWidth', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SE_scatter_ind12RV1 = 4.4*ones((length(time_aveSE_12RV1)),1);
f1 = scatter(SE_scatter_ind12RV1, time_aveSE_12RV1, 'r*');
f2 = scatter(4.4, mean(time_aveSE_12RV1), 'k*', 'LineWidth', 1);

SE_scatter_ind12LNCaP = 4.8*ones((length(time_aveSE_12LNCaP)),1);
f1 = scatter(SE_scatter_ind12LNCaP, time_aveSE_12LNCaP, 'm*');
f2 = scatter(4.8, mean(time_aveSE_12LNCaP), 'k*', 'LineWidth', 1);

SE_scatter_ind12DU145 = 5.2*ones((length(time_aveSE_12DU145)),1);
f1 = scatter(SE_scatter_ind12DU145, time_aveSE_12DU145, 'g*');
f2 = scatter(5.2, mean(time_aveSE_12DU145), 'k*', 'LineWidth', 1);


SE_scatter_ind12PC3 = 5.6*ones((length(time_aveSE_12PC3)),1);
f1 = scatter(SE_scatter_ind12PC3, time_aveSE_12PC3, 'b*');
f2 = scatter(5.6, mean(time_aveSE_12PC3), 'kd', 'LineWidth', 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% SE_scatter_ind50RV1 =6.4*ones((length(time_aveSE_50RV1)),1);
% f1 = scatter(SE_scatter_ind50RV1, time_aveSE_50RV1, 'r^');
% 
% 
% SE_scatter_ind50LNCaP = 6.8*ones((length(time_aveSE_50LNCaP)),1);
% f1 = scatter(SE_scatter_ind50LNCaP, time_aveSE_50LNCaP, 'm^');
% 
% 
% 
% SE_scatter_ind50DU145 = 7.2*ones((length(time_aveSE_50DU145)),1);
% f1 = scatter(SE_scatter_ind50DU145, time_aveSE_50DU145, 'g^');
% 
% 
% SE_scatter_ind50PC3 = 7.6*ones((length(time_aveSE_50PC3)),1);
% f1 = scatter(SE_scatter_ind50PC3, time_aveSE_50PC3, 'b^');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlim([0, 6]);


mlabel1 = scatter(10,1, 'ko'); 
mlabel2 = scatter(10,1, 'kd'); 
mlabel3 = scatter(10,1, 'k*'); 
% mlabel4 = scatter(10,1, 'k^'); 


legend([mlabel1, mlabel2, mlabel3, mlabel4], '1 kPa','3 kPa', '12 kPa', '50 kPa', 'Location', 'northeast')
set(gca, 'FontSize', 15, 'LineWidth', 1);
% ylim([0,16])
ylim([0,10])
saveas(gcf, 'Strain Energy (pJ)', 'tiff')




%%




%% Strain Energy - redundant
% 
% f_box_rmst_softness_compare = figure('Name', 'FBS_PC3/22RV1 SE Stiffness Compare');
% hold on
% ylabel('Strain Energy (pJ)');
% 
% boxplot([time_aveSE_1RV1, time_aveSE_3RV1, time_aveSE_12RV1, time_aveSE_50RV1, ...
%          time_aveSE_1LNCaP, time_aveSE_3LNCaP, time_aveSE_12LNCaP, time_aveSE_50LNCaP,...
%          time_aveSE_1DU145, time_aveSE_3DU145, time_aveSE_12DU145, time_aveSE_50DU145,...
%          time_aveSE_1PC3, time_aveSE_3PC3, time_aveSE_12PC3, time_aveSE_50PC3],...
%          [0.4*ones((length(time_aveSE_1RV1)),1); 0.8*ones((length(time_aveSE_3RV1)),1);... 
%           1.2*ones((length(time_aveSE_12RV1)),1); 1.6*ones((length( time_aveSE_50RV1)),1); ...
%           2.4*ones((length(time_aveSE_1LNCaP)),1); 2.8*ones((length(time_aveSE_3LNCaP)),1); ...
%           3.2*ones((length(time_aveSE_12LNCaP)),1); 3.6*ones((length( time_aveSE_50LNCaP)),1); ...
%           4.4*ones((length(time_aveSE_1DU145)),1); 4.8*ones((length(time_aveSE_3DU145)),1);...
%           5.2*ones((length(time_aveSE_12DU145)),1); 5.6*ones((length(time_aveSE_50DU145)),1);...
%           6.4*ones((length( time_aveSE_1PC3)),1); 6.8*ones((length( time_aveSE_3PC3)),1);......
%           7.2*ones((length(time_aveSE_12PC3)),1); 7.6*ones((length(time_aveSE_50PC3)),1)],...
%           'Positions', ...
%           [0.4*ones((length(time_aveSE_1RV1)),1); 0.8*ones((length(time_aveSE_3RV1)),1);... 
%           1.2*ones((length(time_aveSE_12RV1)),1); 1.6*ones((length( time_aveSE_50RV1)),1); ...
%           2.4*ones((length(time_aveSE_1LNCaP)),1); 2.8*ones((length(time_aveSE_3LNCaP)),1); ...
%           3.2*ones((length(time_aveSE_12LNCaP)),1); 3.6*ones((length(time_aveSE_50LNCaP)),1); ...
%           4.4*ones((length(time_aveSE_1DU145)),1); 4.8*ones((length(time_aveSE_3DU145)),1);...
%           5.2*ones((length(time_aveSE_12DU145)),1); 5.6*ones((length(time_aveSE_50DU145)),1);...
%           6.4*ones((length( time_aveSE_1PC3)),1); 6.8*ones((length( time_aveSE_3PC3)),1);...
%           7.2*ones((length(time_aveSE_12PC3)),1); 7.6*ones((length(time_aveSE_50PC3)),1)], ...
%           'Labels', ...
%           {'22RV1, 1 kPa', '22RV1, 3 kPa', '22RV1, 12 kPa', '22RV1, 50 kPa',...
%           'LNCaP, 1 kPa', 'LNCaP, 3 kPa', 'LNCaP, 12 kPa', 'LNCaP, 50 kPa',...
%           'DU145, 1 kPa', 'DU145, 3 kPa', 'DU145, 12 kPa', 'DU145, 50 kPa',...
%           'PC3, 1 kPa', 'PC3, 3 kPa', 'PC3, 12 kPa', 'PC3, 50 kPa'}, 'widths', 0.3)
% set(gca, 'FontSize', 15);
%     
% delete(findobj(gca,'type','text'))  % wipe out the boxplot labels that are problematic
% set(gca,'xtick',1:2:7,'xticklabel',{ '22RV1' 'LNCaP' 'DU145' 'PC3'})  % and write the labels wanted
% % xlim([0.5,4.5])
% SE_scatter_ind1RV1 = 0.4*ones((length(time_aveSE_1RV1)),1);
% f1 = scatter(SE_scatter_ind1RV1, time_aveSE_1RV1, 'ro');
% 
% SE_scatter_ind3RV1 = 0.8*ones((length(time_aveSE_3RV1)),1);
% f1 = scatter(SE_scatter_ind3RV1, time_aveSE_3RV1, 'rd');
% 
% SE_scatter_ind12RV1 = 1.2*ones((length(time_aveSE_12RV1)),1);
% f1 = scatter(SE_scatter_ind12RV1, time_aveSE_12RV1, 'r*');
% 
% SE_scatter_ind50RV1 =1.6*ones((length(time_aveSE_50RV1)),1);
% f1 = scatter(SE_scatter_ind50RV1, time_aveSE_50RV1, 'r^');
% %%%%%%%%%%%%%%%%%%%
% 
% SE_scatter_ind1LNCaP = 2.4*ones((length(time_aveSE_1LNCaP)),1);
% f1 = scatter(SE_scatter_ind1LNCaP, time_aveSE_1LNCaP, 'mo');
% 
% SE_scatter_ind3LNCaP = 2.8*ones((length(time_aveSE_3LNCaP)),1);
% f1 = scatter(SE_scatter_ind3LNCaP, time_aveSE_3LNCaP, 'md');
% 
% SE_scatter_ind12LNCaP = 3.2*ones((length(time_aveSE_12LNCaP)),1);
% f1 = scatter(SE_scatter_ind12LNCaP, time_aveSE_12LNCaP, 'm*');
% 
% SE_scatter_ind50LNCaP = 3.6*ones((length(time_aveSE_50LNCaP)),1);
% f1 = scatter(SE_scatter_ind50LNCaP, time_aveSE_50LNCaP, 'm^');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% SE_scatter_ind1DU145 = 4.4*ones((length(time_aveSE_1DU145)),1);
% f1 = scatter(SE_scatter_ind1DU145, time_aveSE_1DU145, 'go');
% 
% SE_scatter_ind3DU145 = 4.8*ones((length(time_aveSE_3DU145)),1);
% f1 = scatter(SE_scatter_ind3DU145, time_aveSE_3DU145, 'gd');
% 
% SE_scatter_ind12DU145 = 5.2*ones((length(time_aveSE_12DU145)),1);
% f1 = scatter(SE_scatter_ind12DU145, time_aveSE_12DU145, 'g*');
% 
% SE_scatter_ind50DU145 = 5.6*ones((length(time_aveSE_50DU145)),1);
% f1 = scatter(SE_scatter_ind50DU145, time_aveSE_50DU145, 'g^');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% SE_scatter_ind1PC3 = 6.4*ones((length(time_aveSE_1PC3)),1);
% f1 = scatter(SE_scatter_ind1PC3, time_aveSE_1PC3, 'bo');
% 
% SE_scatter_ind3PC3 = 6.8*ones((length(time_aveSE_3PC3)),1);
% f1 = scatter(SE_scatter_ind3PC3, time_aveSE_3PC3, 'bd');
% 
% SE_scatter_ind12PC3 = 7.2*ones((length(time_aveSE_12PC3)),1);
% f1 = scatter(SE_scatter_ind12PC3, time_aveSE_12PC3, 'b*');
% 
% SE_scatter_ind50PC3 = 7.6*ones((length(time_aveSE_50PC3)),1);
% f1 = scatter(SE_scatter_ind50PC3, time_aveSE_50PC3, 'b^');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% xlim([0, 8.2]);
% 
% 
% mlabel1 = scatter(10,1, 'ko'); 
% mlabel2 = scatter(10,1, 'kd'); 
% mlabel3 = scatter(10,1, 'k*'); 
% mlabel4 = scatter(10,1, 'k^'); 
% 
% 
% legend([mlabel1, mlabel2, mlabel3, mlabel4], '1 kPa','3 kPa', '12 kPa', '50 kPa', 'Location', 'northwest')
% set(gca, 'FontSize', 15, 'LineWidth', 1);
% % ylim([0,16])
% ylim([0,10])
% 
% 
% saveas(gcf, 'SE_Summary', 'tiff')

%%
% 
% f_box_rmst_softness_compare = figure('Name', 'FBS_PC3/22RV1 Displacement Compare');
% hold on
% ylabel('Average Displacement (\mum)');
% 
% boxplot([time_aveAD_1RV1, time_aveAD_3RV1, time_aveAD_12RV1, time_aveAD_50RV1, ...
%          time_aveAD_1LNCaP, time_aveAD_3LNCaP, time_aveAD_12LNCaP, time_aveAD_50LNCaP,...
%          time_aveAD_1DU145, time_aveAD_3DU145, time_aveAD_12DU145, time_aveAD_50DU145,...
%          time_aveAD_1PC3, time_aveAD_3PC3, time_aveAD_12PC3, time_aveAD_50PC3],...
%          [0.4*ones((length(time_aveAD_1RV1)),1); 0.8*ones((length(time_aveAD_3RV1)),1);... 
%           1.2*ones((length(time_aveAD_12RV1)),1); 1.6*ones((length( time_aveAD_50RV1)),1); ...
%           2.4*ones((length(time_aveAD_1LNCaP)),1); 2.8*ones((length(time_aveAD_3LNCaP)),1); ...
%           3.2*ones((length(time_aveAD_12LNCaP)),1); 3.6*ones((length( time_aveAD_50LNCaP)),1); ...
%           4.4*ones((length(time_aveAD_1DU145)),1); 4.8*ones((length(time_aveAD_3DU145)),1);...
%           5.2*ones((length(time_aveAD_12DU145)),1); 5.6*ones((length(time_aveAD_50DU145)),1);...
%           6.4*ones((length( time_aveAD_1PC3)),1); 6.8*ones((length( time_aveAD_3PC3)),1);......
%           7.2*ones((length(time_aveAD_12PC3)),1); 7.6*ones((length(time_aveAD_50PC3)),1)],...
%           'Positions', ...
%           [0.4*ones((length(time_aveAD_1RV1)),1); 0.8*ones((length(time_aveAD_3RV1)),1);... 
%           1.2*ones((length(time_aveAD_12RV1)),1); 1.6*ones((length( time_aveAD_50RV1)),1); ...
%           2.4*ones((length(time_aveAD_1LNCaP)),1); 2.8*ones((length(time_aveAD_3LNCaP)),1); ...
%           3.2*ones((length(time_aveAD_12LNCaP)),1); 3.6*ones((length(time_aveAD_50LNCaP)),1); ...
%           4.4*ones((length(time_aveAD_1DU145)),1); 4.8*ones((length(time_aveAD_3DU145)),1);...
%           5.2*ones((length(time_aveAD_12DU145)),1); 5.6*ones((length(time_aveAD_50DU145)),1);...
%           6.4*ones((length( time_aveAD_1PC3)),1); 6.8*ones((length( time_aveAD_3PC3)),1);...
%           7.2*ones((length(time_aveAD_12PC3)),1); 7.6*ones((length(time_aveAD_50PC3)),1)], ...
%           'Labels', ...
%           {'22RV1, 1 kPa', '22RV1, 3 kPa', '22RV1, 12 kPa', '22RV1, 50 kPa',...
%           'LNCaP, 1 kPa', 'LNCaP, 3 kPa', 'LNCaP, 12 kPa', 'LNCaP, 50 kPa',...
%           'DU145, 1 kPa', 'DU145, 3 kPa', 'DU145, 12 kPa', 'DU145, 50 kPa',...
%           'PC3, 1 kPa', 'PC3, 3 kPa', 'PC3, 12 kPa', 'PC3, 50 kPa'}, 'widths', 0.3)
% 
% set(gca, 'FontSize', 15);
%     
% delete(findobj(gca,'type','text'))  % wipe out the boxplot labels that are problematic
% set(gca,'xtick',1:2:7,'xticklabel',{ '22RV1' 'LNCaP' 'DU145' 'PC3'})  % and write the labels wanted
% % xlim([0.5,4.5])
% AD_scatter_ind1RV1 = 0.4*ones((length(time_aveAD_1RV1)),1);
% f1 = scatter(AD_scatter_ind1RV1, time_aveAD_1RV1, 'ro');
% 
% AD_scatter_ind3RV1 = 0.8*ones((length(time_aveAD_3RV1)),1);
% f1 = scatter(AD_scatter_ind3RV1, time_aveAD_3RV1, 'rd');
% 
% AD_scatter_ind12RV1 = 1.2*ones((length(time_aveAD_12RV1)),1);
% f1 = scatter(AD_scatter_ind12RV1, time_aveAD_12RV1, 'r*');
% 
% AD_scatter_ind50RV1 =1.6*ones((length(time_aveAD_50RV1)),1);
% f1 = scatter(AD_scatter_ind50RV1, time_aveAD_50RV1, 'r^');
% %%%%%%%%%%%%%%%%%%%
% 
% AD_scatter_ind1LNCaP = 2.4*ones((length(time_aveAD_1LNCaP)),1);
% f1 = scatter(AD_scatter_ind1LNCaP, time_aveAD_1LNCaP, 'mo');
% 
% AD_scatter_ind3LNCaP = 2.8*ones((length(time_aveAD_3LNCaP)),1);
% f1 = scatter(AD_scatter_ind3LNCaP, time_aveAD_3LNCaP, 'md');
% 
% AD_scatter_ind12LNCaP = 3.2*ones((length(time_aveAD_12LNCaP)),1);
% f1 = scatter(AD_scatter_ind12LNCaP, time_aveAD_12LNCaP, 'm*');
% 
% AD_scatter_ind50LNCaP = 3.6*ones((length(time_aveAD_50LNCaP)),1);
% f1 = scatter(AD_scatter_ind50LNCaP, time_aveAD_50LNCaP, 'm^');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% AD_scatter_ind1DU145 = 4.4*ones((length(time_aveAD_1DU145)),1);
% f1 = scatter(AD_scatter_ind1DU145, time_aveAD_1DU145, 'go');
% 
% AD_scatter_ind3DU145 = 4.8*ones((length(time_aveAD_3DU145)),1);
% f1 = scatter(AD_scatter_ind3DU145, time_aveAD_3DU145, 'gd');
% 
% AD_scatter_ind12DU145 = 5.2*ones((length(time_aveAD_12DU145)),1);
% f1 = scatter(AD_scatter_ind12DU145, time_aveAD_12DU145, 'g*');
% 
% AD_scatter_ind50DU145 = 5.6*ones((length(time_aveAD_50DU145)),1);
% f1 = scatter(AD_scatter_ind50DU145, time_aveAD_50DU145, 'g^');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% AD_scatter_ind1PC3 = 6.4*ones((length(time_aveAD_1PC3)),1);
% f1 = scatter(AD_scatter_ind1PC3, time_aveAD_1PC3, 'bo');
% 
% AD_scatter_ind3PC3 = 6.8*ones((length(time_aveAD_3PC3)),1);
% f1 = scatter(AD_scatter_ind3PC3, time_aveAD_3PC3, 'bd');
% 
% AD_scatter_ind12PC3 = 7.2*ones((length(time_aveAD_12PC3)),1);
% f1 = scatter(AD_scatter_ind12PC3, time_aveAD_12PC3, 'b*');
% 
% AD_scatter_ind50PC3 = 7.6*ones((length(time_aveAD_50PC3)),1);
% f1 = scatter(AD_scatter_ind50PC3, time_aveAD_50PC3, 'b^');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% xlim([0, 8.2]);
% 
% 
% mlabel1 = scatter(10,1, 'ko'); 
% mlabel2 = scatter(10,1, 'kd'); 
% mlabel3 = scatter(10,1, 'k*'); 
% mlabel4 = scatter(10,1, 'k^'); 
% 
% 
% legend([mlabel1, mlabel2, mlabel3, mlabel4], '1 kPa','3 kPa', '12 kPa', '50 kPa', 'Location', 'northeast')
% set(gca, 'FontSize', 15, 'LineWidth', 1);
% % ylim([0,16])
% ylim([0,5])
% saveas(gcf, 'Average Strain_Energy_PerStiffness', 'tiff')

%% RMST


f_box_rmst_softness_compare = figure('Name', 'Clustered_stiffnesses_RMST');
hold on
ylabel('RMST (Pa)');

boxplot([time_aveRMST_1RV1,  time_aveRMST_1LNCaP, time_aveRMST_1DU145,time_aveRMST_1PC3, ...
         time_aveRMST_3RV1, time_aveRMST_3LNCaP, time_aveRMST_3DU145,  time_aveRMST_3PC3,...
         time_aveRMST_12RV1, time_aveRMST_12LNCaP, time_aveRMST_12DU145,  time_aveRMST_12PC3],... time_aveRMST_50RV1, time_aveRMST_50LNCaP, time_aveRMST_50DU145, time_aveRMST_50PC3
         [0.4*ones((length(time_aveRMST_1RV1)),1); 0.8*ones((length(time_aveRMST_1LNCaP)),1);... 
          1.2*ones((length(time_aveRMST_1DU145)),1); 1.6*ones((length( time_aveRMST_1PC3)),1); ...
          2.4*ones((length(time_aveRMST_3RV1)),1); 2.8*ones((length( time_aveRMST_3LNCaP)),1); ...
          3.2*ones((length(time_aveRMST_3DU145)),1); 3.6*ones((length( time_aveRMST_3PC3)),1); ...
          4.4*ones((length( time_aveRMST_12RV1)),1); 4.8*ones((length( time_aveRMST_12LNCaP)),1);...
          5.2*ones((length(time_aveRMST_12DU145)),1); 5.6*ones((length(time_aveRMST_12PC3)),1);...
%           6.4*ones((length( time_aveRMST_50RV1)),1); 6.8*ones((length(time_aveRMST_50LNCaP)),1);...
%           7.2*ones((length(time_aveRMST_50DU145)),1); 7.6*ones((length(time_aveRMST_50PC3)),1)
            ],...
          'Positions', ...
           [0.4*ones((length(time_aveRMST_1RV1)),1); 0.8*ones((length(time_aveRMST_1LNCaP)),1);... 
          1.2*ones((length(time_aveRMST_1DU145)),1); 1.6*ones((length(time_aveRMST_1PC3)),1); ...
          2.4*ones((length(time_aveRMST_3RV1)),1); 2.8*ones((length(time_aveRMST_3LNCaP)),1); ...
          3.2*ones((length(time_aveRMST_3DU145)),1); 3.6*ones((length(time_aveRMST_3PC3)),1); ...
          4.4*ones((length( time_aveRMST_12RV1)),1); 4.8*ones((length(time_aveRMST_12LNCaP)),1);...
          5.2*ones((length(time_aveRMST_12DU145)),1); 5.6*ones((length(time_aveRMST_12PC3)),1);...
%           6.4*ones((length(time_aveRMST_50RV1)),1); 6.8*ones((length(time_aveRMST_50LNCaP)),1);......
%           7.2*ones((length(time_aveRMST_50DU145)),1); 7.6*ones((length(time_aveRMST_50PC3)),1)
           ],...
          'Labels', ...
          {'22RV1, 1 kPa', '22RV1, 3 kPa', '22RV1, 12 kPa'... , '22RV1, 50 kPa',
          'LNCaP, 1 kPa', 'LNCaP, 3 kPa', 'LNCaP, 12 kPa', ...
          'DU145, 1 kPa', 'DU145, 3 kPa', 'DU145, 12 kPa',...
          'PC3, 1 kPa', 'PC3, 3 kPa', 'PC3, 12 kPa'}, 'widths', 0.3);

set(gca, 'FontSize', 15);
    
delete(findobj(gca,'type','text'))  % wipe out the boxplot labels that are problematic
set(gca,'xtick',1:2:7,'xticklabel',{ '1 kPa',  '3 kPa' '12 kPa'})  % and write the labels wanted
% xlim([0.5,4.5])
RMST_scatter_ind1RV1 = 0.4*ones((length(time_aveRMST_1RV1)),1);
f1 = scatter(RMST_scatter_ind1RV1, time_aveRMST_1RV1, 'ro');
f2 = scatter(0.4, mean(time_aveRMST_1RV1), 'ko', 'LineWidth', 1);

RMST_scatter_ind1LNCaP = 0.8*ones((length(time_aveRMST_1LNCaP)),1);
f1 = scatter(RMST_scatter_ind1LNCaP, time_aveRMST_1LNCaP, 'mo');
f2 = scatter(0.8, mean(time_aveRMST_1LNCaP), 'ko', 'LineWidth', 1);


RMST_scatter_ind1DU145 = 1.2*ones((length(time_aveRMST_1DU145)),1);
f1 = scatter(RMST_scatter_ind1DU145, time_aveRMST_1DU145, 'go');
f2 = scatter(1.2, mean(time_aveRMST_1DU145), 'ko', 'LineWidth', 1);

RMST_scatter_ind1PC3 = 1.6*ones((length(time_aveRMST_1PC3)),1);
f1 = scatter(RMST_scatter_ind1PC3, time_aveRMST_1PC3, 'bo');
f2 = scatter(1.6, mean(time_aveRMST_1PC3), 'ko', 'LineWidth', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RMST_scatter_ind3RV1 = 2.4*ones((length(time_aveRMST_3RV1)),1);
f1 = scatter(RMST_scatter_ind3RV1, time_aveRMST_3RV1, 'rd');
f2 = scatter(2.4, mean(time_aveRMST_3RV1), 'kd', 'LineWidth', 1);

RMST_scatter_ind3LNCaP = 2.8*ones((length(time_aveRMST_3LNCaP)),1);
f1 = scatter(RMST_scatter_ind3LNCaP, time_aveRMST_3LNCaP, 'md');
f2 = scatter(2.8, mean(time_aveRMST_3LNCaP), 'kd', 'LineWidth', 1);

RMST_scatter_ind3DU145 = 3.2*ones((length(time_aveRMST_3DU145)),1);
f1 = scatter(RMST_scatter_ind3DU145, time_aveRMST_3DU145, 'gd');
f2 = scatter(3.2, mean(time_aveRMST_3DU145), 'kd', 'LineWidth', 1);

RMST_scatter_ind3PC3 = 3.6*ones((length(time_aveRMST_3PC3)),1);
f1 = scatter(RMST_scatter_ind3PC3, time_aveRMST_3PC3, 'bd');
f2 = scatter(3.6, mean(time_aveRMST_3PC3), 'kd', 'LineWidth', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RMST_scatter_ind12RV1 = 4.4*ones((length(time_aveRMST_12RV1)),1);
f1 = scatter(RMST_scatter_ind12RV1, time_aveRMST_12RV1, 'r*');
f2 = scatter(4.4, mean(time_aveRMST_12RV1), 'k*', 'LineWidth', 1);

RMST_scatter_ind12LNCaP = 4.8*ones((length(time_aveRMST_12LNCaP)),1);
f1 = scatter(RMST_scatter_ind12LNCaP, time_aveRMST_12LNCaP, 'm*');
f2 = scatter(4.8, mean(time_aveRMST_12LNCaP), 'k*', 'LineWidth', 1);

RMST_scatter_ind12DU145 = 5.2*ones((length(time_aveRMST_12DU145)),1);
f1 = scatter(RMST_scatter_ind12DU145, time_aveRMST_12DU145, 'g*');
f2 = scatter(5.2, mean(time_aveRMST_12DU145), 'k*', 'LineWidth', 1);


RMST_scatter_ind12PC3 = 5.6*ones((length(time_aveRMST_12PC3)),1);
f1 = scatter(RMST_scatter_ind12PC3, time_aveRMST_12PC3, 'b*');
f2 = scatter(5.6, mean(time_aveRMST_12PC3), 'kd', 'LineWidth', 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% RMST_scatter_ind50RV1 =6.4*ones((length(time_aveRMST_50RV1)),1);
% f1 = scatter(RMST_scatter_ind50RV1, time_aveRMST_50RV1, 'r^');
% 
% 
% RMST_scatter_ind50LNCaP = 6.8*ones((length(time_aveRMST_50LNCaP)),1);
% f1 = scatter(RMST_scatter_ind50LNCaP, time_aveRMST_50LNCaP, 'm^');
% 
% 
% 
% RMST_scatter_ind50DU145 = 7.2*ones((length(time_aveRMST_50DU145)),1);
% f1 = scatter(RMST_scatter_ind50DU145, time_aveRMST_50DU145, 'g^');
% 
% 
% RMST_scatter_ind50PC3 = 7.6*ones((length(time_aveRMST_50PC3)),1);
% f1 = scatter(RMST_scatter_ind50PC3, time_aveRMST_50PC3, 'b^');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlim([0, 6]);


mlabel1 = scatter(10,1, 'ko'); 
mlabel2 = scatter(10,1, 'kd'); 
mlabel3 = scatter(10,1, 'k*'); 
% mlabel4 = scatter(10,1, 'k^'); 


legend([mlabel1, mlabel2, mlabel3, mlabel4], '1 kPa','3 kPa', '12 kPa', '50 kPa', 'Location', 'northwest')
set(gca, 'FontSize', 15, 'LineWidth', 1);

ylim([0,400])
saveas(gcf, 'Clustered_Stiffness_RMST', 'tiff')


%%
% Stats_array_RMST = [time_aveRMST_1RV1, time_aveRMST_3RV1, time_aveRMST_12RV1, ... time_aveRMST_50PC3, ...
%          time_aveRMST_1LNCaP, time_aveRMST_3LNCaP, time_aveRMST_12LNCaP, ... time_aveRMST_50PC3, ...
%          time_aveRMST_1DU145, time_aveRMST_3DU145, time_aveRMST_12DU145, ... time_aveRMST_50PC3, ...
%          time_aveRMST_1PC3, time_aveRMST_3PC3, time_aveRMST_12PC3]; % time_aveRMST_50PC3, ...
% 
% groups = {[repmat({'RV1_1'}, 1, length(time_aveRMST_1RV1)), repmat({'RV1_3'}, 1, length(time_aveRMST_3RV1)), ...
%     repmat({'RV1_12'}, 1, length(time_aveRMST_12RV1)), ... repmat({'RV1_50'}, 1, length(time_aveRMST_50RV1)), ... 
%     repmat({'LNCaP_1'}, 1, length(time_aveRMST_1LNCaP)), repmat({'LNCaP_3'}, 1, length(time_aveRMST_3LNCaP)), ...
%     repmat({'LNCaP_12'}, 1, length(time_aveRMST_12LNCaP)), ... repmat({'LNCaP_50'}, 1, length(time_aveRMST_50LNCaP))]...
%     repmat({'DU145_1'}, 1, length(time_aveRMST_1DU145)), repmat({'DU145_3'}, 1, length(time_aveRMST_3DU145)), ...
%     repmat({'DU145_12'}, 1, length(time_aveRMST_12DU145)), ... repmat({'DU145_50'}, 1, length(time_aveRMST_50DU145)), ... 
%     repmat({'PC3_1'}, 1, length(time_aveRMST_1PC3)), repmat({'PC3_3'}, 1, length(time_aveRMST_3PC3)), ...
%     repmat({'PC3_12'}, 1, length(time_aveRMST_12PC3)), ... repmat({'PC3_50'}, 1, length(time_aveRMST_50PC3)), ... 
%     ]};
% 
%  [P, tbl, stats] = kruskalwallis(Stats_array_RMST, groups{1});
%  [c1, ~,~,gnames] = multcompare(stats, "CType","bonferroni"); 
% 
%  table1 = array2table(c1, "VariableNames", ... 
%      ["Group A", "Group B", "Lower Limit", "A-B", "Upper Limit", "P-value"]);
%  table1.("Group A") = gnames(table1.("Group A"));
%  table1.("Group B") = gnames(table1.("Group B"));
% disp(table1)
% 
% 
% Stats_array_SE = [time_aveSE_1RV1, time_aveSE_3RV1, time_aveSE_12RV1, ... time_aveSE_50PC3, ...
%          time_aveSE_1LNCaP, time_aveSE_3LNCaP, time_aveSE_12LNCaP, ... time_aveSE_50PC3, ...
%          time_aveSE_1DU145, time_aveSE_3DU145, time_aveSE_12DU145, ... time_aveSE_50PC3, ...
%          time_aveSE_1PC3, time_aveSE_3PC3, time_aveSE_12PC3]; % time_aveSE_50PC3, ...
% 
% groups = {[repmat({'RV1_1'}, 1, length(time_aveSE_1RV1)), repmat({'RV1_3'}, 1, length(time_aveSE_3RV1)), ...
%     repmat({'RV1_12'}, 1, length(time_aveSE_12RV1)), ... repmat({'RV1_50'}, 1, length(time_aveSE_50RV1)), ... 
%     repmat({'LNCaP_1'}, 1, length(time_aveSE_1LNCaP)), repmat({'LNCaP_3'}, 1, length(time_aveSE_3LNCaP)), ...
%     repmat({'LNCaP_12'}, 1, length(time_aveSE_12LNCaP)), ... repmat({'LNCaP_50'}, 1, length(time_aveSE_50LNCaP))]...
%     repmat({'DU145_1'}, 1, length(time_aveSE_1DU145)), repmat({'DU145_3'}, 1, length(time_aveSE_3DU145)), ...
%     repmat({'DU145_12'}, 1, length(time_aveSE_12DU145)), ... repmat({'DU145_50'}, 1, length(time_aveSE_50DU145)), ... 
%     repmat({'PC3_1'}, 1, length(time_aveSE_1PC3)), repmat({'PC3_3'}, 1, length(time_aveSE_3PC3)), ...
%     repmat({'PC3_12'}, 1, length(time_aveSE_12PC3)), ... repmat({'PC3_50'}, 1, length(time_aveSE_50PC3)), ... 
%     ]};
% 
%  [P, tbl, stats] = kruskalwallis(Stats_array_SE, groups{1});
%  [c1, ~,~,gnames] = multcompare(stats, "CType","dunn-sidak"); 
% 
%  table1 = array2table(c1, "VariableNames", ... 
%      ["Group A", "Group B", "Lower Limit", "A-B", "Upper Limit", "P-value"]);
%  table1.("Group A") = gnames(table1.("Group A"));
%  table1.("Group B") = gnames(table1.("Group B"));
% disp(table1)
% Stats_array_SE = [time_aveRMST_1RV1,  time_aveRMST_1DU145, time_aveRMST_1DU145,time_aveRMST_1PC3, ...[
%          time_aveRMST_3RV1, time_aveRMST_3LNCaP, time_aveRMST_3DU145,  time_aveRMST_3PC3,...
%          time_aveRMST_12RV1, time_aveRMST_12LNCaP, time_aveRMST_12DU145,  time_aveRMST_12PC3);

%% PAIRWISE MANN-WHITNEY 
% SE 
Cell_array_SE = {time_aveSE_1RV1, time_aveSE_3RV1, time_aveSE_12RV1, ... time_aveSE_50PC3, ...
         time_aveSE_1LNCaP, time_aveSE_3LNCaP, time_aveSE_12LNCaP, ... time_aveSE_50PC3, ...
         time_aveSE_1DU145, time_aveSE_3DU145, time_aveSE_12DU145, ... time_aveSE_50PC3, ...
         time_aveSE_1PC3, time_aveSE_3PC3, time_aveSE_12PC3};

p_mat_SE = zeros(length(Cell_array_SE));
for i = 1:length(Cell_array_SE)
    for j = 1:length(Cell_array_SE) 
    p_mat_SE(i,j) = ranksum(Cell_array_SE{i}, Cell_array_SE{j});
    end 
end 

% RMST
Cell_array_RMST = {time_aveRMST_1RV1, time_aveRMST_3RV1, time_aveRMST_12RV1, ... time_aveRMST_50PC3, ...
         time_aveRMST_1LNCaP, time_aveRMST_3LNCaP, time_aveRMST_12LNCaP, ... time_aveRMST_50PC3, ...
         time_aveRMST_1DU145, time_aveRMST_3DU145, time_aveRMST_12DU145, ... time_aveRMST_50PC3, ...
         time_aveRMST_1PC3, time_aveRMST_3PC3, time_aveRMST_12PC3};

p_mat_RMST = zeros(length(Cell_array_RMST));
for i = 1:length(Cell_array_RMST)
    for j = 1:length(Cell_array_RMST) 
    p_mat_RMST(i,j) = ranksum(Cell_array_RMST{i}, Cell_array_RMST{j});
    end 
end 

%%% Average Displacement 

Cell_array_AD = {time_aveAD_1RV1, time_aveAD_3RV1, time_aveAD_12RV1, ... time_aveAD_50PC3, ...
         time_aveAD_1LNCaP, time_aveAD_3LNCaP, time_aveAD_12LNCaP, ... time_aveAD_50PC3, ...
         time_aveAD_1DU145, time_aveAD_3DU145, time_aveAD_12DU145, ... time_aveAD_50PC3, ...
         time_aveAD_1PC3, time_aveAD_3PC3, time_aveAD_12PC3};

p_mat_AD = zeros(length(Cell_array_AD));
for i = 1:length(Cell_array_AD)
    for j = 1:length(Cell_array_AD) 
    p_mat_AD(i,j) = ranksum(Cell_array_AD{i}, Cell_array_AD{j});
    end 
end 



% f_box_SE_softness_compare = figure('Name', 'Clustered_stiffnesses SE');
% hold on
% ylabel('Strain Energy (pJ)');
% 
% boxplot([time_aveSE_1RV1,  time_aveSE_1PC3, time_aveSE_1DU145,time_aveSE_1PC3, ...
%          time_aveSE_3RV1, time_aveSE_3LNCaP, time_aveSE_3DU145,  time_aveSE_3PC3,...
%          time_aveSE_12RV1, time_aveSE_12LNCaP, time_aveSE_12DU145,  time_aveSE_12PC3,...
%          time_aveSE_50RV1, time_aveSE_50LNCaP, time_aveSE_50DU145, time_aveSE_50PC3],...
%          [0.4*ones((length(time_aveSE_1RV1)),1); 0.8*ones((length(time_aveSE_1LNCaP)),1);... 
%           1.2*ones((length(time_aveSE_1DU145)),1); 1.6*ones((length( time_aveSE_1PC3)),1); ...
%           2.4*ones((length(time_aveSE_3RV1)),1); 2.8*ones((length( time_aveSE_3LNCaP)),1); ...
%           3.2*ones((length(time_aveSE_3DU145)),1); 3.6*ones((length( time_aveSE_3PC3)),1); ...
%           4.4*ones((length( time_aveSE_12RV1)),1); 4.8*ones((length( time_aveSE_12LNCaP)),1);...
%           5.2*ones((length(time_aveSE_12DU145)),1); 5.6*ones((length(time_aveSE_12PC3)),1);...
%           6.4*ones((length( time_aveSE_50RV1)),1); 6.8*ones((length(time_aveSE_50LNCaP)),1);...
%           7.2*ones((length(time_aveSE_50DU145)),1); 7.6*ones((length(time_aveSE_50PC3)),1)],...
%           'Positions', ...
%            [0.4*ones((length(time_aveSE_1RV1)),1); 0.8*ones((length(time_aveSE_1LNCaP)),1);... 
%           1.2*ones((length(time_aveSE_1DU145)),1); 1.6*ones((length(time_aveSE_1PC3)),1); ...
%           2.4*ones((length(time_aveSE_3RV1)),1); 2.8*ones((length(time_aveSE_3LNCaP)),1); ...
%           3.2*ones((length(time_aveSE_3DU145)),1); 3.6*ones((length(time_aveSE_3PC3)),1); ...
%           4.4*ones((length( time_aveSE_12RV1)),1); 4.8*ones((length(time_aveSE_12LNCaP)),1);...
%           5.2*ones((length(time_aveSE_12DU145)),1); 5.6*ones((length(time_aveSE_12PC3)),1);...
%           6.4*ones((length(time_aveSE_50RV1)),1); 6.8*ones((length(time_aveSE_50LNCaP)),1);......
%           7.2*ones((length(time_aveSE_50DU145)),1); 7.6*ones((length(time_aveSE_50PC3)),1)],...
%           'Labels', ...
%           {'22RV1, 1 kPa', '22RV1, 3 kPa', '22RV1, 12 kPa', '22RV1, 50 kPa',...
%           'LNCaP, 1 kPa', 'LNCaP, 3 kPa', 'LNCaP, 12 kPa', 'LNCaP, 50 kPa',...
%           'DU145, 1 kPa', 'DU145, 3 kPa', 'DU145, 12 kPa', 'DU145, 50 kPa',...
%           'PC3, 1 kPa', 'PC3, 3 kPa', 'PC3, 12 kPa', 'PC3, 50 kPa'}, 'widths', 0.3);
% 
% set(gca, 'FontSize', 15);
%     
% delete(findobj(gca,'type','text'))  % wipe out the boxplot labels that are problematic
% set(gca,'xtick',1:2:7,'xticklabel',{ '1 kPa' '3 kPa' ' 12 kPa' '50 kPa'})  % and write the labels wanted
% % xlim([0.5,4.5])
% SE_scatter_ind1RV1 = 0.4*ones((length(time_aveSE_1RV1)),1);
% f1 = scatter(SE_scatter_ind1RV1, time_aveSE_1RV1, 'ro');
% 
% SE_scatter_ind1LNCaP = 0.8*ones((length(time_aveSE_1LNCaP)),1);
% f1 = scatter(SE_scatter_ind1LNCaP, time_aveSE_1LNCaP, 'mo');
% 
% SE_scatter_ind1DU145 = 1.2*ones((length(time_aveSE_1DU145)),1);
% f1 = scatter(SE_scatter_ind1DU145, time_aveSE_1DU145, 'go');
% 
% 
% SE_scatter_ind1PC3 = 1.6*ones((length(time_aveSE_1PC3)),1);
% f1 = scatter(SE_scatter_ind1PC3, time_aveSE_1PC3, 'bo');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% SE_scatter_ind3RV1 = 2.4*ones((length(time_aveSE_3RV1)),1);
% f1 = scatter(SE_scatter_ind3RV1, time_aveSE_3RV1, 'rd');
% 
% 
% SE_scatter_ind3LNCaP = 2.8*ones((length(time_aveSE_3LNCaP)),1);
% f1 = scatter(SE_scatter_ind3LNCaP, time_aveSE_3LNCaP, 'md');
% 
% SE_scatter_ind3DU145 = 3.2*ones((length(time_aveSE_3DU145)),1);
% f1 = scatter(SE_scatter_ind3DU145, time_aveSE_3DU145, 'gd');
% 
% 
% SE_scatter_ind3PC3 = 3.6*ones((length(time_aveSE_3PC3)),1);
% f1 = scatter(SE_scatter_ind3PC3, time_aveSE_3PC3, 'bd');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% SE_scatter_ind12RV1 = 4.4*ones((length(time_aveSE_12RV1)),1);
% f1 = scatter(SE_scatter_ind12RV1, time_aveSE_12RV1, 'r*');
% 
% 
% SE_scatter_ind12LNCaP = 4.8*ones((length(time_aveSE_12LNCaP)),1);
% f1 = scatter(SE_scatter_ind12LNCaP, time_aveSE_12LNCaP, 'm*');
% 
% SE_scatter_ind12DU145 = 5.2*ones((length(time_aveSE_12DU145)),1);
% f1 = scatter(SE_scatter_ind12DU145, time_aveSE_12DU145, 'g*');
% 
% 
% SE_scatter_ind12PC3 = 5.6*ones((length(time_aveSE_12PC3)),1);
% f1 = scatter(SE_scatter_ind12PC3, time_aveSE_12PC3, 'b*');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% SE_scatter_ind50RV1 =6.4*ones((length(time_aveSE_50RV1)),1);
% f1 = scatter(SE_scatter_ind50RV1, time_aveSE_50RV1, 'r^');
% 
% 
% SE_scatter_ind50LNCaP = 6.8*ones((length(time_aveSE_50LNCaP)),1);
% f1 = scatter(SE_scatter_ind50LNCaP, time_aveSE_50LNCaP, 'm^');
% 
% 
% 
% SE_scatter_ind50DU145 = 7.2*ones((length(time_aveSE_50DU145)),1);
% f1 = scatter(SE_scatter_ind50DU145, time_aveSE_50DU145, 'g^');
% 
% 
% SE_scatter_ind50PC3 = 7.6*ones((length(time_aveSE_50PC3)),1);
% f1 = scatter(SE_scatter_ind50PC3, time_aveSE_50PC3, 'b^');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% xlim([0, 6.2]);
% 
% 
% mlabel1 = plot(10,1, 'r-'); 
% mlabel2 = plot(10,1, 'm-'); 
% mlabel3 = plot(10,1, 'g-'); 
% mlabel4 = plot(10,1, 'b-'); 
% 
% 
% legend([mlabel1, mlabel2, mlabel3, mlabel4], '1 kPa','3 kPa', '12 kPa', '50 kPa', 'Location', 'northeast', 'Linewidth', 1)
% set(gca, 'FontSize', 15, 'LineWidth', 1);
% ylim([0,7])
% saveas(gcf, 'SE_Stiff_Clustered', 'tiff')
% 

% RMST 
Cell_array_SE = {time_aveSE_1RV1, time_aveSE_3RV1, time_aveSE_12RV1, ... time_aveSE_50PC3, ...
         time_aveSE_1LNCaP, time_aveSE_3LNCaP, time_aveSE_12LNCaP, ... time_aveSE_50PC3, ...
         time_aveSE_1DU145, time_aveSE_3DU145, time_aveSE_12DU145, ... time_aveSE_50PC3, ...
         time_aveSE_1PC3, time_aveSE_3PC3, time_aveSE_12PC3};

p_mat_SE = zeros(length(Cell_array_SE));
for i = 1:length(Cell_array_SE)
    for j = 1:length(Cell_array_SE) 
    p_mat_SE(i,j) = ranksum(Cell_array_SE{i}, Cell_array_SE{j});
    end 
end 


Cell_array_SE = {time_aveSE_1RV1, time_aveSE_3RV1, time_aveSE_12RV1, ... time_aveSE_50PC3, ...
         time_aveSE_1LNCaP, time_aveSE_3LNCaP, time_aveSE_12LNCaP, ... time_aveSE_50PC3, ...
         time_aveSE_1DU145, time_aveSE_3DU145, time_aveSE_12DU145, ... time_aveSE_50PC3, ...
         time_aveSE_1PC3, time_aveSE_3PC3, time_aveSE_12PC3};

p_mat_SE = zeros(length(Cell_array_SE));
for i = 1:length(Cell_array_SE)
    for j = 1:length(Cell_array_SE) 
    p_mat_SE(i,j) = ranksum(Cell_array_SE{i}, Cell_array_SE{j});
    end 
end






%%
%  
% f_box_SE_stiffness_compare = figure('Name','FBS_PC3/22RV1 SE Stiffness Compare');
% hold on
% ylabel('Strain Energy (pJ)');
% boxplot([time_aveSE_3RV1, time_aveSE_12RV1,time_aveSE_3LNCaP, time_aveSE_12LNCaP,   time_aveSE_3DU145, time_aveSE_12DU145 time_aveSE_3PC3, time_aveSE_12PC3], ...
%         [0.8*ones((length(time_aveSE_3RV1)),1);... 
%           1.2*ones((length(time_aveSE_12RV1)),1); 1.8*ones((length(time_aveSE_3LNCaP)),1);...
%           2.2*ones((length(time_aveSE_12LNCaP)),1);2.8*ones((length(time_aveSE_3DU145)),1);...
%           3.2*ones((length(time_aveSE_12DU145)),1); 3.8*ones((length(time_aveSE_3PC3)),1);...
%           4.2*ones((length(time_aveSE_12PC3)),1)],'Positions', [0.8*ones((length(time_aveSE_3RV1)),1);... 
%           1.2*ones((length(time_aveSE_12RV1)),1); 1.8*ones((length(time_aveSE_3LNCaP)),1);...
%           2.2*ones((length(time_aveSE_12LNCaP)),1);2.8*ones((length(time_aveSE_3DU145)),1);...
%           3.2*ones((length(time_aveSE_12DU145)),1); 3.8*ones((length(time_aveSE_3PC3)),1);...
%           4.2*ones((length(time_aveSE_12PC3)),1)], 'Labels', {'22RV1, 3 kPa', '22RV1, 12 kPa',...
%            'LNCaP, 3 kPa', 'LNCaP, 12 kPa', 'DU145, 3 kPa', 'DU145, 12 kPa', 'PC3, 3 kPa', '22RV1, 12 kPa'}, 'widths', 0.3);
% set(gca, 'FontSize', 15);
%       
% delete(findobj(gca,'type','text'))  % wipe out the boxplot labels that are problematic
% set(gca,'xtick',1:4,'xticklabel',{ '22RV1' 'LNCaP' 'DU145' 'PC3'})  % and write the labels wanted
% xlim([0.5,4.5])
% 
% SE_scatter_ind3RV1 = 0.8*ones((length(time_aveSE_3RV1)),1);
% f1 = scatter(SE_scatter_ind3RV1, time_aveSE_3RV1, 'ro');
% 
% SE_scatter_ind12RV1 = 1.2*ones((length(time_aveSE_12RV1)),1);
% f1 = scatter(SE_scatter_ind12RV1, time_aveSE_12RV1, 'r');
% 
% SE_scatter_ind3LNCaP = 1.8*ones((length(time_aveSE_3LNCaP)),1);
% f1 = scatter(SE_scatter_ind3LNCaP, time_aveSE_3LNCaP, 'mo');
% 
% SE_scatter_ind12LNCaP = 2.2*ones((length(time_aveSE_12LNCaP)),1);
% f1 = scatter(SE_scatter_ind12LNCaP, time_aveSE_12LNCaP, 'md');
% 
% SE_scatter_ind3DU145 = 2.8*ones((length(time_aveSE_3DU145)),1);
% f1 = scatter(SE_scatter_ind3DU145, time_aveSE_3DU145, 'go');
% 
% SE_scatter_ind12DU145 = 3.2*ones((length(time_aveSE_12DU145)),1);
% f1 = scatter(SE_scatter_ind12DU145, time_aveSE_12DU145, 'gd');
% 
% 
% SE_scatter_ind3PC3 = 3.8*ones((length(time_aveSE_3PC3)),1);
% f1 = scatter(SE_scatter_ind3PC3, time_aveSE_3PC3, 'bo');
% 
% SE_scatter_ind12PC3 = 4.2*ones((length(time_aveSE_12PC3)),1);
% f1 = scatter(SE_scatter_ind12PC3, time_aveSE_12PC3, 'bd');
% 
% xlim([0.25,4.75]);
% mlabel1 = scatter(10,1, 'ko'); 
% mlabel2 = scatter(10,1, 'kd'); 
% 
% legend([mlabel1, mlabel2], '3 kPa', '12 kPa', 'Location', 'northwest')
% set(gca, 'FontSize', 15, 'LineWidth', 1);
% ylim([0,16])
% % ylim([0,500])
% 
% saveas(gcf, 'FBS_f_box_SE_stiffness_compare', 'tiff')




%% CSS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% SECTION 3: 22RV1 vs. PC3 vs. DU145 vs. LNCaP as a Function of Stiffness in CSS 
% f_stiff_rmst = figure('Name', 'RMST vs. Time, 22RV1 & PC3 12kPa');
% ylabel('RMST (Pa)');
% xlabel('Time (Minutes)'); 
% 
% f_stiff_SE = figure('Name', 'SE vs. Time 22RV1 vs. PC3 12kPa');
% ylabel('Strain  Energy (pJ)');
% xlabel('Time (Minutes)'); 
%   
% % Plot RMST and SE of 12 kPa 22RV1
% 
% % Identify Relevant Exps
% to_open = text_cats(RV1_12kPa_72_CSS_ind, 9);
% 
% % Identify tie interval
% time_interval = cell2mat(text_cats(RV1_12kPa_72_CSS_ind,10))/60;
% frame_add = cell2mat(text_cats(RV1_12kPa_72_CSS_ind,6));
% 
% % Identify bad positions
% series_omit_array = text_cats(RV1_12kPa_72_CSS_ind,7);
% 
% % Run through each position and plot 
% q=1;
% 
% for i = 1:length(to_open) 
%     % Open file and extract values rxc = pos x time
%   data_temp = open(['Mats_So_Far/' to_open{i}]);
%   rmst = data_temp.rmst';
%   strain_energy = data_temp.strain_energy'*10^12; % Adjust to pJ
%  average_displacement = data_temp.average_displacement';
%   
%   % Identify positions to omit and do so
%   if series_omit_array{i} ~=0 
%     try 
%     omit_indices = str2double(strsplit(series_omit_array{i}, ','));
%     catch
%     omit_indices = series_omit_array{i}; 
%     end 
%     rmst(omit_indices,:) = []; 
%     strain_energy(omit_indices, :) = [];
%     average_displacement(omit_indices,:) = [];
%   end     
% 
% 
%   % Identify the frame at which the treatment was applied 
%   [pos_num, time_num] = size(rmst);
%   if frame_add(i) ~= 0
%       time_num = frame_add(i) - 1;
%   end 
%   t = ([1:time_num]-1)*time_interval(i); 
% 
%   for j = 1:pos_num  
%   figure(f_stiff_rmst)  
%   hold on
%   p1 = plot(t, rmst(j,1:length(t)), 'r'); 
%   time_aveRMST_12RV1(q) = mean(rmst(j, 1:length(t)));
% % text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
%   figure(f_stiff_SE) 
%   hold on
%   q1 =  plot(t, strain_energy(j,1:length(t)), 'r'); 
%   time_aveSE_12RV1(q) = mean(strain_energy(j, 1:length(t)));
%   q = q+1;
%   end 
% end 
% 
% % Plot RMST and SE of 12 kPa PC3
% % Identify Relevant Exps
% to_open = text_cats(PC3_12kPa_72_CSS_ind , 9);
% 
% % Identify time interval
% time_interval = cell2mat(text_cats(PC3_12kPa_72_CSS_ind ,10))/60;
% 
% % Identify when the frame for adding 
% frame_add = cell2mat(text_cats(PC3_12kPa_72_CSS_ind ,6));
% 
% % Identify bad positions
% series_omit_array = text_cats(PC3_12kPa_72_CSS_ind ,7);
% 
% % Run through each position and plot 
% q=1;
% for i = 1:length(to_open) 
%     % Open file and extract values rxc = pos x time
%   data_temp = open(['Mats_So_Far/' to_open{i}]);
%   rmst = data_temp.rmst';
%   strain_energy = data_temp.strain_energy'*10^12;
%    average_displacement = data_temp.average_displacement';
%   
%   % Identify positions to omit and do so
%   if series_omit_array{i} ~=0 
%     try 
%     omit_indices = str2double(strsplit(series_omit_array{i}, ','));
%     catch
%     omit_indices = series_omit_array{i}; 
%     end 
%     rmst(omit_indices,:) = []; 
%     strain_energy(omit_indices, :) = [];
%     average_displacement(omit_indices,:) = [];
%   end     
% 
% 
%   [pos_num, time_num] = size(rmst);
%   if frame_add(i) ~= 0
%       time_num = frame_add(i) - 1;
%   end 
%   t = ([1:time_num]-1)*time_interval(i); 
%   
%   for j = 1:pos_num  
%       figure(f_stiff_rmst)  
%       hold on
%       p2 = plot(t, rmst(j,1:length(t)), 'b');
%       time_aveRMST_12PC3(q) = mean(rmst(j, 1:length(t)));
% % text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
%       figure(f_stiff_SE)
%       hold on
%       q2 = plot(t, strain_energy(j,1:length(t)), 'b');
%       time_aveSE_12PC3(q) = mean(strain_energy(j, 1:length(t)));
%   q = q+1;
%   end 
% end 
%    
% % Plot RMST and SE of 12 kPa DU145 in CSS 
% % Identify Relevant Exps
% to_open = text_cats(DU145_12kPa_72_CSS_ind, 9);
% 
% % Identify time interval
% time_interval = cell2mat(text_cats(DU145_12kPa_72_CSS_ind ,10))/60;
% 
% % Identify when the frame for adding 
% frame_add = cell2mat(text_cats(DU145_12kPa_72_CSS_ind ,6));
% 
% % Identify bad positions
% series_omit_array = text_cats(DU145_12kPa_72_CSS_ind ,7);
% 
% % Run through each position and plot 
% q=1;
% for i = 1:length(to_open) 
%     % Open file and extract values rxc = pos x time
%   data_temp = open(['Mats_So_Far/' to_open{i}]);
%   rmst = data_temp.rmst';
%   strain_energy = data_temp.strain_energy'*10^12;
%    average_displacement = data_temp.average_displacement';
%   
%   % Identify positions to omit and do so
%   if series_omit_array{i} ~=0 
%     try 
%     omit_indices = str2double(strsplit(series_omit_array{i}, ','));
%     catch
%     omit_indices = series_omit_array{i}; 
%     end 
%     rmst(omit_indices,:) = []; 
%     strain_energy(omit_indices, :) = [];
%     average_displacement(omit_indices,:) = [];
%   end     
% 
% 
%   [pos_num, time_num] = size(rmst);
%   if frame_add(i) ~= 0
%       time_num = frame_add(i) - 1;
%   end 
%   t = ([1:time_num]-1)*time_interval(i); 
%   
%   for j = 1:pos_num  
%       figure(f_stiff_rmst)  
%       hold on
%       p3 = plot(t, rmst(j,1:length(t)), 'g');
%       time_aveRMST_12DU145(q) = mean(rmst(j, 1:length(t)));
% % text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
%       figure(f_stiff_SE)
%       hold on
%       q3 = plot(t, strain_energy(j,1:length(t)), 'g');
%       time_aveSE_12DU145(q) = mean(strain_energy(j, 1:length(t)));
%   q = q+1;
%   end 
%   
% end 
% 
%   
% % Plot RMST and SE of 12 kPa DU145 in CSS 
% % Identify Relevant Exps
% to_open = text_cats(LNCaP_12kPa_72_CSS_ind, 9);
% 
% % Identify time interval
% time_interval = cell2mat(text_cats(LNCaP_12kPa_72_CSS_ind ,10))/60;
% 
% % Identify when the frame for adding 
% frame_add = cell2mat(text_cats(LNCaP_12kPa_72_CSS_ind ,6));
% 
% % Identify bad positions
% series_omit_array = text_cats(LNCaP_12kPa_72_CSS_ind ,7);
% 
% % Run through each position and plot 
% q=1;
% for i = 1:length(to_open) 
%     % Open file and extract values rxc = pos x time
%   data_temp = open(['Mats_So_Far/' to_open{i}]);
%   rmst = data_temp.rmst';
%   strain_energy = data_temp.strain_energy'*10^12;
%    average_displacement = data_temp.average_displacement';
%   
%   % Identify positions to omit and do so
%   if series_omit_array{i} ~=0 
%     try 
%     omit_indices = str2double(strsplit(series_omit_array{i}, ','));
%     catch
%     omit_indices = series_omit_array{i}; 
%     end 
%     rmst(omit_indices,:) = []; 
%     strain_energy(omit_indices, :) = [];
%     average_displacement(omit_indices,:) = [];
%   end     
% 
% 
%   [pos_num, time_num] = size(rmst);
%   if frame_add(i) ~= 0
%       time_num = frame_add(i) - 1;
%   end 
%   t = ([1:time_num]-1)*time_interval(i); 
%   
%   for j = 1:pos_num  
%       figure(f_stiff_rmst)  
%       hold on
%       p4 = plot(t, rmst(j,1:length(t)), 'm');
%       time_aveRMST_12LNCaP(q) = mean(rmst(j, 1:length(t)));
% % text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
%       figure(f_stiff_SE)
%       hold on
%       q4 = plot(t, strain_energy(j,1:length(t)), 'm');
%       time_aveSE_12LNCaP(q) = mean(strain_energy(j, 1:length(t)));
%   q = q+1;
%   end 
%   
% end 
% 
% 
% hold off 
% figure(f_stiff_rmst)
% hold on
% legend([p1,p2], '22RV1, 12 kPa', 'PC3, 12 kPa')
% saveas(gcf, 'stiff_rmst.tiff');
% 
% figure(f_stiff_SE)
% hold on
% legend([q1,q2], '22RV1, 12 kPa', 'PC3, 12 kPa')
% saveas(gcf, 'stiff_SE.tiff');
% %% RMST and Strain Energy Boxplots STIFF
% f_box_stiff_rmst = figure('Name', 'RMST boxplot 22RV1 vs. PC3 12kPa');
% ylabel('RMST (Pa)');
% 
% f_box_stiff_SE = figure('Name', 'SE boxplot 22RV1 vs. PC3 12kPa');
% ylabel('Strain  Energy (pJ)');
% 
% figure(f_box_stiff_rmst)
% hold on 
% boxplot([time_aveRMST_12RV1,time_aveRMST_12LNCaP ,time_aveRMST_12DU145, time_aveRMST_12PC3], [1*ones((length(time_aveRMST_12RV1)),1); 2*ones((length(time_aveRMST_12LNCaP)),1); ...
%           3*ones((length(time_aveRMST_12DU145)),1); 4*ones((length(time_aveRMST_12PC3)),1)], 'Positions',[1*ones((length(time_aveRMST_12RV1)),1); 2*ones((length(time_aveRMST_12LNCaP)),1); ...
%           3*ones((length(time_aveRMST_12DU145)),1); 4*ones((length(time_aveRMST_12PC3)),1)],  'Labels', {'22RV1 12 kPa', 'LNCaP', 'DU145', 'PC3 12 kPa'});
%       
% RMST_scatter_ind12RV1 = 1*ones((length(time_aveRMST_12RV1)),1);
% f1 = scatter(RMST_scatter_ind12RV1, time_aveRMST_12RV1, 'r');
% 
% RMST_scatter_ind12DU145 = 2*ones((length(time_aveRMST_12LNCaP)),1);
% f1 = scatter(RMST_scatter_ind12DU145, time_aveRMST_12LNCaP, 'r');
% 
% RMST_scatter_ind12DU145 = 3*ones((length(time_aveRMST_12DU145)),1);
% f1 = scatter(RMST_scatter_ind12DU145, time_aveRMST_12DU145, 'r');
% 
% RMST_scatter_ind12PC3 = 4*ones((length(time_aveRMST_12PC3)),1);
% f1 = scatter(RMST_scatter_ind12PC3, time_aveRMST_12PC3, 'b');
%       
% % saveas('f_box_stiff_rmst.tiff');
% 
% figure(f_box_stiff_SE)
% hold on 
% boxplot([time_aveSE_12RV1,time_aveSE_12LNCaP ,time_aveSE_12DU145, time_aveSE_12PC3], [1*ones((length(time_aveSE_12RV1)),1); 2*ones((length(time_aveSE_12LNCaP)),1); ...
%           3*ones((length(time_aveSE_12DU145)),1); 4*ones((length(time_aveSE_12PC3)),1)], 'Positions',[1*ones((length(time_aveSE_12RV1)),1); 2*ones((length(time_aveSE_12LNCaP)),1); ...
%           3*ones((length(time_aveSE_12DU145)),1); 4*ones((length(time_aveSE_12PC3)),1)],  'Labels', {'22RV1 12 kPa', 'LNCaP', 'DU145', 'PC3 12 kPa'});
%       
% SE_scatter_ind12RV1 = 1*ones((length(time_aveSE_12RV1)),1);
% f1 = scatter(SE_scatter_ind12RV1, time_aveSE_12RV1, 'r');
% 
% SE_scatter_ind12LNCaP = 2*ones((length(time_aveSE_12LNCaP)),1);
% f1 = scatter(SE_scatter_ind12LNCaP, time_aveSE_12LNCaP, 'r');
% 
% SE_scatter_ind12DU145 = 3*ones((length(time_aveSE_12DU145)),1);
% f1 = scatter(SE_scatter_ind12DU145, time_aveSE_12DU145, 'r');
% 
% SE_scatter_ind12PC3 = 4*ones((length(time_aveRMST_12PC3)),1);
% f1 = scatter(SE_scatter_ind12PC3, time_aveSE_12PC3, 'b');
% 
% % saveas(gcf, 'f_box_stiff_SE.tiff');
% %% SOFT TIME PLOTS 
% f_soft_rmst = figure('Name', 'RMST vs. Time, 22RV1 & PC3 3kPa');
% ylabel('RMST (Pa)');
% xlabel('Time (Minutes)'); 
% 
% f_soft_SE = figure('Name', 'SE vs. Time 22RV1 vs. PC3 3kPa');
% ylabel('Strain  Energy (pJ)');
% xlabel('Time (Minutes)'); 
%   
% % Plot RMST and SE of 3 kPa 22RV1
% 
% % Identify Relevant Exps
% to_open = text_cats(RV1_3kPa_72_CSS_ind, 9);
% 
% % Identify tie interval
% time_interval = cell2mat(text_cats(RV1_3kPa_72_CSS_ind,10))/60;
% frame_add = cell2mat(text_cats(RV1_3kPa_72_CSS_ind,6));
% 
% % Identify bad positions
% series_omit_array = text_cats(RV1_3kPa_72_CSS_ind,7);
% 
% % Run through each position and plot 
% q=1;
% 
% for i = 1:length(to_open) 
%     % Open file and extract values rxc = pos x time
%   data_temp = open(['Mats_So_Far/' to_open{i}]);
%   rmst = data_temp.rmst';
%   strain_energy = data_temp.strain_energy'*10^12; % Adjust to pJ
%    average_displacement = data_temp.average_displacement';
%   
%   % Identify positions to omit and do so
%   if series_omit_array{i} ~=0 
%     try 
%     omit_indices = str2double(strsplit(series_omit_array{i}, ','));
%     catch
%     omit_indices = series_omit_array{i}; 
%     end 
%     rmst(omit_indices,:) = []; 
%     strain_energy(omit_indices, :) = [];
%     average_displacement(omit_indices,:) = [];
%   end     
% 
% 
%   % Identify the frame at which the treatment was applied 
%   [pos_num, time_num] = size(rmst);
%   if frame_add(i) ~= 0
%       time_num = frame_add(i) - 1;
%   end 
%   t = ([1:time_num]-1)*time_interval(i); 
% 
%   for j = 1:pos_num  
%   figure(f_soft_rmst)  
%   hold on
%   p1 = plot(t, rmst(j,1:length(t)), 'r'); 
%   time_aveRMST_3RV1(q) = mean(rmst(j, 1:length(t)));
% % text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
%   figure(f_soft_SE) 
%   hold on
%   q1 =  plot(t, strain_energy(j,1:length(t)), 'r'); 
%   time_aveSE_3RV1(q) = mean(strain_energy(j, 1:length(t)));
%   q = q+1;
%   end 
% end 
% 
% % Plot RMST and SE of 3 kPa PC3
% % Identify Relevant Exps
% to_open = text_cats(PC3_3kPa_72_CSS_ind , 9);
% 
% % Identify time interval
% time_interval = cell2mat(text_cats(PC3_3kPa_72_CSS_ind ,10))/60;
% 
% % Identify when the frame for adding 
% frame_add = cell2mat(text_cats(PC3_3kPa_72_CSS_ind ,6));
% 
% % Identify bad positions
% series_omit_array = text_cats(PC3_3kPa_72_CSS_ind ,7);
% 
% % Run through each position and plot 
% q=1;
% for i = 1:length(to_open) 
%     % Open file and extract values rxc = pos x time
%   data_temp = open(['Mats_So_Far/' to_open{i}]);
%   rmst = data_temp.rmst';
%   strain_energy = data_temp.strain_energy'*10^12;
%    average_displacement = data_temp.average_displacement';
%   
%   % Identify positions to omit and do so
%   if series_omit_array{i} ~=0 
%     try 
%     omit_indices = str2double(strsplit(series_omit_array{i}, ','));
%     catch
%     omit_indices = series_omit_array{i}; 
%     end 
%     rmst(omit_indices,:) = []; 
%     strain_energy(omit_indices, :) = [];
%     average_displacement(omit_indices,:) = [];
%   end     
% 
%   [pos_num, time_num] = size(rmst);
%   if frame_add(i) ~= 0
%       time_num = frame_add(i) - 1;
%   end 
%   t = ([1:time_num]-1)*time_interval(i); 
%   
%   for j = 1:pos_num  
%       figure(f_soft_rmst)  
%       hold on
%       p2 = plot(t, rmst(j,1:length(t)), 'b');
%       time_aveRMST_3PC3(q) = mean(rmst(j, 1:length(t)));
% % text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
%       figure(f_soft_SE)
%       hold on
%       q2 = plot(t, strain_energy(j,1:length(t)), 'b');
%       time_aveSE_3PC3(q) = mean(strain_energy(j, 1:length(t)));
%   q = q+1;
%   end 
% end 
%    
% % Plot RMST and SE of 3 kPa DU145 in CSS 
% % Identify Relevant Exps
% to_open = text_cats(DU145_3kPa_72_CSS_ind, 9);
% 
% % Identify time interval
% time_interval = cell2mat(text_cats(DU145_3kPa_72_CSS_ind ,10))/60;
% 
% % Identify when the frame for adding 
% frame_add = cell2mat(text_cats(DU145_3kPa_72_CSS_ind ,6));
% 
% % Identify bad positions
% series_omit_array = text_cats(DU145_3kPa_72_CSS_ind ,7);
% 
% % Run through each position and plot 
% q=1;
% for i = 1:length(to_open) 
%     % Open file and extract values rxc = pos x time
%   data_temp = open(['Mats_So_Far/' to_open{i}]);
%   rmst = data_temp.rmst';
%   strain_energy = data_temp.strain_energy'*10^12;
%    average_displacement = data_temp.average_displacement';
%   
%   % Identify positions to omit and do so
%   if series_omit_array{i} ~=0 
%     try 
%     omit_indices = str2double(strsplit(series_omit_array{i}, ','));
%     catch
%     omit_indices = series_omit_array{i}; 
%     end 
%     rmst(omit_indices,:) = []; 
%     strain_energy(omit_indices, :) = [];
%     average_displacement(omit_indices,:) = [];
%   end     
% 
%   [pos_num, time_num] = size(rmst);
%   if frame_add(i) ~= 0
%       time_num = frame_add(i) - 1;
%   end 
%   t = ([1:time_num]-1)*time_interval(i); 
%   
%   for j = 1:pos_num  
%       figure(f_soft_rmst)  
%       hold on
%       p2 = plot(t, rmst(j,1:length(t)), 'g');
%       time_aveRMST_3DU145(q) = mean(rmst(j, 1:length(t)));
% % text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
%       figure(f_soft_SE)
%       hold on
%       q2 = plot(t, strain_energy(j,1:length(t)), 'g');
%       time_aveSE_3DU145(q) = mean(strain_energy(j, 1:length(t)));
%   q = q+1;
%   end 
%   
% end 
% 
%   
% % Plot RMST and SE of 3 kPa DU145 in CSS 
% % Identify Relevant Exps
% to_open = text_cats(LNCaP_3kPa_72_CSS_ind, 9);
% 
% % Identify time interval
% time_interval = cell2mat(text_cats(LNCaP_3kPa_72_CSS_ind ,10))/60;
% 
% % Identify when the frame for adding 
% frame_add = cell2mat(text_cats(LNCaP_3kPa_72_CSS_ind ,6));
% 
% % Identify bad positions
% series_omit_array = text_cats(LNCaP_3kPa_72_CSS_ind ,7);
% 
% % Run through each position and plot 
% q=1;
% for i = 1:length(to_open) 
%     % Open file and extract values rxc = pos x time
%   data_temp = open(['Mats_So_Far/' to_open{i}]);
%   rmst = data_temp.rmst';
%   strain_energy = data_temp.strain_energy'*10^12;
%    average_displacement = data_temp.average_displacement';
%   
%   % Identify positions to omit and do so
%   if series_omit_array{i} ~=0 
%     try 
%     omit_indices = str2double(strsplit(series_omit_array{i}, ','));
%     catch
%     omit_indices = series_omit_array{i}; 
%     end 
%     rmst(omit_indices,:) = []; 
%     strain_energy(omit_indices, :) = [];
%     average_displacement(omit_indices,:) = [];
%   end     
% 
%   [pos_num, time_num] = size(rmst);
%   if frame_add(i) ~= 0
%       time_num = frame_add(i) - 1;
%   end 
%   t = ([1:time_num]-1)*time_interval(i); 
%   
%   for j = 1:pos_num  
%       figure(f_soft_rmst)  
%       hold on
%       p2 = plot(t, rmst(j,1:length(t)), 'm');
%       time_aveRMST_3LNCaP(q) = mean(rmst(j, 1:length(t)));
% % text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
%       figure(f_soft_SE)
%       hold on
%       q2 = plot(t, strain_energy(j,1:length(t)), 'm');
%       time_aveSE_3LNCaP(q) = mean(strain_energy(j, 1:length(t)));
%   q = q+1;
%   end 
%   
% end 
% 
% 
% hold off 
% figure(f_soft_rmst)
% hold on
% legend([p1,p2], '22RV1, 3 kPa', 'PC3, 3 kPa')
% saveas(gcf, 'soft_rmst.tiff');
% 
% figure(f_soft_SE)
% hold on
% legend([q1,q2], '22RV1, 3 kPa', 'PC3, 3 kPa')
% saveas(gcf, 'soft_SE.tiff');
% %% RMST and Strain Energy Boxplots Soft
% f_box_soft_rmst = figure('Name', 'RMST boxplot 22RV1 vs. PC3 3kPa');
% ylabel('RMST (Pa)');
% 
% f_box_soft_SE = figure('Name', 'SE boxplot 22RV1 vs. PC3 3kPa');
% ylabel('Strain  Energy (pJ)');
% 
% figure(f_box_soft_rmst)
% hold on 
% boxplot([time_aveRMST_3RV1,time_aveRMST_3LNCaP ,time_aveRMST_3DU145, time_aveRMST_3PC3], [1*ones((length(time_aveRMST_3RV1)),1); 2*ones((length(time_aveRMST_3LNCaP)),1); ...
%           3*ones((length(time_aveRMST_3DU145)),1); 4*ones((length(time_aveRMST_3PC3)),1)], 'Positions',[1*ones((length(time_aveRMST_3RV1)),1); 2*ones((length(time_aveRMST_3LNCaP)),1); ...
%           3*ones((length(time_aveRMST_3DU145)),1); 4*ones((length(time_aveRMST_3PC3)),1)],  'Labels', {'22RV1 3 kPa', 'LNCaP', 'DU145', 'PC3 3 kPa'});
%       
% RMST_scatter_ind3RV1 = 1*ones((length(time_aveRMST_3RV1)),1);
% f1 = scatter(RMST_scatter_ind3RV1, time_aveRMST_3RV1, 'r');
% 
% RMST_scatter_ind3LNCaP = 2*ones((length(time_aveRMST_3LNCaP)),1);
% f1 = scatter(RMST_scatter_ind3LNCaP, time_aveRMST_3LNCaP, 'r');
% 
% RMST_scatter_ind3DU145 = 3*ones((length(time_aveRMST_3DU145)),1);
% f1 = scatter(RMST_scatter_ind3DU145, time_aveRMST_3DU145, 'r');
% 
% RMST_scatter_ind3PC3 = 4*ones((length(time_aveRMST_3PC3)),1);
% f1 = scatter(RMST_scatter_ind3PC3, time_aveRMST_3PC3, 'b');
%       
% % saveas('f_box_soft_rmst.tiff');
% 
% figure(f_box_soft_SE)
% hold on 
% boxplot([time_aveSE_3RV1,time_aveSE_3LNCaP ,time_aveSE_3DU145, time_aveSE_3PC3], [1*ones((length(time_aveSE_3RV1)),1); 2*ones((length(time_aveSE_3LNCaP)),1); ...
%           3*ones((length(time_aveSE_3DU145)),1); 4*ones((length(time_aveSE_3PC3)),1)], 'Positions',[1*ones((length(time_aveSE_3RV1)),1); 2*ones((length(time_aveSE_3LNCaP)),1); ...
%           3*ones((length(time_aveSE_3DU145)),1); 4*ones((length(time_aveSE_3PC3)),1)],  'Labels', {'22RV1 3 kPa', 'LNCaP', 'DU145', 'PC3 3 kPa'});
%       
% SE_scatter_ind3RV1 = 1*ones((length(time_aveSE_3RV1)),1);
% f1 = scatter(SE_scatter_ind3RV1, time_aveSE_3RV1, 'r');
% 
% SE_scatter_ind3LNCaP = 2*ones((length(time_aveSE_3LNCaP)),1);
% f1 = scatter(SE_scatter_ind3LNCaP, time_aveSE_3LNCaP, 'r');
% 
% SE_scatter_ind3DU145 = 3*ones((length(time_aveSE_3DU145)),1);
% f1 = scatter(SE_scatter_ind3DU145, time_aveSE_3DU145, 'r');
% 
% SE_scatter_ind3PC3 = 4*ones((length(time_aveRMST_3PC3)),1);
% f1 = scatter(SE_scatter_ind3PC3, time_aveSE_3PC3, 'b');
% 
% % saveas(gcf, 'f_box_soft_SE.tiff');
% 
% 
% %%
% % RMST boxplots 
% 
% f_box_soft_rmst = figure('Name', 'RMST boxplot, 22RV1 vs. PC3 3 kPa');
% ylabel('RMST (Pa)');
% % xlabel('Time (Minutes)'); 
% 
% f_box_soft_SE = figure('Name', 'SE boxplot, 22RV1 vs. PC3 3 kPa');
% ylabel('Strain  Energy (pJ)');
% % xlabel('Time (Minutes)'); 
% 
% 
% figure(f_box_soft_rmst)
% hold on 
% boxplot([time_aveRMST_3RV1, time_aveRMST_3PC3], [1*ones((length(time_aveRMST_3RV1)),1);...
%           2*ones((length(time_aveRMST_3PC3)),1)], 'Positions',[1*ones((length(time_aveRMST_3RV1)),1);...
%           2*ones((length(time_aveRMST_3PC3)),1)],  'Labels', {'22RV1, 3 kPa', 'PC3, 3 kPa'});
%       
% RMST_scatter_ind3RV1 = 1*ones((length(time_aveRMST_3RV1)),1);
% f1 = scatter(RMST_scatter_ind3RV1, time_aveRMST_3RV1, 'r');
%   
% RMST_scatter_ind3PC3 = 2*ones((length(time_aveRMST_3PC3)),1);
% f1 = scatter(RMST_scatter_ind3PC3, time_aveRMST_3PC3, 'b');
% 
% saveas(gcf, 'f_box_soft_rmst.tiff')
%  
% figure(f_box_soft_SE)
% hold on 
% boxplot([time_aveSE_3RV1, time_aveSE_3PC3], [1*ones((length(time_aveSE_3RV1)),1);...
%           2*ones((length(time_aveSE_3PC3)),1)], 'Positions',[1*ones((length(time_aveSE_3RV1)),1);...
%           2*ones((length(time_aveSE_3PC3)),1)],  'Labels', {'22RV1, 3 kPa', 'PC3, 3 kPa'});
%       
% SE_scatter_ind3RV1 = 1*ones((length(time_aveSE_3RV1)),1);
% f1 = scatter(SE_scatter_ind3RV1, time_aveSE_3RV1, 'r');
%   
% SE_scatter_ind3PC3 = 2*ones((length(time_aveSE_3PC3)),1);
% f1 = scatter(SE_scatter_ind3PC3, time_aveSE_3PC3, 'b');
% 
% saveas(gcf, 'f_box_soft_SE.tiff')
%  
% %%  Stiffness Sensitivity Summary 
% f_box_rmst_softness_compare = figure('Name', 'CSS PC3/22RV1 RMST Stiffness Compare');
% hold on
% ylabel('RMST (Pa)');
% boxplot([time_aveRMST_3RV1, time_aveRMST_12RV1,time_aveRMST_3LNCaP, time_aveRMST_12LNCaP,   time_aveRMST_3DU145, time_aveRMST_12DU145 time_aveRMST_3PC3, time_aveRMST_12PC3], ...
%         [0.8*ones((length(time_aveRMST_3RV1)),1);... 
%           1.2*ones((length(time_aveRMST_12RV1)),1); 1.8*ones((length(time_aveRMST_3LNCaP)),1);...
%           2.2*ones((length(time_aveRMST_12LNCaP)),1);2.8*ones((length(time_aveRMST_3DU145)),1);...
%           3.2*ones((length(time_aveRMST_12DU145)),1); 3.8*ones((length(time_aveRMST_3PC3)),1);...
%           4.2*ones((length(time_aveRMST_12PC3)),1)],'Positions', [0.8*ones((length(time_aveRMST_3RV1)),1);... 
%           1.2*ones((length(time_aveRMST_12RV1)),1); 1.8*ones((length(time_aveRMST_3LNCaP)),1);...
%           2.2*ones((length(time_aveRMST_12LNCaP)),1);2.8*ones((length(time_aveRMST_3DU145)),1);...
%           3.2*ones((length(time_aveRMST_12DU145)),1); 3.8*ones((length(time_aveRMST_3PC3)),1);...
%           4.2*ones((length(time_aveRMST_12PC3)),1)], 'Labels', {'22RV1, 3 kPa', '22RV1, 12 kPa',...
%            'LNCaP, 3 kPa', 'LNCaP, 12 kPa', 'DU145, 3 kPa', 'DU145, 12 kPa', 'PC3, 3 kPa', '22RV1, 12 kPa'}, 'widths', 0.3);
% set(gca, 'FontSize', 15);
%       
% delete(findobj(gca,'type','text'))  % wipe out the boxplot labels that are problematic
% set(gca,'xtick',1:4,'xticklabel',{ '22RV1' 'LNCaP' 'DU145' 'PC3'})  % and write the labels wanted
% xlim([0.5,4.5])
% 
% RMST_scatter_ind3RV1 = 0.8*ones((length(time_aveRMST_3RV1)),1);
% f1 = scatter(RMST_scatter_ind3RV1, time_aveRMST_3RV1, 'ro');
% 
% RMST_scatter_ind12RV1 = 1.2*ones((length(time_aveRMST_12RV1)),1);
% f1 = scatter(RMST_scatter_ind12RV1, time_aveRMST_12RV1, 'rd');
% 
% RMST_scatter_ind3LNCaP = 1.8*ones((length(time_aveRMST_3LNCaP)),1);
% f1 = scatter(RMST_scatter_ind3LNCaP, time_aveRMST_3LNCaP, 'mo');
% 
% RMST_scatter_ind12LNCaP = 2.2*ones((length(time_aveRMST_12LNCaP)),1);
% f1 = scatter(RMST_scatter_ind12LNCaP, time_aveRMST_12LNCaP, 'md');
% 
% RMST_scatter_ind3DU145 = 2.8*ones((length(time_aveRMST_3DU145)),1);
% f1 = scatter(RMST_scatter_ind3DU145, time_aveRMST_3DU145, 'go');
% 
% RMST_scatter_ind12DU145 = 3.2*ones((length(time_aveRMST_12DU145)),1);
% f1 = scatter(RMST_scatter_ind12DU145, time_aveRMST_12DU145, 'gd');
% 
% 
% RMST_scatter_ind3PC3 = 3.8*ones((length(time_aveRMST_3PC3)),1);
% f1 = scatter(RMST_scatter_ind3PC3, time_aveRMST_3PC3, 'bo');
% 
% RMST_scatter_ind12PC3 = 4.2*ones((length(time_aveRMST_12PC3)),1);
% f1 = scatter(RMST_scatter_ind12PC3, time_aveRMST_12PC3, 'bd');
% 
% xlim([0.25,4.75]);
% mlabel1 = scatter(10,1, 'ko'); 
% mlabel2 = scatter(10,1, 'kd'); 
% 
% legend([mlabel1, mlabel2], '3 kPa', '12 kPa', 'Location', 'northwest')
% set(gca, 'FontSize', 15, 'LineWidth', 1);
% % ylim([0,16])
% ylim([0,500])
% 
% saveas(gcf, 'CSS_f_box_RMST_stiffness_compare', 'tiff')
%  
% f_box_SE_stiffness_compare = figure('Name','CSS PC3/22RV1 SE Stiffness Compare');
% hold on
% ylabel('Strain Energy (pJ)');
% boxplot([time_aveSE_3RV1, time_aveSE_12RV1,time_aveSE_3LNCaP, time_aveSE_12LNCaP,   time_aveSE_3DU145, time_aveSE_12DU145 time_aveSE_3PC3, time_aveSE_12PC3], ...
%         [0.8*ones((length(time_aveSE_3RV1)),1);... 
%           1.2*ones((length(time_aveSE_12RV1)),1); 1.8*ones((length(time_aveSE_3LNCaP)),1);...
%           2.2*ones((length(time_aveSE_12LNCaP)),1);2.8*ones((length(time_aveSE_3DU145)),1);...
%           3.2*ones((length(time_aveSE_12DU145)),1); 3.8*ones((length(time_aveSE_3PC3)),1);...
%           4.2*ones((length(time_aveSE_12PC3)),1)],'Positions', [0.8*ones((length(time_aveSE_3RV1)),1);... 
%           1.2*ones((length(time_aveSE_12RV1)),1); 1.8*ones((length(time_aveSE_3LNCaP)),1);...
%           2.2*ones((length(time_aveSE_12LNCaP)),1);2.8*ones((length(time_aveSE_3DU145)),1);...
%           3.2*ones((length(time_aveSE_12DU145)),1); 3.8*ones((length(time_aveSE_3PC3)),1);...
%           4.2*ones((length(time_aveSE_12PC3)),1)], 'Labels', {'22RV1, 3 kPa', '22RV1, 12 kPa',...
%            'LNCaP, 3 kPa', 'LNCaP, 12 kPa', 'DU145, 3 kPa', 'DU145, 12 kPa', 'PC3, 3 kPa', '22RV1, 12 kPa'}, 'widths', 0.3);
% set(gca, 'FontSize', 15);
%       
% delete(findobj(gca,'type','text'))  % wipe out the boxplot labels that are problematic
% set(gca,'xtick',1:4,'xticklabel',{ '22RV1' 'LNCaP' 'DU145' 'PC3'})  % and write the labels wanted
% xlim([0.5,4.5])
% 
% SE_scatter_ind3RV1 = 0.8*ones((length(time_aveSE_3RV1)),1);
% f1 = scatter(SE_scatter_ind3RV1, time_aveSE_3RV1, 'ro');
% 
% SE_scatter_ind12RV1 = 1.2*ones((length(time_aveSE_12RV1)),1);
% f1 = scatter(SE_scatter_ind12RV1, time_aveSE_12RV1, 'rd');
% 
% SE_scatter_ind3LNCaP = 1.8*ones((length(time_aveSE_3LNCaP)),1);
% f1 = scatter(SE_scatter_ind3LNCaP, time_aveSE_3LNCaP, 'mo');
% 
% SE_scatter_ind12LNCaP = 2.2*ones((length(time_aveSE_12LNCaP)),1);
% f1 = scatter(SE_scatter_ind12LNCaP, time_aveSE_12LNCaP, 'md');
% 
% SE_scatter_ind3DU145 = 2.8*ones((length(time_aveSE_3DU145)),1);
% f1 = scatter(SE_scatter_ind3DU145, time_aveSE_3DU145, 'go');
% 
% SE_scatter_ind12DU145 = 3.2*ones((length(time_aveSE_12DU145)),1);
% f1 = scatter(SE_scatter_ind12DU145, time_aveSE_12DU145, 'gd');
% 
% 
% SE_scatter_ind3PC3 = 3.8*ones((length(time_aveSE_3PC3)),1);
% f1 = scatter(SE_scatter_ind3PC3, time_aveSE_3PC3, 'bo');
% 
% SE_scatter_ind12PC3 = 4.2*ones((length(time_aveSE_12PC3)),1);
% f1 = scatter(SE_scatter_ind12PC3, time_aveSE_12PC3, 'bd');
% 
% xlim([0.25,4.75]);
% mlabel1 = scatter(10,1, 'ko'); 
% mlabel2 = scatter(10,1, 'kd'); 
% 
% 
% legend([mlabel1, mlabel2], '3 kPa', '12 kPa', 'Location', 'northwest')
% set(gca, 'FontSize', 15, 'LineWidth', 1);
% ylim([0,16])
% % ylim([0,500])
% saveas(gcf, 'CSS_f_box_SE_stiffness_compare', 'tiff')
% 
% %% OLD CONTENT 
% for i = 1
%     %% Old content
% 
% f_soft_rmst = figure('Name', 'RMST vs. Time, 22RV1 & PC3 3 kPa');
% ylabel('RMST (Pa)');
% xlabel('Time (Minutes)'); 
% 
% f_soft_SE = figure('Name', 'SE vs. Time, 22RV1 & PC3 3 kPa');
% ylabel('Strain  Energy (pJ)');
% xlabel('Time (Minutes)'); 
% 
% % Identify Relevant Exps
% to_open = text_cats(RV1_3kPa_72_FBS_ind, 9);
% 
% % Identify tie interval
% time_interval = cell2mat(text_cats(RV1_3kPa_72_CSS_ind,10))/60;
% frame_add = cell2mat(text_cats(RV1_3kPa_72_CSS_ind,6));
% 
% % Identify bad positions
% series_omit_array = text_cats(RV1_3kPa_72_CSS_ind,7);
% 
% % Run through each position and plot 
% q=1;
% 
% 
% for i = 1:length(to_open) 
%     % Open file and extract values rxc = pos x time
%   data_temp = open(['Mats_So_Far/' to_open{i}]);
%   rmst = data_temp.rmst';
%   strain_energy = data_temp.strain_energy'*10^12;
%   
%   % Identify positions to omit and do so
%   if series_omit_array{i} ~=0 
%     omit_indices = str2double(strsplit(series_omit_array{i}, ','));   
%     rmst(omit_indices,:) = []; 
%     strain_energy(omit_indices, :) = [];
%   end     
% 
%   [pos_num, time_num] = size(rmst);
%   if frame_add(i) ~= 0
%       time_num = frame_add(i) - 1;
%   end 
%   t = ([1:time_num]-1)*time_interval(i); 
%   
%   for j = 1:pos_num  
%   figure(f_soft_rmst)  
%   hold on
%   p1 = plot(t, rmst(j,1:length(t)), 'r'); 
%   time_aveRMST_3RV1(q) = mean(rmst(j, 1:length(t)));
% % text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
%   figure(f_soft_SE) 
%   hold on
%   q1 =  plot(t, strain_energy(j,1:length(t)), 'r'); 
% 
%   time_aveSE_3RV1(q) = mean(strain_energy(j, 1:length(t)));
%   q = q+1;
%   end 
% end 
% 
% % Identify Relevant Exps
% to_open = text_cats(PC3_3kPa_72_CSS_ind , 9);
% 
% % Identify tie interval
% time_interval = cell2mat(text_cats(PC3_3kPa_72_CSS_ind ,10))/60;
% frame_add = cell2mat(text_cats(PC3_3kPa_72_CSS_ind ,6));
% 
% % Identify bad positions
% series_omit_array = text_cats(PC3_3kPa_72_CSS_ind ,7);
% 
% % Run through each position and plot 
% q=1;
% for i = 1:length(to_open) 
%     % Open file and extract values rxc = pos x time
%   data_temp = open(['Mats_So_Far/' to_open{i}]);
%   rmst = data_temp.rmst';
%   strain_energy = data_temp.strain_energy'*10^12;
%   
%   % Identify positions to omit and do so
%   if series_omit_array{i} ~=0
%     omit_indices = str2double(strsplit(num2str(series_omit_array{i})', ','));   
%     rmst(omit_indices,:) = []; 
%     strain_energy(omit_indices, :) = [];
%   end     
% 
%   [pos_num, time_num] = size(rmst);
%   if frame_add(i) ~= 0
%       time_num = frame_add(i) - 1;
%   end 
%   t = ([1:time_num]-1)*time_interval(i); 
%   
%   for j = 1:pos_num  
%       figure(f_soft_rmst)  
%       hold on
%       p2 = plot(t, rmst(j,1:length(t)), 'b');
%       time_aveRMST_3PC3(q) = mean(rmst(j, 1:length(t)));
% % text(t(end),StrainEnergy(i,end),['Cell ' num2str(i)]);
%       figure(f_soft_SE)
%       hold on
%       q2 = plot(t, strain_energy(j,1:length(t)), 'b');
%       time_aveSE_3PC3(q) = mean(strain_energy(j, 1:length(t)));
%   q = q+1;
%   end 
%   
% end 
% 
% figure(f_soft_rmst)
% hold on
% legend([p1,p2], '22RV1, 3 kPa', 'PC3, 3 kPa')
% saveas(gcf, 'f_soft_rmst.tiff')
% 
% figure(f_soft_SE)
% hold on
% legend([q1,q2], '22RV1, 3 kPa', 'PC3, 3 kPa')
% saveas(gcf, 'f_soft_SE.tiff')
% end 
% 
% % Media comparisons 
% for i = 1
% old_CSS = open('Mats_So_Far/BR_PC3_CSS_3_72.mat'); 
% new_CSS = open('Mats_So_Far/BS_PC3_12_CSS0.7CSSFRAME6_72.mat');
% 
% time_interval_old = 1914.33/60;
% Frame_treatment_old = 3;
% 
% time_interval_new = 1704.11/60; 
% Frame_treatment_new = 6;
% 
% old_rmst = old_CSS.rmst'; 
% new_rmst = new_CSS.rmst'; 
% 
% old_SE = old_CSS.strain_energy'*10^12 
% new_SE = new_CSS.strain_energy'*10^12 
% 
% omit_indices_old = 1; 
% omit_indices_new = 2;
% 
% old_SE(omit_indices_old, :) = [];
% new_SE(omit_indices_new, :) = [];
% 
% old_rmst(omit_indices_old, :) = [];
% new_rmst(omit_indices_new, :) = [];
% 
%  [pos_num_old, time_num_old] = size(old_rmst); 
%  [pos_num_new, time_num_new] = size(new_rmst);
%  t_old =([1:time_num_old]-1)*time_interval_old;
%  t_new =([1:time_num_new]-1)*time_interval_new;
%  
%  f_rmst_CSS_age = figure('Name', 'CSS Age Effect - RMST vs. Time');
%  hold on
%       
%  % RMST old and new
%   for j = 1:pos_num_old  
%       p1 = plot(t_old, old_rmst(j,1:length(t_old)), 'r');
%       %time_aveRMST_3PC3(q) = mean(rmst(j, 1:length(t)));
%   end 
%       x1 = xline(t_new(Frame_treatment_old), 'r--');
%       
%  q = 1;
%   for k = 1:pos_num_new
%       hold on
%       p2 = plot(t_new, new_rmst(k,1:length(t_new)), 'g');
%       time_aveRMST_CSSPC3new(q) = mean(new_rmst(k, 1:length(t_new)));
%   q = q+1;
%   end 
% x2 = xline(t_new(Frame_treatment_new-1), 'g--');  
% 
% legend([p1,x1,p2,x2], 'PC3 3kPa || (Old) CSS', 'FBS Added to Old CSS', 'PC3 12 kPa || (New) CSS', 'FBS Added to New CSS', 'Location', 'East')
% xlabel('Time (minutes)'); 
% ylabel('RMST (Pa)');
% 
% saveas(gcf, 'f_RMST_CSS_age.tiff'); 
% 
%  f_SE_CSS_age = figure('Name', 'CSS Age Effect - SE vs. Time');
%  hold on
%       
%  % RMST old and new
%   for j = 1:pos_num_old  
%       p1 = plot(t_old, old_SE(j,1:length(t_old)), 'r');
%       %time_aveSE_3PC3(q) = mean(rmst(j, 1:length(t)));
%   end 
%       x1 = xline(t_new(Frame_treatment_old), 'r--');
%       
%  q = 1;
%   for k = 1:pos_num_new
%       hold on
%       p2 = plot(t_new, new_SE(k,1:length(t_new)), 'g');
%       time_aveSE_CSSPC3new(q) = mean(new_rmst(k, 1:length(t_new)));
%   q = q+1;
%   end 
% x2 = xline(t_new(Frame_treatment_new-1), 'g--');  
% 
% legend([p1,x1,p2,x2], 'PC3 3kPa || Old CSS', 'FBS Added to Old CSS', 'PC3 12 kPa || New CSS', 'FBS Added to New CSS')
% xlabel('Time (minutes)'); 
% ylabel('Strain Energy (pJ)');
% saveas(gcf, 'f_SE_CSS_age.tiff'); 
% 
% 
% %% 72 hour casodex treatment Time course 
% % 3 kPa multi-treatment 
% 
% CSST_12 = open('Mats_So_Far/BQ_22RV1_FBS_12_72_CDx3.mat'); 
% CSST_3 = open('Mats_So_Far/BQ_22RV1_FBS_3_72.mat');
% % CSS_treat = open(
% 
% time_interval = 1588.24/60;
% % Frame_treatment_old = 3;
% 
% % time_interval_new = 1704.11/60; 
% % Frame_treatment_new = 6;
% 
% CSST_12_rmst = CSST_12.rmst';
% CSST_12_SE = CSST_12.strain_energy'*10^12; 
% % Timepoint 7 is out of focu --> Take averages: 
% CSST_12_rmst(:,7) = (CSST_12_rmst(:,8)+ CSST_12_rmst(:,6))/2;
% CSST_12_SE(:,7) = (CSST_12_SE(:,8)+ CSST_12_SE(:,6))/2;
% 
% CSST_3_rmst = CSST_3.rmst';
% CSST_3_SE = CSST_3.strain_energy'*10^12; 
% 
% [pos_num_CSST_12, time_num_CSST_12] = size(CSST_12_rmst); 
% [pos_num_CSST_3, time_num_CSST_3] = size(CSST_3_rmst); 
% 
% t_CSST_12 =([1:time_num_CSST_12]-1)*time_interval;
% t_CSST_3 =([1:time_num_CSST_3]-1)*time_interval;
% 
% % ARpos3kPa_ind = [1:6] ;
% f_CSST3_rmst = figure('Name', '22RV1 CSS 3kPa - RMST vs. Time');
% hold on 
% xlabel('Time (minutes)');
% ylabel('RMST (Pa)');
% 
% f_CSST3_SE = figure('Name', '22RV1 CSS 3kPa - SE vs. Time');
% hold on 
% xlabel('Time (minutes)');
% ylabel('Strain Energy (pJ)');
% 
% f_CSST12_rmst = figure('Name', '22RV1 CSS 12kPa - RMST vs. Time');
% hold on 
% xlabel('Time (minutes)','FontWeight','bold','FontName','Arial');
% ylabel('RMST (Pa)');
% 
% f_CSST12_SE = figure('Name', '22RV1 CSS 12kPa - SE vs. Time');
% hold on 
% xlabel('Time (minutes)','FontWeight','bold','FontName','Arial');
% ylabel('Strain Energy (pJ)');
% 
% % CSST_12
% q = 1; 
% for i = 1:pos_num_CSST_12
% figure(f_CSST12_rmst);
% hold on
% p1 = plot(t_CSST_12, CSST_12_rmst(i,:), 'b'); 
% time_aveRMST_CSST_12(q) = mean(CSST_12_rmst(i));
% x1 = xline(t_CSST_12(3), 'g-'); 
% 
% figure(f_CSST12_SE);
% hold on 
% q1 = plot(t_CSST_12, CSST_12_SE(i,:),'b'); 
%  time_aveSE_CSST_12(q) = mean(CSST_12_SE(i));
% x1 = xline(t_CSST_12(3), 'g-'); 
% 
% q = q+1;
% end 
% 
% % CSST3
% q = 1; 
% for j = 1:pos_num_CSST_3
% figure(f_CSST3_rmst);
% hold on
% p2 = plot(t_CSST_3, CSST_3_rmst(j,:), 'r'); 
% time_aveCSST_3_FBS(q) = mean(CSST_3_rmst(j));
% x1 = xline(t_CSST_3(3), 'g-'); 
% 
% figure(f_CSST3_SE);
% hold on 
% q2 = plot(t_CSST_3, CSST_3_SE(j,:),'r'); 
% time_aveSE_CSST_3(q) = mean(CSST_3_SE(j));
% x1 = xline(t_CSST_3(3), 'g-'); 
% q = q+1;
% end 
% 
% figure(f_CSST12_rmst);
% hold on 
% ylabel('RMST (Pa)','FontWeight','bold','FontName','Arial');
% xlabel('Time (minutes)','FontWeight','bold','FontName','Arial');
% legend([p1,x1], 'FBS + CDX, 22RV1 12 kPa', 'CDX Added');
% saveas(gcf,'RMST12 CDXTreat.tiff')
% 
% figure(f_CSST3_rmst);
% hold on 
% ylabel('RMST (Pa)');
% xlabel('Time (minutes)');
% legend([p2,x1], 'FBS + CDX, 22RV1 3 kPa', 'CDX Added');
% saveas(gcf,'RMST3 CDXTreat.tiff')
% 
% figure(f_CSST12_SE);
% hold on 
% ylabel('Strain Energy (pJ)');
% xlabel('Time (minutes)');
% legend([q1,x1], 'FBS + CDX , 22RV1 12 kPa', 'CDX Added');
% hold off 
% saveas(gcf,'SE12 CSSTreat.tiff')
% 
% figure(f_CSST3_SE);
% hold on 
% ylabel('Strain Energy (pJ)');
% xlabel('Time (minutes)');
% legend([q2,x1], 'FBS + CDX, 22RV1 3 kPa', 'CDX Added');
% hold off 
% saveas(gcf,'SE3 CSSTreat.tiff')
% %%
% f_SE_CSSvFBS22RV1_12 = figure('Name', '22RV1 CSSvFBS 12kPa - SE Boxplot');
% hold on 
% ylabel('Strain Energy(Pa)');
% %  check = StrainEnergy;
% %  check = ~check == 0;
% % StrainEnergy = StrainEnergy(check);
% boxplot([time_aveSE_CSS'; time_aveSE_FBS'], [1*ones((length(time_aveSE_CSS)),1);...
%          2*ones((length(time_aveSE_FBS)),1)], 'Positions',[1*ones((length(time_aveSE_CSS)),1);...
%          2*ones((length(time_aveSE_FBS)),1)],  'Labels', {'CSS, 22RV1 12 kPa', 'FBS, 22RV1 12 kPa'});
%     
% SE_scatter_indCSS = 1*ones((length(time_aveSE_CSS)),1);
% f1 = scatter(SE_scatter_indCSS, time_aveSE_CSS, 'b');
% 
% SE_scatter_indFBS = 2*ones((length(time_aveSE_FBS)),1);
% f1 = scatter(SE_scatter_indFBS, time_aveSE_FBS, 'r');
% 
% saveas(gcf,'SE bx 22RV1 CSSvFBS.tiff');
% 
% figure('Name', '22RV1 CSSvFBS 12kPa - RMST Boxplot');
% hold on 
% ylabel('RMST(Pa)');
% %  check = StrainEnergy;
% %  check = ~check == 0;
% % StrainEnergy = StrainEnergy(check);
% boxplot([time_aveRMST_CSS'; time_aveRMST_FBS'], [1*ones((length(time_aveRMST_CSS)),1);...
%          2*ones((length(time_aveRMST_FBS)),1)], 'Positions',[1*ones((length(time_aveRMST_CSS)),1);...
%          2*ones((length(time_aveRMST_FBS)),1)],  'Labels', {'CSS, 22RV1 12 kPa', 'FBS 22RV1 12 kPa',});
%     
% RMST_scatter_indCSS = 1*ones((length(time_aveRMST_CSS)),1);
% f1 = scatter(RMST_scatter_indCSS, time_aveRMST_CSS, 'b');
% 
% RMST_scatter_indFBS = 2*ones((length(time_aveRMST_FBS)),1);
% f1 = scatter(RMST_scatter_indFBS, time_aveRMST_FBS, 'r');
% 
% saveas(gcf,'RMST box 22RV1 CSSvFBS.tiff');
% 
% 
% end 
