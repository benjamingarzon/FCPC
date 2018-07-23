%% plot weights from prediction model

clear all
close all

%% define paths

TR = 2.42;
addpath ~/Software/BCT
addpath /usr/local/freesurfer/matlab/
addpath /usr/local/freesurfer/fsfast/toolbox/
addpath /home/benjamin.garzon/Software/FSLnets_pack/FSLNets/
addpath /home/benjamin.garzon/Software/FSLnets_pack/NSPN/
workdir = '/home/benjamin.garzon/Data/NSPN/ICA_ME/nspn_ME/ica200.gica';
%atlas = '/home/benjamin.garzon/Data/NSPN/templates/MD.nii.gz';
bad_nets_file = 'bad_nets200.txt'
ncomp = 200;

cd(workdir)
group_maps=fullfile(workdir, 'groupmelodic.ica/melodic_IC');
ts_dir=fullfile(workdir, 'grot');


%% load data
ts=nets_load(ts_dir,TR,1);
bad_nets = load(bad_nets_file)';
ts.DD=setdiff([1:ncomp], bad_nets + 1);

ts=nets_tsclean(ts,1);

%% compute adjacencies and save

%
% netmats1 = nets_netmats(ts,1,'corr');       % full correlation (normalised covariances)
% netmats5 = nets_netmats(ts,1,'ridgep');     % Ridge Regression partial, with rho=0.01
%
% [Znet1,Mnet1]=nets_groupmean(netmats1,0); % full corr
% [Znet5,Mnet5]=nets_groupmean(netmats5,0); % ridge
%
% subjects = load(fullfile(workdir, 'subjects.txt'));
% netmats1_table = reshape_nets(netmats1, subjects);
% dlmwrite(fullfile(workdir, 'netmats_full.csv'), netmats1_table);
% netmats5_table = reshape_nets(netmats5, subjects);
% dlmwrite(fullfile(workdir, 'netmats_ridge.csv'), netmats5_table);
%

%% represent coefs

coefs = load('results/decAc_coefs.mat.txt');
showN = sum(coefs(:)~=0)/2;
[netmat] = nets_edgepics(ts, group_maps, abs(coefs), coefs, showN);
saveas(gcf, 'results/decAc_connections.png')

%coefs_full = zeros(ncomp,ncomp);
%coefs_full(ts.DD, ts.DD) = coefs;
%nets_netweb(Znet1, coefs_full, ts.DD, group_maps,'netweb');

