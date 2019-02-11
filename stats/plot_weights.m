%% plot weights from prediction model

clear all
close all

%% define paths

TR = 2.42;
addpath ~/Software/BCT
addpath /usr/local/freesurfer/matlab/
addpath /usr/local/freesurfer/fsfast/toolbox/
addpath /home/benjamin.garzon/Software/FSLnets_pack/FSLNets/

workdir = '/home/benjamin.garzon/Data/NSPN/ICA_ME/nspn_ME/ica200.gica';
coefs_dir = '/home/benjamin.garzon/Data/NSPN/module_stats/mediation_data_2k_7mods_100_decAc/';
coefs_file = fullfile(coefs_dir, 'decAc_coefs_mat.txt');

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

%% represent coefs
coefs = load(coefs_file);
showN = sum(coefs(:)~=0)/2;
showN = 9;
[netmat] = nets_edgepics(ts, group_maps, abs(coefs), coefs, showN);fig = gcf;

%%
set(gcf, 'Position',  [100, 100, 2000, 1600])
set(gcf, 'Position',  [100, 100, 1000, 1300])
print(fullfile(coefs_dir, 'decAc_connections'),'-dpng','-r500')

figure
hist(coefs(tril(coefs)~=0), 100)
%saveas(gcf, 'results/decAc_connections.png')

%coefs_full = zeros(ncomp,ncomp);
%coefs_full(ts.DD, ts.DD) = coefs;
%nets_netweb(Znet1, coefs_full, ts.DD, group_maps,'netweb');

