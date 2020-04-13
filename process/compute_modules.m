%% compute modularity and graph metrics that can be used for predictive modelling

clear all
close all

%% define parameters
N_gammas =  30; % number of gamma values to try
MINMODSIZE = 0; % min size acceptable for a modules (smaller modules will be removed)
REPS = 2000; % repetitions of modularity algorithm
REPSCONS = 2000; % repetitions of 
TAU = 0.2; % min strength for consensus matrix
TRYGAMMAS = 1; % 1 = try different gamma values to find the best partition
THR = 5; % thershold for ICA maps (don't show voxels with value below this)
gamma = 1.3; % fixed gamma value

%% define paths

TR = 2.42;
addpath ~/Software/BCT
addpath /usr/local/freesurfer/matlab/
addpath /usr/local/freesurfer/fsfast/toolbox/
addpath /home/benjamin.garzon/Software/FSLnets_pack/FSLNets/
addpath /home/benjamin.garzon/Software/FSLnets_pack/NSPN/
workdir = '/home/benjamin.garzon/Data/NSPN/ICA_ME/nspn_ME/ica200.gica';

bad_nets_file = 'bad_nets200.txt'
ncomp = 200;

cd(workdir)
group_maps=fullfile(workdir, 'groupmelodic.ica/melodic_IC');
ts_dir=fullfile(workdir, 'grot');


%% load data
ts=nets_load(ts_dir,TR,1);
bad_nets = load(bad_nets_file)';
ts.DD=setdiff([1:ncomp], bad_nets + 1);

%% compute adjacencies and save
ts=nets_tsclean(ts,1);

netmats1 = nets_netmats(ts,1,'corr');       % full correlation (normalised covariances)
netmats5 = nets_netmats(ts,1,'ridgep');     % Ridge Regression partial, with rho=0.01

[Znet1,Mnet1]=nets_groupmean(netmats1,0); % full corr
[Znet5,Mnet5]=nets_groupmean(netmats5,0); % ridge

subjects = load(fullfile(workdir, 'subjects.txt'));
netmats1_table = reshape_nets(netmats1, subjects);
dlmwrite(fullfile(workdir, 'netmats_full.csv'), netmats1_table);
netmats5_table = reshape_nets(netmats5, subjects);
dlmwrite(fullfile(workdir, 'netmats_ridge.csv'), netmats5_table);

%% compute modularity

adj = Znet5;
hist(adj(:), 100)
adj(adj < 0) = 0;

gammas = linspace(1, 4, N_gammas);
consensus = adj*0;

allM = zeros(size(adj,1), REPS);

if TRYGAMMAS
    % try different gammas
    for i = 1:N_gammas
        display(i)
        allM = zeros(size(adj,1), REPS);
        gamma = gammas(i);
        for j = 1 : REPS
            [M Q(i)] = community_louvain(adj, gamma);
            allM(:,j ) = M;
        end
        
        D = agreement(allM)/REPS;
        C = consensus_und(D, TAU, REPSCONS);
        maxM(i) = numel(unique(M));
        allC(:,i) = C;
        maxC(i) = numel(unique(C));
        freqs = tabulate(C);
        minC(i) = min(freqs(:,2));
    end
    
    NMI = zeros(N_gammas);
    VIn = zeros(N_gammas);
    
    % compute normalized mutual information (symmetric uncertainty)
    for i = 1:N_gammas
        for j = 1:N_gammas
            [VIn(i, j) NMI(i, j)] = partition_distance(allC(:, i), allC(:, j));
        end
    end
    
    meanNMI = (sum(NMI) - 1)/(N_gammas - 1);
    meanVIn = (sum(VIn) - 1)/(N_gammas - 1);
    minNMI = min(NMI);
    
    [m, index] = max(meanNMI);
    %[m, index] = min(meanVIn);
    
    % plot NMI and Vin    
    figure
    subplot(1,4,1)
    plot(gammas, maxM, 'b.-')
    hold on
    plot(gammas, maxC, 'b.-', 'LineWidth', 2)
    ylabel('N clusters')
    subplot(1,4,2)
    plot(gammas, meanNMI, '.-', 'LineWidth', 2)
    hold on
    plot(gammas, minNMI, 'r.-')
    plot(gammas, 4*meanVIn, 'g.-')
    ylabel('Mean NMI')
    
    subplot(1,4,3)
    imagesc(NMI)
    subplot(1,4,4)
    plot(gammas, minC, '.-')
    ylabel('Min Cluster Size')
    
    partition_strict = allC(:, index);
    gamma = gammas(index);
else
    
    % fix value for gamma and compute modularity
    for j = 1 : REPS
        [M Q] = community_louvain(adj, gamma);
        
        allM(:,j ) = M;
    end
    
    D = agreement(allM)/REPS;
    C = consensus_und(D, TAU, REPSCONS);
    
    partition_strict = C;
end


%% output images for each module

partition = partition_strict*0;
freqs = tabulate(partition_strict);
display(freqs)

plot(freqs(:))

%% remove modules with less than MINNODSIZE nodes
[x, sortind] = sort(freqs(:,2), 'descend');
for i = freqs(freqs(:, 2) < MINMODSIZE)'
    display(sum(partition_strict==i))
    partition_strict(partition_strict==i) = 0;
end

j = 1;
for i = freqs(sortind, 1)'
    partition(partition_strict == i) = j;
    j = j + 1;
end

%%
mri = MRIread(group_maps);
parc = mri.vol(:,:,:,ts.DD) > THR;
nvols = size(parc, 4);
s = size(parc);
s(4) = max(partition);
module_parc = zeros(s);
for i = 1:max(partition)
    module_parc(:,:,:,i) = mean(parc(:,:,:,partition == i), 4);
end

mri.vol = module_parc;

MRIwrite(mri, 'modules_optimal4.nii.gz')
dlmwrite('modules_optimal4.txt', partition)

%MRIwrite(mri, 'modules_1.3.nii.gz')
%dlmwrite('modules_1.3.txt', partition)
