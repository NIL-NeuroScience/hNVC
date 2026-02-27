function data = f_loadDataset()

data_dir = utils.f_getPackagePath;
data_dir = fullfile(data_dir,'datasets','data.mat');

load_data = load(data_dir);

subjects = fieldnames(load_data);
N = numel(subjects);

data = struct;

for i = 1:N
    data.(subjects{i}) = struct;

    MR_labels = {'v1', 'v3', 'v4', 'v5', 'a1', ...
        'sommot', 'amygdata', 'caudate', 'pallidum', 'thalamus'};

    data.(subjects{i}).MR_labels = MR_labels;
    
    data.(subjects{i}).BOLD = NaN(numel( ...
        load_data.(subjects{i}).(MR_labels{1})), ...
        numel(MR_labels));

    for c = 1:numel(MR_labels)
        data.(subjects{i}).BOLD(:,c) = ...
            load_data.(subjects{i}).(MR_labels{c});
    end
end

% %%
% tmp = load(data_dir,subject);
% tmp = tmp.(subject);
% 
% data = struct;
% 
% % gather BOLD labels
% data.MR_labels = {'v1', 'v3', 'v4', 'v5', 'a1', 'sommot', 'amygdata', 'caudate', 'pallidum', 'thalamus'};
% 
% % compile BOLD data
% data.BOLD = NaN(numel(tmp.(data.MR_labels{1})),numel(data.MR_labels));
% for i = 1:numel(data.MR_labels); data.BOLD(:,i) = tmp.(data.MR_labels{i});
% 
% data.pupil = tmp.eye';
% 
% data.EEG_labels = {tmp.eeg_chnlocs.labels};
% data.EEG = tmp.EEG';
% 
% data.settings.MR_fs = 1 / tmp.fMRI_fs;
% data.settings.pupil_fs = 1 / tmp.eye_fs;
% data.settings.EEG_fs = tmp.EEG_fs;

end