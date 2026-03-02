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

    data.(subjects{i}).pupil = load_data.(subjects{i}).eye';
    data.(subjects{i}).EEG_labels = {load_data.(subjects{i}).eeg_chnlocs.labels};
    data.(subjects{i}).EEG = load_data.(subjects{i}).EEG';

    data.(subjects{i}).MR_fs = 1 / load_data.(subjects{i}).fMRI_fs;
    data.(subjects{i}).pupil_fs = 1 / load_data.(subjects{i}).eye_fs;
    data.(subjects{i}).EEG_fs = load_data.(subjects{i}).EEG_fs;

end

end