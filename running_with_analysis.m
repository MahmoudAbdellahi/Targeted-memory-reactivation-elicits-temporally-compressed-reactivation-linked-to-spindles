analysis_types = {'behavioural_analyses', 'classification', 'erp_tf', 'high_spindle_pw', 'low_spindle_pw', 'temporal_compression'};

for i = 1:length(analysis_types)
    code_copy_jitt(analysis_types{i});
end
