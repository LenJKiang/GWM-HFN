function [results] = calculate_seguin_metrics_three_sessions(data1, data2, data3, network_type)
sim_matrix_12 = calculate_similarity_matrix_seguin(data1, data2);
sim_matrix_13 = calculate_similarity_matrix_seguin(data1, data3);
sim_matrix_23 = calculate_similarity_matrix_seguin(data2, data3);

metrics_12 = calculate_seguin_metrics(sim_matrix_12, network_type, 'BNU3_Session1-2');
metrics_13 = calculate_seguin_metrics(sim_matrix_13, network_type, 'BNU3_Session1-3');
metrics_23 = calculate_seguin_metrics(sim_matrix_23, network_type, 'BNU3_Session2-3');

all_intra = [metrics_12.intra_similarities; metrics_13.intra_similarities; metrics_23.intra_similarities];
all_inter = [metrics_12.inter_similarities; metrics_13.inter_similarities; metrics_23.inter_similarities];

overall_reliability = mean(all_intra);
overall_uniformity = mean(all_inter);
s_intra_overall = std(all_intra);
s_inter_overall = std(all_inter);

n_intra_overall = length(all_intra);
n_inter_overall = length(all_inter);
s_pooled_overall = sqrt(((n_intra_overall-1)*s_intra_overall^2 + (n_inter_overall-1)*s_inter_overall^2) / (n_intra_overall + n_inter_overall - 2));
overall_identifiability = (overall_reliability - overall_uniformity) / s_pooled_overall;

results = struct();
results.network_type = network_type;
results.dataset_name = 'BNU-3';
results.n_sessions = 3;
results.reliability = overall_reliability;
results.uniformity = overall_uniformity;
results.identifiability = overall_identifiability;
results.s_pooled = s_pooled_overall;
results.session_pairs.metrics_12 = metrics_12;
results.session_pairs.metrics_13 = metrics_13;
results.session_pairs.metrics_23 = metrics_23;
results.intra_similarities = all_intra;
results.inter_similarities = all_inter;
end