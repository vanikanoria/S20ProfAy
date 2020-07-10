function find_sync_score(mh1)
% %mh1 should have dimensions = number of cells x number of timesteps
% 
% double post_sync (sim_data& sd, con_levels& cl, int con, int start, int end) {
% 	double comp_cell[end - start + 1];
% 	double cur_cell[end - start + 1];
% 	
% 	int middle_cell = (sd.height / 2) * sd.width_total + (sd.width_current / 2);
% 	
% 	for (int j = start; j < end; j++) {
% 		comp_cell[j - start] = cl.cons[con][j][middle_cell];
% 	}
% 	
% 	double pearson_sum = 0;
% 	for (int x = 0; x < sd.height; x++) {
% 		for (int y = 0; y < sd.width_initial; y++) {
% 			int cell = x * sd.width_total + y;
% 			
% 			if (cell != middle_cell) {
% 				for (int j = start; j < end; j++) {
% 					cur_cell[j - start] = cl.cons[con][j][cell];
% 				}
% 				pearson_sum += pearson_correlation(comp_cell, cur_cell, 0, end - start);
% 			}
% 		}
% 	}
% 	
% 	return pearson_sum / ((sd.height * sd.width_initial) - 1);
% }
% 
% 
% end