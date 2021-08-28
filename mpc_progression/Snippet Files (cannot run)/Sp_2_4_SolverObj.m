opt_prob = struct('f', objfun, 'x', w, 'g', g, 'p', P);
config = struct;
config.print_time = 0;
config.ipopt.print_level = 0;  %0-3 min-max print
config.ipopt.warm_start_init_point = 'yes';

solver = nlpsol('ctrl_solver', 'ipopt', opt_prob, config);
