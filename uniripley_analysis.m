% Univariate Ripley's K estimator (second-order function)
function calc_k = uniripley_analysis(x1, y1, points, xmin, ymin, xmax, ...
                                        ymax, xcol, ycol, sheet_name, ...
                                        sort, label, t, max_step, ...
                                        interval, n_simul, limit, ...
                                        k_method, bygroups, group, ...
                                        max_group, strow, Optional)
    testname = "Univariate Ripley's K";
    
    bins = AsymUp(max_step / t) + 1;
    if ((bins * t) > (max_step + t))
        bins = bins - 1; % Make sure don't exceed max step ....
    end
    
    % iniitlize all matrices to be zeros 
    kt_trans = zeros(bins);
    kt_theo = zeros(bins);
    hi_trans = zeros(bins);
    lo_trans = zeros(bins);

    if sizeWarn == true
        checkscale(max_step, xmin, xmax, ymin, ymax);
    end

    if sort == true && bygroups == false 
         sort_points(x1, y1, points, xcol, ycol, 1);
    elseif sort == true && bygroups == true 
         sort_points(x1, y1, points, xcol, 1, strow);
    end

    if label == true   
        Label_Points(xcol, points);
    end

    area = (xmax - xmin) * (ymax - ymin);

    % Empirical estimation of k(t) first
    pctDone = 0;
    
    if k_method == 1 
        calc_k(k_t, x1, y1, t, max_step, bins, points, area)
    elseif k_method == 2   
        calc_k_ew(k_t, x1, y1, t, max_step, bins, points, xmin, ymin, xmax, ymax, area)
    elseif k_method == 3   
        calc_k_toroid(k_t, x1, y1, t, max_step, bins, points, xmax, ymax, area)
    elseif k_method == 4   
        calc_k_buffer(k_t, x1, y1, t, max_step, bins, points, xmin, ymin, xmax, ymax, area);
    end
    
    % Calculate theoretical values && l(t) transformation
    i = 1;
    t_incr = 0;
    while (i <= bins)
        kt_trans(i) = Sqr(k_t(i) / (Pi)) - t_incr;
        kt_theo(i) = Pi * t_incr ^ 2;
        i = i + 1; t_incr = t_incr + t;
    end
    
    % Now the Monte Carlo simualtions
    if (n_simul > 0 && interval == true) 
        Dim  = buffer_t(); 
        ReDim  = buffer_t(bins);
        Initialise(buffer_t(), 0, bins);
        
        Dim  = all_rep_buffer(); 
        ReDim = all_rep_buffer(bins, n_simul);
        Initialise2d(all_rep_buffer, -999999999, bins, n_simul);
        
        Dim  = cramer_list(); 
        ReDim  = cramer_list(n_simul); % List for reps
        Dim  = cramer_p;
        cramer = monte_cramer(k_t, bins, t, 0, "K"); % Cramer from data
        
        Dim  = supremum_list(); 
        ReDim  = supremum_list(n_simul); % List for reps
        Dim  = supremum_p;
        supremum = monte_supremum(k_t, bins, t, 0, "K"); % Supremum from data
        
        Dim =  x_rnd(); 
        ReDim  = x_rnd(points);
        Dim = y_rnd(); 
        ReDim  = y_rnd(points);

        for rep = 1:n_simul
            pctDone = Int((rep / n_simul) * 100);

            if bygroups == false 
                Application.StatusBar = "Calculating " & testname & " CIs " & pctDone & "% completed."
            else
                Application.StatusBar = "Calculating " & testname & " CIs for group " & group & "; " & pctDone & "% completed."
            end

             generate_csr(x_rnd, y_rnd, xmin, ymin, xmax, ymax, points)
             Initialise(buffer_t(), 0, bins)
            
            if k_method == 1   
                calc_k(buffer_t, x_rnd, y_rnd, t, max_step, bins, points, area);
            elseif k_method == 2   
                calc_k_ew(buffer_t, x_rnd, y_rnd, t, max_step, bins, points, xmin, ymin, xmax, ymax, area);
            elseif k_method == 3   
                calc_k_toroid(buffer_t, x_rnd, y_rnd, t, max_step, bins, points, xmax, ymax, area);
            elseif k_method == 4   
                calc_k_buffer(buffer_t, x_rnd, y_rnd, t, max_step, bins, points, xmin, ymin, xmax, ymax, area);
            end
                
            for i = 1:bins
                all_rep_buffer(i, rep) = buffer_t(i);
            end
                        
            cramer_list(rep) = monte_cramer(buffer_t, bins, t, 0, "K");
            supremum_list(rep) = monte_supremum(buffer_t, bins, t, 0, "K");
            
        end

        sort_rows_array(all_rep_buffer, bins, n_simul);

        Dim lo_trans() As Single; ReDim lo_trans(bins)
        Dim hi_trans() As Single; ReDim hi_trans(bins)
        Dim low_list() As Single; ReDim low_list(bins)
        Dim high_list() As Single; ReDim high_list(bins)

        extract_limits(all_rep_buffer, low_list, high_list, bins, n_simul, limit);
        
        cramer_p = monte_p2t(dist_posn(cramer_list, cramer, n_simul, false), n_simul)
        Dim cramer_ave As Single; cramer_ave = Ave_fl_Array(cramer_list, n_simul, 1)
        Dim cramer_var As Single; cramer_var = Var_fl_Array(cramer_list, n_simul, 1)
        
        supremum_p = monte_p2t(dist_posn(supremum_list, supremum, n_simul, false), n_simul)
        Dim supremum_ave As Single; supremum_ave = Ave_fl_Array(supremum_list, n_simul, 1)
        Dim supremum_var As Single; supremum_var = Var_fl_Array(supremum_list, n_simul, 1)
        
        i = 1; 
        t_incr = 0;
        while i <= bins
            lo_trans(i) = Sqr(low_list(i) / (Pi)) - t_incr;
            hi_trans(i) = Sqr(high_list(i) / (Pi)) - t_incr;
            i = i + 1; 
            t_incr = t_incr + t;
        end
    end

    if bygroups == false 

         write_kt(k_t, kt_theo, low_list, high_list, kt_trans, lo_trans, hi_trans, (bins), interval, n_simul, t);

        if (n_simul > 0 && interval == true) 
             write_unisig(cramer, cramer_p, cramer_ave, cramer_var, supremum, supremum_p, supremum_ave, supremum_var);
        end
    elseif bygroups == true 
         write_kt_grp(k_t, kt_theo, low_list, high_list, kt_trans, lo_trans, hi_trans, ...
                           (bins), interval, n_simul, t, points, group, max_group, "K by Groups");

        if (n_simul > 0 && interval == true) 
             write_unisig_grpK(cramer, cramer_p, cramer_ave, cramer_var, ...
                                    supremum, supremum_p, supremum_ave, supremum_var, points, group, ...
                                    max_group, "K by Groups MC");
        end
    end
end