%%% Univariate Ripley's K estimator (second-order function)

% input data (hard coded)
x1; % x column vector
y1; % y column vector
points = size(x1); % number of values
t; % step size
max_step; % max distances we are performing the analysis to 
interval; % length of the interval 
n_simul;
limit;
k_method; % type of ripley k analysis performed 
bygroups;
group;
max_group;
strow;
sizeWarn = true;

% parameters
t_incr;
rep; % Counters
cramer;
supremum;
pctDone;
xmax = max(x1);
xmin = min(x1);
ymax = max(y1);
ymin = min(y1);
area = (ymax-ymin)*(xmax-xmin);

bins = ceil(max_step / t) + 1;
if (bins * t > (max_step + t))
    bins = bins - 1; % Make sure don't exceed max step ....
end

k_t = (bins + 1); 
kt_trans =  (bins + 1); 
kt_theo =  (bins + 1);

if (sizeWarn == true)
    checkscale(max_step, xmin, xmax, ymin, ymax);
end

% Empirical estimation of k(t) first
pctDone = 0;

if (k_method == 1)
    k_t = calc_k(x1, y1, t, max_step, bins, points, area);
elseif (k_method == 2)
    k_t = calc_k_ew(x1, y1, t, max_step, bins, points, xmin, ymin, xmax, ymax, area);
elseif (k_method == 3)
    k_t = calc_k_toroid(x1, y1, t, max_step, bins, points, xmax, ymax, area);
elseif (k_method == 4)
    k_t = calc_k_buffer(x1, y1, t, max_step, bins, points, xmin, ymin, xmax, ymax, area);
end

% save data