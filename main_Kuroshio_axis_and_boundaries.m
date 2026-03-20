% Kuroshio axis and boundaries.m
clc; clear;

%% ---------------- 1. Parameters & Initialization ----------------
data_path = '2024.mat';
mdata = matfile(data_path,'Writable',false);
lat   = double(mdata.lat);
lon   = double(mdata.lon);
depth = double(mdata.depth);
 
[LonGrid, LatGrid] = meshgrid(lon, lat);

kuroshio_data = struct();
load('mask.mat');

pre = 0.125;
domain = [115, 132; 15, 32];
start_longitudes = domain(1,1):pre:domain(1,2);
start_latitudes  = domain(2,1):pre:domain(2,2);
[xg, yg] = meshgrid(start_longitudes, start_latitudes);
starts = [xg(:), yg(:)];

valid_idx = (starts(:,2) < 24 & starts(:,1) >= 115 & starts(:,1) <= 124) | ...
            (starts(:,2) >= 24 & starts(:,1) >= 120 & starts(:,1) <= 132);

starts = starts(valid_idx,:);

[~, lon_idx] = min(abs(lon(:)' - starts(:,1)), [], 2);
[~, lat_idx] = min(abs(lat(:)' - starts(:,2)), [], 2);

lon_idx = max(1, min(lon_idx, length(lon)));
lat_idx = max(1, min(lat_idx, length(lat)));

linear_idx = sub2ind(size(mask_grid), lat_idx, lon_idx);
is_sea = mask_grid(linear_idx) >= 1;
starts = starts(is_sea, :);

min_loop_points = 8;
distance_threshold = 0.15;

tol = 0.5;
cs=30;
for yearIdx = 1993:2024
     data_path = sprintf('%d.mat', yearIdx);
     mdata = matfile(data_path,'Writable',false);
  days = sum(eomday(yearIdx,1:12));

if ~isfile(sprintf('kuroshio_%d.mat', yearIdx))
    kuroshio_data = repmat(struct('depth', []), 1, days);

    for d = 1:days
        kuroshio_data(d).depth = repmat(struct( ...
            'main_axis', [], ...
            'left_bdy', [], ...
            'right_bdy', [], ...
            'date', [], ...
            'depth_level', []), 1, cs);
    end

    save(sprintf('kuroshio_%d.mat', yearIdx), 'kuroshio_data', '-v7.3');
end

mf = matfile(sprintf('kuroshio_%d.mat', yearIdx), 'Writable', true);
kuroshio_data = mf.kuroshio_data;

last_done = 0;
for d = days:-1:1
    day_data = kuroshio_data(d).depth;
    if ~isempty(day_data) && ~isempty(day_data(1).main_axis)
        last_done = d;
        break;
    end
end

fprintf('Completed to day %d, continuing from day %d\n', last_done, last_done+1);

for dayIdx = last_done+1:days

   parfor dep = 1:cs
         try  
        u0 = squeeze(mdata.u_all(:, :,dep, dayIdx));
    v0 = squeeze(mdata.v_all(:, :, dep, dayIdx));
    sub_u = single(u0);  sub_u(isnan(sub_u)|mask_grid < 1) = 0;
    sub_v = single(v0);  sub_v(isnan(sub_v)|mask_grid < 1) = 0;
        Fu = griddedInterpolant(LatGrid, LonGrid, sub_u, 'linear');
        Fv = griddedInterpolant(LatGrid, LonGrid, sub_v, 'linear');

        xy_pool = cell(size(starts,1),1);

        for i = 1:size(starts,1)
            xp = starts(i,1); 
            yp = starts(i,2);

            fwd = stream2(LonGrid, LatGrid, sub_u, sub_v, xp, yp, 1);
            bwd = stream2(LonGrid, LatGrid, -sub_u, -sub_v, xp, yp, 1);

            xy1 = []; xy2 = [];
            if ~isempty(bwd), xy1 = flipud(bwd{1}); end
            if ~isempty(fwd), xy2 = fwd{1}; end

            if isempty(xy1) && isempty(xy2)
                xy_full = [];
            elseif isempty(xy1)
                xy_full = xy2;
            elseif isempty(xy2)
                xy_full = xy1;
            else
                if norm(xy1(end,:) - xy2(1,:)) < 1e-5
                    xy_full = [xy1; xy2(2:end,:)];
                else
                    xy_full = [xy1; xy2];
                end
            end

        xy_pool{i} = xy_full;
        end

        xy_pool = xy_pool(~cellfun('isempty',xy_pool));

seg_def = [ ...
    struct('name','A', 'cond', @(lat,lon) lat<=19, 'lat_lim',[16,19], ...
           'south_sec',[122,16,124,16], 'north_sec',[120.85,19,123,19], ...
           'lon_range',[120.85,124], 'mid_sec',[122.22,17.5,123.5,17.5]), ...
    struct('name','B', 'cond', @(lat,lon) lat>=19 & lat<=22, 'lat_lim',[19,22], ...
           'south_sec',[120.85,19,123,19], 'north_sec',[120.85,22,123,22], ...
           'lon_range',[118,123], 'mid_sec',[118,20.5,123,20.5]), ...
    struct('name','C', 'cond', @(lat,lon) lat>=22 & lat<=25, 'lat_lim',[22,25], ...
           'south_sec',[120.85,22,123,22], 'north_sec',[122,25,124,25], ...
           'lon_range',[120.85,124], 'mid_sec',[121.49,23.5,123.5,23.5]), ...
    struct('name','D', 'cond', @(lat,lon) lat>=25& lat<=29, 'lat_lim',[25,29], ...
           'south_sec',[122,25,124,25], 'north_sec',[126,29,127.85,27.77], ...
           'lon_range',[122,127.85], 'mid_sec',[125,27.375,126.85,26.145]), ...
    struct('name','E', 'cond', @(lat,lon) lat>=27.77& lat<=31.25, 'lat_lim',[27.77,31.25], ...
           'south_sec',[126,29,127.85,27.77], 'north_sec',[130.38,28.25,130.38,31.25], ...
           'lon_range',[126,130.38], 'mid_sec',[127,30.625,128.85,29.395]) ...
]; 

xy_seg = cell(size(seg_def,2),1); 
xy_segbei = cell(2,1);
prev_end = [];

for seg_idx = 1:size(seg_def,2)
    seg = seg_def(seg_idx);
    cond_func = seg.cond;
    lat_lim   = seg.lat_lim;
    name      = seg.name;
    lon_range = seg.lon_range;
    south_section = seg.south_sec;
    north_section = seg.north_sec;
    mid_section = seg.mid_sec;

    candidates = struct('xy',{},'start_dist',{},'avg_speed',{});

    for k = 1:numel(xy_pool)
        xy = xy_pool{k};
        lat_xy = xy(:,2);  lon_xy = xy(:,1);

        idx_s = find(is_on_section_fast(lon_xy,lat_xy, south_section(1),south_section(2),south_section(3),south_section(4),tol));
        idx_n = find(is_on_section_fast(lon_xy,lat_xy, north_section(1),north_section(2),north_section(3),north_section(4),tol));

        if numel(idx_s)<1 || numel(idx_n)<1, continue; end

        for ii = numel(idx_s):-1:1
            s = idx_s(ii);
            idx = idx_n(idx_n > s);
            if isempty(idx), continue; end
            e = idx(1);
            seg_xy = xy(s:e,:);

            if ~validate_streamline(seg_xy, seg, min_loop_points, distance_threshold)
                continue;
            end

            if isempty(prev_end)
                start_dist = 0;
            else
                start_dist = hypot(seg_xy(1,1)-prev_end(1), seg_xy(1,2)-prev_end(2));
            end
            uu = Fu(seg_xy(:,2),seg_xy(:,1));
            vv = Fv(seg_xy(:,2),seg_xy(:,1));
            sp = mean(hypot(uu,vv),'omitnan');

            candidates(end+1).xy         = seg_xy;
            candidates(end).start_dist   = start_dist;
            candidates(end).avg_speed    = sp;
            break;
        end
    end

    if ~isempty(candidates)
        [best_idx, best_seg] = pick_best_candidate(candidates);
        xy_seg{seg_idx} = best_seg.xy;
        prev_end        = best_seg.xy(end,:);
    else
        prev_end_local = prev_end;

        candidates_south_mid = struct('xy',{},'start_dist',{},'avg_speed',{});
        for k = 1:numel(xy_pool)
            xy = xy_pool{k};
            lat_xy = xy(:,2); lon_xy = xy(:,1);

            idx_s = find(is_on_section_fast(lon_xy,lat_xy, south_section(1),south_section(2),south_section(3),south_section(4),tol));
            idx_m = find(is_on_section_fast(lon_xy,lat_xy, mid_section(1),mid_section(2),mid_section(3),mid_section(4),tol));

            if isempty(idx_s) || isempty(idx_m), continue; end

            for ii = 1:numel(idx_s)
                s_temp = idx_s(ii);
                m_idx = idx_m(idx_m > s_temp);
                if isempty(m_idx), continue; end

                seg_xy = xy(s_temp:m_idx(1),:);
                if ~validate_streamline(seg_xy, seg, min_loop_points, distance_threshold), continue; end

                if isempty(prev_end_local)
                    start_dist = 0;
                else
                    start_dist = hypot(seg_xy(1,1)-prev_end_local(1), seg_xy(1,2)-prev_end_local(2));
                end
                uu = Fu(seg_xy(:,2),seg_xy(:,1));
                vv = Fv(seg_xy(:,2),seg_xy(:,1));
                sp = mean(hypot(uu,vv),'omitnan');

                candidates_south_mid(end+1).xy = seg_xy;
                candidates_south_mid(end).start_dist = start_dist;
                candidates_south_mid(end).avg_speed = sp;
            end
        end

        if ~isempty(candidates_south_mid)
            [best_idx, best_seg] = pick_best_candidate(candidates_south_mid);
            xy_segbei{1} = best_seg.xy;
            prev_end_local = best_seg.xy(end,:);
        end

        candidates_mid_north = struct('xy',{},'start_dist',{},'avg_speed',{});
        for k = 1:numel(xy_pool)
            xy = xy_pool{k};
            lat_xy = xy(:,2); lon_xy = xy(:,1);

            idx_m = find(is_on_section_fast(lon_xy,lat_xy, mid_section(1),mid_section(2),mid_section(3),mid_section(4),tol));
            idx_n = find(is_on_section_fast(lon_xy,lat_xy, north_section(1),north_section(2),north_section(3),north_section(4),tol));

            if isempty(idx_m) || isempty(idx_n), continue; end

            for ii = 1:numel(idx_m)
                m_temp = idx_m(ii);
                n_idx = idx_n(idx_n > m_temp);
                if isempty(n_idx), continue; end

                seg_xy = xy(m_temp:n_idx(1),:);
                if ~validate_streamline(seg_xy, seg, min_loop_points, distance_threshold), continue; end
                if isempty(prev_end_local)
                    start_dist = 0;
                else
                    start_dist = hypot(seg_xy(1,1)-prev_end_local(1), seg_xy(1,2)-prev_end_local(2));
                end

                uu = Fu(seg_xy(:,2),seg_xy(:,1));
                vv = Fv(seg_xy(:,2),seg_xy(:,1));
                sp = mean(hypot(uu,vv),'omitnan');

                candidates_mid_north(end+1).xy = seg_xy;
                candidates_mid_north(end).start_dist = start_dist;
                candidates_mid_north(end).avg_speed = sp;
            end
        end

        if ~isempty(candidates_mid_north)
            [best_idx, best_seg] = pick_best_candidate(candidates_mid_north);
            xy_segbei{2} = best_seg.xy;
            prev_end_local = best_seg.xy(end,:);
        end

        xy_seg{seg_idx} = vertcat(xy_segbei{:});
    end
end

reference_axis = vertcat(xy_seg{:});

[reference_axis,removed_idx] = smooth_path_bd(reference_axis,'gauss',recommend_window_radius(reference_axis,55),Fmask);
reference_axis(removed_idx, :) = [];

main_axis= correct_axis_transect(Fu, Fv, reference_axis, Fmask);

start_sec = [122, 16, 124, 16];
end_sec = [130.38, 28.25, 130.38, 31.25];

main_axis = cut_reference_axis(main_axis, start_sec, end_sec);

[left_bdy, right_bdy,removed_idx1,removed_idx2] = find_kuroshio_boundaries(main_axis, Fu, Fv,Fmask, 0.3);
left_bdy(removed_idx1,:) = nan;
right_bdy(removed_idx2,:) = nan;

tmp = struct( ...
    'main_axis',      main_axis, ...
    'left_bdy',       left_bdy, ...
    'right_bdy',      right_bdy, ...
    'date',           dayIdx, ...
    'depth_level',    dep);

     catch
        tmp = struct( ...
            'main_axis',      [], ...
            'left_bdy',       [], ...
            'right_bdy',      [], ...
            'date',           dayIdx, ...
            'depth_level',    dep);
     end
temp_result{dep} = tmp; 
  end

    kuroshio_data(dayIdx).depth = [temp_result{:}];
mf.kuroshio_data(1,dayIdx) = kuroshio_data(dayIdx);
    clear temp_result xy_pool tmp;
    
        fid = fopen('progress.log','a');
        fprintf(fid,'Completed day %d\n', dayIdx);
        fclose(fid);
end
end

function [best_idx, best_seg] = pick_best_candidate(candidates)
% Select best candidate streamline based on distance and speed
alpha = 0.5;

dist_vec = [candidates.start_dist];
spd_vec  = [candidates.avg_speed];

if numel(dist_vec) > 1
    dist_norm = (dist_vec - min(dist_vec)) ./ (max(dist_vec) - min(dist_vec) + eps);
else
    dist_norm = 0;
end

if numel(spd_vec) > 1
    spd_norm = (spd_vec - min(spd_vec)) ./ (max(spd_vec) - min(spd_vec) + eps);
else
    spd_norm = 1;
end

score = alpha*dist_norm + (1-alpha)*(1 - spd_norm);

[~, best_idx] = min(score);
best_seg = candidates(best_idx);
end

function is_valid = validate_streamline(seg_xy, seg, min_loop_points, distance_threshold)
% Validate streamline against constraints
    if any(seg_xy(:,1) < seg.lon_range(1) | seg_xy(:,1) > seg.lon_range(2)) || ...
       any(seg_xy(:,2) < seg.lat_lim(1) | seg_xy(:,2) > seg.lat_lim(2))
        is_valid = false;
        return;
    end
    
    has_loop = detect_loops(seg_xy, min_loop_points, distance_threshold);
    if has_loop
        is_valid = false;
        return;
    end
    
    is_valid = true;
end

function mask = is_on_section_fast(lon, lat, lon1, lat1, lon2, lat2, tol)
% Fast check if points are on section line
    v = [lon2-lon1, lat2-lat1];
    w = [lon(:)-lon1, lat(:)-lat1];
    b = (w * v.') ./ sum(v.^2);
    b = max(0,min(1,b));
    pb = [lon1,lat1] + b .* v;
    dist = hypot(lon(:)-pb(:,1), lat(:)-pb(:,2));
    mask = (dist <= tol);
end

function ref_axis_cut = cut_reference_axis(reference_axis, start_sec, end_sec)
% Cut reference axis between start and end sections
lon = reference_axis(:,1);
lat = reference_axis(:,2);

north_idx = find(lat >= start_sec(2), 1);
if isempty(north_idx)
    north_idx = length(lat);
end
idx_start = north_idx;

post_idx = (idx_start+1):length(lon);
west_idx = find(lon(post_idx) <= end_sec(1));

if ~isempty(west_idx)
    idx_end = post_idx(west_idx(end));
else
    [~, min_idx] = min(abs(lon(post_idx) - end_sec(1)));
    idx_end = post_idx(min_idx);
end

idx_end = max(idx_end, idx_start);
ref_axis_cut = reference_axis(idx_start:idx_end, :);
end

function [corrected_path, removed_idx] = correct_axis_transect(Fu, Fv, ref_path, Fmask)
% Correct axis using cross-section velocity analysis
    half_len_km = 70; n_sect = 20;
    N = size(ref_path,1);
    corrected_path = nan(N,2);

    transition_radius = recommend_window_radius(ref_path,110);

    for i = 1:N
        u0 = Fu(ref_path(i,2), ref_path(i,1));
        v0 = Fv(ref_path(i,2), ref_path(i,1));
        if isnan(u0)||isnan(v0)||(u0==0&&v0==0)
            corrected_path(i,:) = ref_path(i,:); 
            continue;
        end

 alpha = 1;
u0 = Fu(ref_path(i,2), ref_path(i,1));
v0 = Fv(ref_path(i,2), ref_path(i,1));

if isnan(u0) || isnan(v0) || (u0==0 && v0==0)
    theta_ref = NaN;
else
    theta_ref = atan2(v0, u0);
end

        deg_lat = half_len_km / 111;
        deg_lon = half_len_km / (111*cosd(ref_path(i,2)));

        perp_theta0 = theta_ref + pi/2;
        delta_lon = deg_lon * cos(perp_theta0);
        delta_lat = deg_lat * sin(perp_theta0);
        lon_sect = linspace(ref_path(i,1)-delta_lon, ref_path(i,1)+delta_lon, n_sect);
        lat_sect = linspace(ref_path(i,2)-delta_lat, ref_path(i,2)+delta_lat, n_sect);

        u_sect = Fu(lat_sect, lon_sect);
        v_sect = Fv(lat_sect, lon_sect);
        speed_sect = hypot(u_sect, v_sect);

        perp_theta = theta_ref + pi/2;
        delta_lon2 = deg_lon * cos(perp_theta);
        delta_lat2 = deg_lat * sin(perp_theta);
        lon_sect2 = linspace(ref_path(i,1)-delta_lon2, ref_path(i,1)+delta_lon2, n_sect);
        lat_sect2 = linspace(ref_path(i,2)-delta_lat2, ref_path(i,2)+delta_lat2, n_sect);

        u_sect2 = Fu(lat_sect2, lon_sect2);
        v_sect2 = Fv(lat_sect2, lon_sect2);
        V_proj = u_sect2 .* cos(theta_ref) + v_sect2 .* sin(theta_ref);

mask_sect = Fmask(lat_sect2, lon_sect2);
V_proj_masked = V_proj;
V_proj_masked(mask_sect < 1) = 0;

[~, idx_max] = max(V_proj_masked, [], 'omitnan'); 
if isempty(idx_max) || isnan(idx_max)
    corrected_path(i,:) = ref_path(i,:);
else
    target = [lon_sect2(idx_max), lat_sect2(idx_max)];
    corrected_path(i,:) = alpha*target + (1-alpha)*ref_path(i,:);
end
    end

    [corrected_path,removed_idx] = smooth_path_bd(corrected_path, 'gauss', transition_radius, Fmask);
    corrected_path(removed_idx,:)=[];
end

function [left_boundary, right_boundary, removed_idx1, removed_idx2] = find_kuroshio_boundaries(ref_path, Fu, Fv, Fmask, pct_speed)
% Find Kuroshio boundaries with geometric correction
N = size(ref_path,1);
left_boundary  = nan(N,2);
right_boundary = nan(N,2);

max_dist_km = 70;
step_km = max_dist_km/10;

for i = 1:N
    lon0 = ref_path(i,1); lat0 = ref_path(i,2);
    u0 = Fu(lat0,lon0); v0 = Fv(lat0,lon0);
    if isnan(u0)||isnan(v0)||(u0==0&&v0==0), continue; end

T_ref = [u0, v0];
T_ref = T_ref / norm(T_ref);

    N_left  = [-T_ref(2), T_ref(1)]; N_left  = N_left / norm(N_left);
    N_right = [ T_ref(2),-T_ref(1)]; N_right = N_right / norm(N_right);

    thr = abs(u0*T_ref(1) + v0*T_ref(2)) * pct_speed;

    left_boundary(i,:)  = search_boundary_with_check_km(lon0, lat0, N_left, thr, Fu,Fv,Fmask, max_dist_km, step_km, T_ref, true);
    right_boundary(i,:) = search_boundary_with_check_km(lon0, lat0, N_right,thr, Fu,Fv,Fmask, max_dist_km, step_km, T_ref, true);
end

[left_boundary, removed_idx1]  = smooth_path_bd(left_boundary,'gauss', recommend_window_radius(left_boundary,110), Fmask);
[right_boundary, removed_idx2] = smooth_path_bd(right_boundary,'gauss', recommend_window_radius(right_boundary,110), Fmask);

left_boundary  = drop_decreasing_lat(left_boundary);
right_boundary = drop_decreasing_lat(right_boundary);
[left_boundary, right_boundary] = mask_side_flip(ref_path, left_boundary, right_boundary);
end

function [left, right] = mask_side_flip(ref_path, left, right)
% Ensure boundaries stay on correct sides of reference path
N = size(ref_path,1);

for i = 1:N
    if any(isnan(ref_path(i,:))) || any(isnan(left(i,:))) || any(isnan(right(i,:)))
        continue;
    end
    
    if i < N
        T = ref_path(i+1,:) - ref_path(i,:);
    else
        T = ref_path(i,:) - ref_path(i-1,:);
    end
    T = T / norm(T);
    N_left = [-T(2), T(1)];
    
    side_left = dot(left(i,:) - ref_path(i,:), N_left);
    side_right = dot(right(i,:) - ref_path(i,:), N_left);
    
    if side_left <= 0
        left(i,:) = NaN;
    end
    if side_right >= 0
        right(i,:) = NaN;
    end
end
end

function pts = drop_decreasing_lat(pts)
% Remove points with decreasing latitude and interpolate gaps
N = size(pts,1);
if N<2, return; end

process_mask = pts(:,1) > 120.85;

if ~any(process_mask)
    return;
end

pts_west = pts(~process_mask,:);

validIdx = find(~any(isnan(pts),2) & process_mask);
if isempty(validIdx), return; end

[~, maxPos] = max(pts(validIdx,2));
maxIdx = validIdx(maxPos);

lastLat = pts(1,2);
for i = 2:maxIdx
    if ~process_mask(i) || isnan(pts(i,1))
        continue;
    end
    if pts(i,2) <= lastLat
        pts(i,:) = NaN;
    else
        lastLat = pts(i,2);
    end
end

lastLon = pts(maxIdx,1);
for i = maxIdx+1:N
    if ~process_mask(i) || isnan(pts(i,1))
        continue;
    end
    if pts(i,1) <= lastLon
        pts(i,:) = NaN;
    else
        lastLon = pts(i,1);
    end
end

pts = fill_nan_segments_selective(pts, process_mask);
pts(~process_mask,:) = pts_west;
end

function pts = fill_nan_segments_selective(pts, process_mask)
% Fill NaN segments with linear interpolation in selective regions
N = size(pts,1);
i = 1;

while i <= N
    if process_mask(i) && any(isnan(pts(i,:)))
        gap_start = i;
        
        while i <= N && process_mask(i) && any(isnan(pts(i,:)))
            i = i + 1;
        end
        gap_end = i - 1;
        
        if gap_start > 1 && gap_end < N && ...
           ~any(any(isnan(pts([gap_start-1, gap_end+1],:))))
            
            pre_val = pts(gap_start-1,:);
            post_val = pts(gap_end+1,:);
            
            for j = gap_start:gap_end
                t = (j - gap_start + 1) / (gap_end - gap_start + 2);
                pts(j,:) = (1-t)*pre_val + t*post_val;
            end
        end
    else
        i = i + 1;
    end
end
end

function pos = search_boundary_with_check_km(lon0, lat0, dir_vec, threshold, Fu, Fv, Fmask, max_dist_km, step_km, T, check_direction)
% Vectorized boundary search with direction checking
    pos = [nan, nan];
    last_valid = [nan, nan];
    farthest_directional = [nan, nan];
    farthest_distance = 0;

    steps = floor(max_dist_km / step_km);
    s_idx = 1:steps;
    
    delta_lat_deg = (step_km * s_idx) / 111;
    lat_tests = lat0 + dir_vec(2) * delta_lat_deg;
    km_per_deg_lon = max(111 * cosd(lat_tests), 1e-6);
    delta_lon_deg = (step_km * s_idx .* dir_vec(1)) ./ km_per_deg_lon;
    lon_tests = lon0 + delta_lon_deg;

    u_all = Fu(lat_tests, lon_tests);
    v_all = Fv(lat_tests, lon_tests);
    mask_all = Fmask(lat_tests, lon_tests);

    for s = 1:steps
        u_t = u_all(s);
        v_t = v_all(s);
        mask_t = mask_all(s);
        lon_t = lon_tests(s);
        lat_t = lat_tests(s);

        if isnan(u_t) || isnan(v_t) || isnan(mask_t) || mask_t < 1
            break;
        end

        last_valid = [lon_t, lat_t];

        V_proj = u_t*T(1) + v_t*T(2);

        if check_direction
            speed_mag = hypot(u_t, v_t);
            if speed_mag > 0
                dir_cos = (u_t*T(1) + v_t*T(2)) / speed_mag;
                if dir_cos >0
                    dist = sqrt((lon_t - lon0)^2 + (lat_t - lat0)^2);
                    if dist > farthest_distance
                        farthest_directional = [lon_t, lat_t];
                        farthest_distance = dist;
                    end
                    if abs(V_proj) <= threshold
                        pos = [lon_t, lat_t];
                        return;
                    end
                end
            end
        else
            if abs(V_proj) <= threshold
                pos = [lon_t, lat_t];
                return;
            end
        end
    end

    if isnan(pos(1))
        if ~any(isnan(farthest_directional))
            pos = farthest_directional;
        elseif ~any(isnan(last_valid))
            pos = last_valid;
        else
            pos = [lon0, lat0];
        end
    end
end

function [smoothed_path, removed_idx] = smooth_path_bd(path, method, window_radius, Fmask)
% Smooth path with land/sea checking
N = size(path,1);

is_sea = Fmask(double(path(:,2)), double(path(:,1))) == 1;

sea_segments = {};
start_idx = [];
for i = 1:N
    if is_sea(i)
        if isempty(start_idx), start_idx = i; end
    else
        if ~isempty(start_idx)
            sea_segments{end+1} = start_idx:i-1;
            start_idx = [];
        end
    end
end
if ~isempty(start_idx), sea_segments{end+1} = start_idx:N; end

smoothed_path = nan(size(path));
for s = 1:length(sea_segments)
    idx = sea_segments{s};
    seg = path(idx,:);
    seg_N = size(seg,1);

    if seg_N < 2
        smoothed_path(idx,:) = seg;
        continue;
    end

    seg_smoothed = nan(seg_N,2);
    for dim = 1:2
        vals = seg(:,dim);
        for i = 1:seg_N
            cur_radius = min([window_radius, i-1, seg_N-i]);
            idx_window = (i-cur_radius):(i+cur_radius);

            w_local = gausswin(length(idx_window));
            w_local = w_local / sum(w_local);

            vals_window = vals(idx_window);
            mask_window = is_sea(idx_window);
            w_window = w_local .* mask_window;
            if sum(w_window)==0
                seg_smoothed(i,dim) = NaN;
            else
                w_window = w_window / sum(w_window);
                seg_smoothed(i,dim) = sum(vals_window .* w_window);
            end
        end
    end
    smoothed_path(idx,:) = seg_smoothed;
end

mask_check = Fmask(double(smoothed_path(:,2)), double(smoothed_path(:,1)));
nan_idx  = find(any(isnan(smoothed_path),2));
land_idx = find(mask_check < 1);
removed_idx = unique([nan_idx; land_idx]);
smoothed_path(removed_idx,:) = NaN;
end

function best_radius = recommend_window_radius(path,target_km)
% Calculate optimal window radius for smoothing
    R=6371;
    lat=path(:,2)*pi/180; lon=path(:,1)*pi/180;
    a=sin(diff(lat)/2).^2+cos(lat(1:end-1)).*cos(lat(2:end)).*sin(diff(lon)/2).^2;
    segment_dist=2*R*atan2(sqrt(a),sqrt(1-a));
    mean_spacing=mean(segment_dist);
    best_radius=max(2,ceil(1.5*target_km/mean_spacing));
end

function has_loop = detect_loops(xy, min_loop_points, distance_threshold)
% Detect loop structures in streamline
    n = size(xy, 1);
    has_loop = false;
    
    if n < min_loop_points * 2
        return;
    end
    
    for i = 1:n-min_loop_points
        for j = i+min_loop_points:n
            dist = sqrt((xy(i,1)-xy(j,1))^2 + (xy(i,2)-xy(j,2))^2);
            
            if dist < distance_threshold
                has_loop = true;
                return;
            end
        end
    end
end
