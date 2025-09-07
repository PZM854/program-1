code in dox:

case_to_graph.m:
function G = case_to_graph(case_WT, type)

    define_constants;
    % case_WT:outcome of powerflow analysis

    % type is 'graph': make non-direction graph
    % type is 'diagraph': make non-direction graph

    % case to graph, only produce raw graph
    
    p = inputParser;
    addRequired(p, 'case_WT', @(x) 1);
    addRequired(p, 'type', @(x) ischar(x) || isstring(x));
    parse(p, case_WT, type);

    type = string(lower(type));

    [~, ~, table_branch_ope] = powerflow_to_table(case_WT);

    if strcmp(type,'graph') || strcmpi(type,"graph")
        G = graph(table_branch_ope);
        
    elseif strcmp(type,'digraph') || strcmpi(type,"digraph")
        G = digraph(table_branch_ope);

    else
        error('type must be "graph" or "digraph", got: %s', type);
    end
        G.Nodes.Volt = case_WT.bus(:, BASE_KV);
        G.Nodes.uni_index = [1:size(case_WT.bus,1)]';
        G.Nodes.matpower_index = case_WT.bus(:, BUS_I); %uni_index
        G.Nodes.angle = case_WT.bus(:, VA);

end

color_colorbrewer.m:
function color = colorbrewer(num_color)

% 12 colors
color = [166,206,227
31,120,180
178,223,138
51,160,44
251,154,153
227,26,28
253,191,111
255,127,0
202,178,214
106,61,154
255,255,153
177,89,40];
color = color/255;

color_idx = floor(linspace(1,size(color, 1), num_color));
color = color(color_idx,:);

end


edge_combination.m:
function edge_new = edge_combination(T, idx_weight, node_key)

% take in a table that about to form a digraph, or a matrix that to form a
% table

% T: table or matrix
% idx_weight: the index of column represents the weight
% node_key: the index of column represents fromnodes and tonodes

% combine those edge have the same from/to nodes by adding up the weights

    p = inputParser;
    addRequired(p, 'T', @(x) istable(x) || isnumeric(x));
    addRequired(p, 'idx_weight', @(x) isnumeric(x) || ischar(x) || isstring(x));
    addRequired(p, 'node_key', @(x) iscell(x)&&numel(x)==2 || isnumeric(x)&&numel(x)==2 ...
        || ischar(x) || isstring(x));
    parse(p, T, idx_weight, node_key);

    if istable(p.Results.T)
        T_new = p.Results.T{:,:};
    else
        T_new = p.Results.T;
    end
  %  T = p.Results.T;
    idx_weight = p.Results.idx_weight;
    node_key = p.Results.node_key;

    if ischar(idx_weight) || isstring(idx_weight)
        idx_weight = get_field_col_range(idx_weight, T);
    end

    if iscell(node_key)
        from = node_key{1};
        to = node_key{2};
        from = get_field_col_range(from, T);
        to = get_field_col_range(to, T);
    elseif ischar(node_key) || isstring(node_key)
        from_to = get_field_col_range(node_key, T);
        from = from_to(1);
        to = from_to(2);
    else
        from = node_key(1);
        to = node_key(2);
    end

    T_new = sortrows(T_new, [from to]);

    num_edges = size(T_new,1);
    edge_new = [];
    p_l = 1;
    p_r = 1;

    while(1)
        if p_r == num_edges
            if T_new(p_r,from)==T_new(p_l,from) && T_new(p_r,to)==T_new(p_l,to)
                edge_new = [edge_new; [T_new(p_l,from) T_new(p_l,to) sum(T_new(p_l:p_r,idx_weight))]];
            else
                edge_new = [edge_new; [T_new(p_r,from) T_new(p_r,to) T_new(p_r,idx_weight)]];
            end
            
            break;
        elseif T_new(p_r,from)==T_new(p_l,from) && T_new(p_r,to)==T_new(p_l,to)
            p_r = p_r + 1;
        else
            edge_new = [edge_new; [T_new(p_l,from) T_new(p_l,to) sum(T_new(p_l:p_r-1,idx_weight))]];
            p_l = p_r;
        end
    end
    
    % if istable(T)  % T is a table 
    %     while(1)
    %         if p_r == num_edges
    %             edge_new = [edge_new; [T.(from)(p_l) T.(to)(p_l) sum(T.(idx_weight)(p_l:p_r))]];
    %             break;
    %         elseif T.(from)(p_r)==T.(from)(p_l) && T.(to)(p_r)==T.(to)(p_l)
    %             p_r = p_r + 1;
    %         else
    %             edge_new = [edge_new; [T.(from)(p_l) T.(to)(p_l) sum(T.(idx_weight)(p_l:p_r-1))]];
    %             p_l = p_r;
    %         end
    %     end
    % else % T is a matrix
    %     while(1)
    %         if p_r == num_edges
    %             edge_new = [edge_new; [T(p_l,from) T(p_l,to) sum(T(p_l:p_r,idx_weight))]];
    %             break;
    %         elseif T(p_r,from)==T(p_l,from) && T(p_r,to)==T(p_l,to)
    %             p_r = p_r + 1;
    %         else
    %             edge_new = [edge_new; [T(p_l,from) T(p_l,to) sum(T(p_l:p_r-1,idx_weight))]];
    %             p_l = p_r;
    %         end
    %     end
    % end

    %edge_new(:,4) = [1:size(edge_new,1)]';

end


edge_offset.m:
function T_new = edge_offset(T, idx_weight, node_key)

% offset the weights of branches of reverse directions

% assume that branches have been combined. Each from/to pair have no
% more than one branch.
% T: table or matrix
% idx_weight: the index of column represents the weight
% node_key: the index of column represents fromnodes and tonodes

    p = inputParser;
    addRequired(p, 'T', @(x) istable(x) || isnumeric(x));
    addRequired(p, 'idx_weight', @(x) isnumeric(x) || ischar(x) || isstring(x));
    addRequired(p, 'node_key', @(x) iscell(x)&&numel(x)==2 || isnumeric(x)&&numel(x)==2);
    parse(p, T, idx_weight, node_key);
    
    if istable(p.Results.T)
        T_new = p.Results.T{:,:};
    else
        T_new = p.Results.T;
    end
    %T_new = p.Results.T;
    idx_weight = p.Results.idx_weight;
    node_key = p.Results.node_key;

    if ischar(idx_weight) || isstring(idx_weight)
        idx_weight = get_field_col_range(idx_weight, T);
    end

    if iscell(node_key)
        from = node_key{1};
        to = node_key{2};
        from = get_field_col_range(from, T);
        to = get_field_col_range(to, T);
    elseif ischar(node_key) || isstring(node_key)
        from_to = get_field_col_range(node_key, T);
        from = from_to(1);
        to = from_to(2);
    else
        from = node_key(1);
        to = node_key(2);
    end

    point = 1;

    while(point <= size(T_new,1))
            rev_edge = (T_new(:, from)==T_new(point, to)) &...
                     (T_new(:, to)==T_new(point, from));
            idx_rev = find(rev_edge);
            if ~isempty(idx_rev)
                T_new(point, idx_weight) = T_new(point, idx_weight) - T_new(idx_rev, idx_weight);
                T_new(idx_rev,:) = [];
            end
            point = point + 1;
    end
    for i = 1:size(T_new,1)
        if T_new(i, idx_weight) < 0
            temp = T_new(i,from);
            T_new(i,from) = T_new(i,to);
            T_new(i,to) = temp;
            T_new(i,idx_weight) = -T_new(i,idx_weight);
        end
    end

    % if istable(T_new)  % T is a table 
    %     while(point <= size(T_new,1))
    %         rev_edge = (T_new(:, from)==T_new(point, to)) &...
    %                  (T_new(:, to)==T_new(point, from));
    %         idx_rev = find(rev_edge);
    %         if ~isempty(idx_rev)
    %             T_new(point, idx_weight) = T_new(point, idx_weight) - T_new(idx_rev, idx_weight);
    %             T_new(idx_rev,:) = [];
    %         end
    %         point = point + 1;
    %     end
    %     for i = 1:size(T_new,1)
    %         if T_new(i, idx_weight) < 0
    %             temp = T_new(i,from);
    %             T_new(i,from) = T_new(i,to);
    %             T_new(i,to) = temp;
    %             T_new(i,idx_weight) = -T_new(i,idx_weight);
    %         end
    %     end
    % 
    % 
    % else % T is a matrix
    %     while(point <= size(T_new,1))
    %         rev_edge = (T_new(:, from)==T_new(point, to)) &...
    %                  (T_new(:, to)==T_new(point, from));
    %         idx_rev = find(rev_edge);
    %         if ~isempty(idx_rev)
    %             T_new(point, idx_weight) = T_new(point, idx_weight) - T_new(idx_rev, idx_weight);
    %             T_new(idx_rev,:) = [];
    %         end
    %         point = point + 1;
    %     end
    %     for i = 1:size(T_new,1)
    %         if T_new(i, idx_weight) < 0
    %             temp = T_new(i,from);
    %             T_new(i,from) = T_new(i,to);
    %             T_new(i,to) = temp;
    %             T_new(i,idx_weight) = -T_new(i,idx_weight);
    %         end
    %     end
    % end
end


find_single_edge.m:
function idx_single = find_single_edge(T, from_to)

    define_constants;       


    from = from_to(1);
    to = from_to(2);

    idx_single = [];

    for k = 1:size(T,1)
        idx_rev = find(T(:,to)==T(k,from)&T(:, from)==T(k,to), 1);
        if isempty(idx_rev)
            idx_single = [idx_single k];
        end

    end

end


get_X_data.m:
function x_coord = get_X_data(G, y_coord)
    % input：
    %   G       - digraph object（zone scale）
    %   y_coord - y coordinate of each node (fixed)
    % output：
    %   x_coord - x coordinate of each node

    num_nodes = numnodes(G);    %23
    x_coord = zeros(num_nodes, 1);  %[0, 0...0] 

    % 
    unique_y = sort(unique(y_coord), 'descend');    %[-1,..,-10]
    num_layers = length(unique_y);  %10

    % 
    layer_of_node = zeros(num_nodes, 1);
    for i = 1:num_layers
        layer_of_node(y_coord == unique_y(i)) = i;
    end

    % 
    for i = 2:num_layers
        curr_nodes = find(layer_of_node == i);
        barycenters = zeros(length(curr_nodes), 1);

        for j = 1:length(curr_nodes)
            node = curr_nodes(j);
            preds = predecessors(G, node);  % 
            if ~isempty(preds)
                barycenters(j) = mean(x_coord(preds));
            else
                barycenters(j) = j; % 
            end
        end

        % 
        [~, sort_idx] = sort(barycenters);
        spacing = 10; 
        num_curr = length(curr_nodes);
        center = (num_curr + 1) / 2;
        for k = 1:num_curr
            x_coord(curr_nodes(sort_idx(k))) = (k - center) * spacing;
        end
    end

    % 
    top_nodes = find(layer_of_node == 1);
    spacing = 10;
    num_top = length(top_nodes);
    center = (num_top + 1) / 2;
    for k = 1:num_top
        x_coord(top_nodes(k)) = (k - center) * spacing;
    end
end


get_field_col_range.m:
function col_range = get_field_col_range(target_field, T)

    varnames = T.Properties.VariableNames;
    n_fields = numel(varnames);

    % 统计每个字段的列宽
    field_widths = zeros(1, n_fields);
    for i = 1:n_fields
        field_widths(i) = size(T.(varnames{i}), 2);
    end

    % 计算各字段在数组中的起始、终止列号
    start_col = cumsum([1, field_widths(1:end-1)]);
    end_col   = cumsum(field_widths);

    % 找目标字段下标
    idx = find(strcmp(varnames, target_field));
    if isempty(idx)
        error('字段名 "%s" 不存在！', target_field);
    end

    % 输出目标字段在数组中的范围
    start_col = start_col(idx);
    end_col = end_col(idx);

    if start_col == end_col
        col_range = end_col;
    else
        col_range = [start_col end_col];
    end

end
% width = [1 2 4 3];
% start = [1 2 4 8];
% end  = [1 3 7 10];


ifelse.m:
function value = ifelse(condition, if_true, if_false)
    if condition
        value = if_true;
    else
        value = if_false;                  
    end
end


main.m:
clear;  
clc;
define_constants;       

rng(0);
% ========================= 1. load a case =========================
%CaseName = 'pglib_opf_case73_ieee_rts';
%CaseName = 'case2869pegase';
%CaseName = 'case9241pegase';
CaseName = 'pglib_opf_case162_ieee_dtc';
%CaseName = 'pglib_opf_case300_ieee';
%CaseName = 'pglib_opf_case179_goc';  % unique volt_level
%CaseName = 'pglib_opf_case118_ieee';

CasePath = ['E:\MATLAB\R2018b\bin\matpower6.0\' CaseName '.m'];
% ========================= 2. waterfall technique =========================

%run('E:\MATLAB\R2018b\bin\重要专项任务\学术进展和与教授的沟通\conference paper 1\正式代码\waterfall_technique.m');
[TopPilotnode, case_WT] = waterfall(CaseName);
% ========================= 3. zone defination =========================

[table_islands, graph_nominal, digraph_nominal,...
    index_edge_withtrans, index_nodes_oftrans] = zonelabeller(case_WT);
%run('E:\MATLAB\R2018b\bin\重要专项任务\学术进展和与教授的沟通\conference paper 1\正式代码\zone_labeller.m');

% ========================= 4. zone aggregation and visualization =========================

%run('E:\MATLAB\R2018b\bin\重要专项任务\学术进展和与教授的沟通\conference paper 1\正式代码\Zonal_Aggregation_and_Visualization.m')
[digraph_zonescale_offset, digraph_zonescale_wheel, table_edge_withtrans_temp] = zonal_aggregation(table_islands, digraph_nominal, index_edge_withtrans);

% ========================= 5. current work =========================

% record each zone's nodes and edges
table_zone = digraph_zonescale_offset.Nodes;
cell_zone = table_islands;

for thiszone = 1:size(table_zone,1)
    table_zone.nodelist{thiszone} = cell_zone{thiszone,1};
    table_zone.edgelist{thiszone} = cell_zone{thiszone,2};
end


[volt_branch_ope, bool_branch_trans, table_branch_ope] = powerflow_to_table(case_WT);

table_oriedge_raw = table_branch_ope;

digraph_ori_raw = digraph_nominal;

table_oriedge_combined = edge_combination(table_oriedge_raw, 'SendingMW', 'EndNodes');

%table_oriedge_offsetted = sortrows(edge_offset(table_oriedge_combined, 3, [1 2]), [1 2]);
table_oriedge_offsetted = sortrows(table_oriedge_combined, [1 2]);
table_oriedge_offsetted = [table_oriedge_offsetted [1:size(table_oriedge_offsetted,1)]'];

table_ori_offsetted = table([table_oriedge_offsetted(:,1) table_oriedge_offsetted(:,2)], table_oriedge_offsetted(:,3), table_oriedge_offsetted(:,4), ...
    'VariableNames',["EndNodes", "SendingMW", "ori_index"]);

digraph_ori_offsetted = table_to_graph(table_ori_offsetted, case_WT, 'digraph');

%multi-graph
digraph_zonescale_multi = table_to_multigraph(table_edge_withtrans_temp, table_zone, digraph_nominal);

%wheel
to_picture_raw(digraph_ori_offsetted, CaseName);
to_picture(digraph_zonescale_offset, TopPilotnode, cell_zone, CaseName, 'force');
to_picture(digraph_zonescale_offset, TopPilotnode, cell_zone, CaseName, 'layered');
to_picture(digraph_zonescale_multi, TopPilotnode, cell_zone, CaseName, 'layered');

% to_picture(digraph_zonescale_wheel, TopPilotnode, cell_zone, CaseName)
% to_picture_wheel(digraph_zonescale_wheel, TopPilotnode, cell_zone, CaseName);

% zone_viewer_ui(digraph_ori_offsetted, digraph_zonescale_offset, digraph_zonescale_multi, ...
%     TopPilotnode, cell_zone, CaseName);
% while(1)
%     % answer = inputdlg( ...
%     %     'Please enter the Zone ID to view (enter 0 or Cancel to exit):', ...
%     %     'Zone Subnetwork Viewer', [1 50]);
% 
%     zone_first = 1;
%     zone_last = size(table_zone,1);
% 
%     prompt = { ...
%    sprintf('Available zones: from %d to %d. Enter the Zone ID (0 to exit):', zone_first, zone_last)...
%    };
% 
%     dlgTitle = 'Zone Subnetwork Viewer';
% 
%     boxSize  = [1 50];    % two lines high, 50 characters wide
% 
%     answer = inputdlg(prompt, dlgTitle, boxSize);
%     pick_zone = str2double(answer{1});
%     if 1 <= pick_zone && pick_zone <= zone_last
%         oridx_edges_ofzone = table_zone.edgelist{pick_zone};
%         idx_edges_ofzone_indigraph = find(digraph_nominal.Edges.zone == pick_zone);
%         table_edge_subgrid = digraph_nominal.Edges(idx_edges_ofzone_indigraph,:);
% 
%         idx_nodes_ofzone_indigraph = find(digraph_nominal.Nodes.zone == pick_zone);
%         table_node_subgrid = digraph_nominal.Nodes(idx_nodes_ofzone_indigraph,:);
% 
%         %allnodes = sort(unique([table_edge_subgrid.EndNodes]));
%         allnodes = sort(table_zone.nodelist{pick_zone});
%         map_allnodes = [allnodes [1:size(allnodes,1)]'];
%         table_node_subgrid.idx_subgrid = [1:length(allnodes)]';
%         EN = table_edge_subgrid.EndNodes;
%         EN_flat = EN(:);

%         [~,loc] = ismember(EN_flat, allnodes);
%         EN_new = [loc(1:length(loc)/2) loc(length(loc)/2+1:end)];
%         table_edge_subgrid.EndNodes = EN_new;
% 
%         digraph_subgrid = digraph(table_edge_subgrid);
%         figure;
%         p = plot(digraph_subgrid, Layout="force");
%     else
%         break;
%     end
% end

% pick up a zone and form its digraph
% pick_zone = 3; %this can be changed
% oridx_edges_ofzone = table_zone.edgelist{pick_zone};
% idx_edges_ofzone_indigraph = find(digraph_nominal.Edges.zone == pick_zone);
% table_edge_subgrid = digraph_nominal.Edges(idx_edges_ofzone_indigraph,:);
% 
% idx_nodes_ofzone_indigraph = find(digraph_nominal.Nodes.zone == pick_zone);
% table_node_subgrid = digraph_nominal.Nodes(idx_nodes_ofzone_indigraph,:);
% 
% allnodes = sort(unique([table_edge_subgrid.EndNodes])); 
% map_allnodes = [allnodes [1:size(allnodes,1)]'];
% table_node_subgrid.idx_subgrid = [1:size(allnodes, 1)]';
% EN = table_edge_subgrid.EndNodes;
% EN_flat = EN(:);
% [~,loc] = ismember(EN_flat, allnodes);
% EN_new = [loc(1:length(loc)/2) loc(length(loc)/2+1:end)];
% table_edge_subgrid.EndNodes = EN_new;
% 
% digraph_subgrid = digraph(table_edge_subgrid);
% figure;
% p = plot(digraph_subgrid, Layout="layered");


my_flipedge.m:
function G_new = my_flipedge(G, PowerFieldName)

%   if the power transmition < 0, then flip the from/to nodes

%   G: original digraph/graph
%   PowerFieldName: name of str of power transmition（like 'SendingMW'）
%   G_new: new digraph/graph



    EdgeTable = G.Edges;
    EndNodes = EdgeTable.EndNodes;
    PowerVal = EdgeTable.(PowerFieldName);

    idx_flip = PowerVal < 0;

    % 翻转方向
    EndNodes_flipped = EndNodes;
    EndNodes_flipped(idx_flip,1) = EndNodes(idx_flip,2);
    EndNodes_flipped(idx_flip,2) = EndNodes(idx_flip,1);

    % 功率全取正
    PowerVal_flipped = abs(PowerVal);

    % 构建新EdgeTable
    EdgeTable_new = EdgeTable;
    EdgeTable_new.EndNodes = EndNodes_flipped;
    EdgeTable_new.(PowerFieldName) = PowerVal_flipped;

    % 用原有节点表和新边表重建图
    G_new = digraph(EdgeTable_new, G.Nodes);
end


powerflow_to_table.m:
function [volt_branch_ope, bool_branch_trans, table_branch_ope] = powerflow_to_table(case_WT)

    define_constants;       

    volt_branch_ope = (case_WT.bus(case_WT.branch(:,F_BUS), BASE_KV) + ...
        case_WT.bus(case_WT.branch(:,T_BUS), BASE_KV))/2;
    
    bool_branch_trans = (case_WT.bus(case_WT.branch(:, F_BUS), BASE_KV) ~= ...
        case_WT.bus(case_WT.branch(:, T_BUS), BASE_KV));
    
    table_branch_ope = table([case_WT.branch(:, F_BUS), case_WT.branch(:, T_BUS)], [1:size(case_WT.branch,1)]' , ...
        case_WT.branch(:, RATE_A), abs(case_WT.branch(:, BR_X)), volt_branch_ope, ...
        bool_branch_trans, case_WT.branch(:, PF), ...
        'VariableNames',["EndNodes", "EdgeOrigIndex", "Cap", "BR_X", "OpVolt", "IsTrafo", "SendingMW"]);

end


sort_table.m:
function T_sorted = sort_table(T, key, varargin)


    p = inputParser;
    addRequired(p, 'T', @(x) istable(x));
    addRequired(p, 'key', @(x) iscell(x) || isnumeric(x));
    addParameter(p, 'directions', [], @(x) isempty(x) || isstring(x) || ischar(x));
    parse(p, T, varargin{:});
    
    T = p.Results.T;
    key = p.Results.key;
    directions = p.Results.directions;
    
    if isempty(key)
        error('必须指定排序字段名或列号(key)');
    end
    if isempty(directions)
        directions = 'ascend'; %默认升序
    end
    
    if iscell(key) || isnumeric(key)
        T_sorted = sortrows(T, key, directions);
    else
        error('key类型错误');
    end
end


table_to_graph.m:
function G = table_to_graph(T, case_WT, type)

    % case_WT:outcome of powerflow analysis
    % T is the table that indicate the imformation of edges

    % type is 'graph': make non-direction graph
    % type is 'diagraph': make non-direction graph

    % case to graph, only produce raw graph

    define_constants;

    p = inputParser;
    addRequired(p, 'T', @(x) istable(x));
    addRequired(p, 'case_WT', @(x) 1);
    addRequired(p, 'type', @(x) ischar(x) || isstring(x));
    parse(p, T, case_WT, type);

    type = string(lower(type));

    if strcmp(type,'graph') || strcmpi(type,"graph")
        G = graph(T);
        
    elseif strcmp(type,'digraph') || strcmpi(type,"digraph")
        G = digraph(T);

    else
        error('type must be "graph" or "digraph", got: %s', type);
    end
        G.Nodes.Volt = case_WT.bus(:, BASE_KV);
        G.Nodes.uni_index = [1:size(case_WT.bus,1)]';
        G.Nodes.matpower_index = case_WT.bus(:, BUS_I); %uni_index
        G.Nodes.angle = case_WT.bus(:, VA);


end


table_to_multigraph.m:
function multi_graph = table_to_multigraph(T, table_zone, digraph_nominal)

% from table to multigraph

% multigraph: only realize nodal aggregation, have not do edge combination yet. 

table_islands = [table_zone.nodelist table_zone.edgelist];
cell_zone = table_islands;

table_zonescale_multi = table( ...
   [T.fromzone T.tozone], ...
    T.SendingMW, ...
    (1:height(T))', ...
    'VariableNames', ["EndNodes","SendingMW","ori_index"]);

num_node_zone = [];
num_zone = size(table_islands,1);
for thiszone = 1:num_zone
    num_node_zone = [num_node_zone length(cell_zone{thiszone,1})];
end
volt_zone = [];
for thiszone = 1:num_zone
    volt_zone = [volt_zone table_zone.OpVolt(thiszone)];
end

multi_graph = digraph(table_zonescale_multi);
multi_graph.Nodes.num_node = num_node_zone';
multi_graph.Nodes.OpVolt = volt_zone';

end


to_force_layout.m:
function p = to_force_layout(G, TopPilotnode, cell_zone, CaseName)


define_constants;

% color
level_opvolt = unique(G.Nodes.OpVolt);
num_level_opvolt = length(level_opvolt);
color = color_colorbrewer(num_level_opvolt);

% width of each edge
table_zonescale = G.Edges;
P_edge_zonescale = table_zonescale.SendingMW';
P_log = log10(P_edge_zonescale + 0.1);
width_max = 3;
width_min = 1;
P_log_norm = (P_log - min(P_log))/(max(P_log) - min(P_log));
width_edge_zonescale = width_min + (width_max - width_min)*P_log_norm;

% color of each edge
color_edge_1 = [55,126,184]/255;
color_edge_2 = [228,26,28]/255;%from blue to red
color_edge = [];
for i = 1:size(table_zonescale,1)
    color_edge = [color_edge; color_edge_1 + (color_edge_2 - color_edge_1) * P_log_norm(i)];
end

%size of each supernode

num_node_eachzone = G.Nodes.num_node;
marker_min = 7;
level_size_node = ceil(log10(num_node_eachzone + 1));
size_node_zonescale = marker_min + (level_size_node - 1)*5;

[~,idx_color_eachzone] = ismember(G.Nodes.OpVolt, level_opvolt);

figure
p = plot(G, Layout='layered');
p.LineWidth = 1.0;
p.EdgeColor = 'black';
p.MarkerSize = size_node_zonescale;
p.NodeLabel = cellstr(string(G.Nodes.num_node));
%p.EdgeColor = color_edge;
color_used = [];
for thiszone = 1:size(G.Nodes,1)
    color_used = [color_used; color(idx_color_eachzone(thiszone), :)];
end
p.NodeColor = color_used;
% add a title
title(['Inter-Zonal Power Transfer - ', CaseName], 'Interpreter', 'none');

% add a legend
hold on;
h = gobjects(num_level_opvolt,1);
for k = 1:num_level_opvolt
    h(k) = scatter(NaN, NaN, 60, color(k,:), 'filled', ...
        'DisplayName', sprintf('%d kV', level_opvolt(k)));
end
hold off;
legend(h, 'Location', 'best');

Gz = G;   % easy to read
ax = gca;

% 1 scale
Z   = size(Gz.Nodes, 1);                          % num of zones
E   = size(Gz.Edges, 1);                          % num of zonal edges
if ismember('num_node', Gz.Nodes.Properties.VariableNames)
    Nbus = sum(Gz.Nodes.num_node);          % num of total nodes
else
    Nbus = NaN;
end
if ismember('OpVolt', Gz.Nodes.Properties.VariableNames)
    L = numel(unique(Gz.Nodes.OpVolt(~isnan(Gz.Nodes.OpVolt))));  % num of voltage level
else
    L = NaN;
end

% 2 indegree/outdegree
d_out = outdegree(Gz);
d_in  = indegree(Gz);
[~, k_out_deg] = max(d_out);
[~, k_in_deg ] = max(d_in);


% 3 text
info_lines = {
    sprintf('Number of zones (Z):         %d', Z)
    sprintf('Number of inter-zone edges:  %d', E)
    sprintf('Total number of nodes (N):   %s', ternum(Nbus))
    sprintf('Number of voltage levels:    %s', ternum(L))
    sprintf('Max out-degree:              Zone %d (%d edges)', k_out_deg, d_out(k_out_deg))
    sprintf('Max in-degree:               Zone %d (%d edges)', k_in_deg,  d_in(k_in_deg))
};
info_text = strjoin(info_lines, newline);

% 4 output onto picture
text(ax, 0.02, 0.98, info_text, 'Units','normalized', ...
    'VerticalAlignment','top', 'FontName','Consolas', 'FontSize',10, ...
    'BackgroundColor',[1 1 1 0.75], 'Margin',6, 'EdgeColor',[0.2 0.2 0.2]);

% mark TopPilotnode 
top_zone = -1;
for i = 1:size(cell_zone, 1)
    if ismember(TopPilotnode, cell_zone{i, 1})
        top_zone = i;
        break;
    end
end
highlight(p, [top_zone], 'Marker','s');

end

function s = ternum(val)
    if isnan(val)
        s = '-';
    else
        s = num2str(val);
    end
end


to_layered_layout.m:
function p = to_layered_layout(G, TopPilotnode, cell_zone, CaseName)

    define_constants;

    % color
    level_opvolt = unique(G.Nodes.OpVolt);
    num_level_opvolt = length(level_opvolt);
    color = color_colorbrewer(num_level_opvolt);
    
    % width of each edge
    table_zonescale = G.Edges;
    P_edge_zonescale = table_zonescale.SendingMW';
    P_log = log10(P_edge_zonescale + 0.1);
    width_max = 3;
    width_min = 1;
    P_log_norm = (P_log - min(P_log))/(max(P_log) - min(P_log));
    width_edge_zonescale = width_min + (width_max - width_min)*P_log_norm;
    
    % color of each edge
    color_edge_1 = [55,126,184]/255;
    color_edge_2 = [228,26,28]/255;%from blue to red
    color_edge = [];
    for i = 1:size(table_zonescale,1)
        color_edge = [color_edge; color_edge_1 + (color_edge_2 - color_edge_1) * P_log_norm(i)];
    end
    
    %size of each supernode
    
    num_node_eachzone = G.Nodes.num_node;
    marker_min = 7;
    level_size_node = ceil(log10(num_node_eachzone + 1));
    size_node_zonescale = marker_min + (level_size_node - 1)*5;
    
    [~,idx_color_eachzone] = ismember(G.Nodes.OpVolt, level_opvolt);

    unique_voltages = sort(unique(G.Nodes.OpVolt), 'descend');
    [~, voltage_level] = ismember(G.Nodes.OpVolt, unique_voltages);

    % value of each layer（such as: 1, 2, 3,...）
    Y_layer = voltage_level;

    % get X data by plotting force layout 
    % p_temp = plot(G, 'Layout', "force", 'visible', 'off');
    % X = p_temp.XData;

    % plot
    figure
    p = plot(G);
    p.YData = -Y_layer;  % higher voltage level zones will be above lower ones
    p.XData = get_X_data(G, -Y_layer);

    p.LineWidth = 1.0;
    p.EdgeColor = 'black';
    p.MarkerSize = size_node_zonescale;
    p.NodeLabel = cellstr(string(G.Nodes.num_node));
    %p.EdgeColor = color_edge;
    color_used = [];
    for thiszone = 1:size(G.Nodes,1)
        color_used = [color_used; color(idx_color_eachzone(thiszone), :)];
    end
    p.NodeColor = color_used;
    % add a title
    title(['Inter-Zonal Power Transfer - ', CaseName], 'Interpreter', 'none');
    
    % add a legend
    hold on;
    h = gobjects(num_level_opvolt,1);
    for k = 1:num_level_opvolt
        h(k) = scatter(NaN, NaN, 60, color(k,:), 'filled', ...
            'DisplayName', sprintf('%d kV', level_opvolt(k)));
    end
    hold off;
    legend(h, 'Location', 'best');
    
    Gz = G;   % easy to read
    ax = gca;
    
    % 1 scale
    Z   = size(Gz.Nodes, 1);                          % num of zones
    E   = size(Gz.Edges, 1);                          % num of zonal edges
    if ismember('num_node', Gz.Nodes.Properties.VariableNames)
        Nbus = sum(Gz.Nodes.num_node);          % num of total nodes
    else
        Nbus = NaN;
    end
    if ismember('OpVolt', Gz.Nodes.Properties.VariableNames)
        L = numel(unique(Gz.Nodes.OpVolt(~isnan(Gz.Nodes.OpVolt))));  % num of voltage level
    else
        L = NaN;
    end
    
    % 2 indegree/outdegree
    d_out = outdegree(Gz);
    d_in  = indegree(Gz);
    [~, k_out_deg] = max(d_out);
    [~, k_in_deg ] = max(d_in);
    
    
    % 3 text
    info_lines = {
        sprintf('Number of zones (Z):         %d', Z)
        sprintf('Number of inter-zone edges:  %d', E)
        sprintf('Total number of nodes (N):   %s', ternum(Nbus))
        sprintf('Number of voltage levels:    %s', ternum(L))
        sprintf('Max out-degree:              Zone %d (%d edges)', k_out_deg, d_out(k_out_deg))
        sprintf('Max in-degree:               Zone %d (%d edges)', k_in_deg,  d_in(k_in_deg))
    };
    info_text = strjoin(info_lines, newline);
    
    % 4 output onto picture
    text(ax, 0.02, 0.98, info_text, 'Units','normalized', ...
        'VerticalAlignment','top', 'FontName','Consolas', 'FontSize',10, ...
        'BackgroundColor',[1 1 1 0.75], 'Margin',6, 'EdgeColor',[0.2 0.2 0.2]);
    
    % mark TopPilotnode 
    top_zone = -1;
    for i = 1:size(cell_zone, 1)
        if ismember(TopPilotnode, cell_zone{i, 1})
            top_zone = i;
            break;
        end
    end
    highlight(p, [top_zone], 'Marker','s');

end

function s = ternum(val)
    if isnan(val)
        s = '-';
    else
        s = num2str(val);
    end
end


to_picture.m:
function p = to_picture(digraph_zonescale_offset, TopPilotnode, cell_zone, CaseName, layout)
    
% to make aggregated graph

% the most classic graph showing the big picture of power transmission
% among zones

if nargin < 5 || isempty(layout)
        layout = 'force';
    end
    layout = char(layout);  

define_constants;

% color
level_opvolt = unique(digraph_zonescale_offset.Nodes.OpVolt);
num_level_opvolt = length(level_opvolt);
color = color_colorbrewer(num_level_opvolt);

% width of each edge
table_zonescale = digraph_zonescale_offset.Edges;
P_edge_zonescale = table_zonescale.SendingMW';
P_log = log10(P_edge_zonescale + 0.1);
width_max = 3;
width_min = 1;
P_log_norm = (P_log - min(P_log))/(max(P_log) - min(P_log));
width_edge_zonescale = width_min + (width_max - width_min)*P_log_norm;

% color of each edge
color_edge_1 = [55,126,184]/255;
color_edge_2 = [228,26,28]/255;%from blue to red
color_edge = [];
for i = 1:size(table_zonescale,1)
    color_edge = [color_edge; color_edge_1 + (color_edge_2 - color_edge_1) * P_log_norm(i)];
end

%size of each supernode
% num_node_eachzone = digraph_zonescale_offset.Nodes.num_node;
% num_log = log10(num_node_eachzone + 1);
% marker_min = 5;
% marker_max = 20;
% num_log_norm = (num_log - min(num_log))/(max(num_log) - min(num_log));
% size_node_zonescale = marker_min + (marker_max - marker_min)*num_log_norm;

num_node_eachzone = digraph_zonescale_offset.Nodes.num_node;
marker_min = 7;
level_size_node = ceil(log10(num_node_eachzone + 1));
size_node_zonescale = marker_min + (level_size_node - 1)*5;

[~,idx_color_eachzone] = ismember(digraph_zonescale_offset.Nodes.OpVolt, level_opvolt);

if strcmp(layout, 'force')

    p = to_force_layout(digraph_zonescale_offset, TopPilotnode, cell_zone, CaseName);

else % layout is 'layered'
    
    p = to_layered_layout(digraph_zonescale_offset, TopPilotnode, cell_zone, CaseName);

end

end


function s = ternum(val)
    if isnan(val)
        s = '-';
    else
        s = num2str(val);
    end
end


to_picture_raw.m:
function p = to_picture_raw(digraph, CaseName)
    
    % to make raw graph

    % A raw graph is made from raw data which exists before combination

    G = digraph;
    
    if ismember('Volt', G.Nodes.Properties.VariableNames)
        num_color = length(unique(digraph.Nodes.Volt));
        num_level_opvolt = num_color;
        level_volt = unique(digraph.Nodes.Volt);
        char_volt = 'Volt';
    else
        num_color = length(unique(digraph.Nodes.OpVolt));
        num_level_opvolt = num_color;
        level_volt = unique(digraph.Nodes.OpVolt);
        char_volt = 'OpVolt';
    end

    
    color = color_colorbrewer(num_color);

    %% parameters

    marksize = 5;
       
    figure
    p = plot(G, Layout="force");
    % p.LineWidth = width_edge_zonescale;
    p.LineWidth = 1.0;
    p.EdgeColor = 'black';
    p.MarkerSize = marksize;

    color_to_use = [];
    [~,idx_color_eachnode] = ismember(G.Nodes.(char_volt), level_volt);
    for i = 1:size(G.Nodes,1)
        color_to_use = [color_to_use; color(idx_color_eachnode(i), :)];
    end
    p.NodeColor = color_to_use;
    % add a title
    title(['original distribution network - ', CaseName], 'Interpreter', 'none');
    
    % add a legend
    hold on;
    h = gobjects(num_level_opvolt,1);
    for k = 1:num_level_opvolt
        h(k) = scatter(NaN, NaN, 60, color(k,:), 'filled', ...
            'DisplayName', sprintf('%d kV', level_volt(k)));
    end
    hold off;
    legend(h, 'Location', 'best');
    
    edge_trans = [];
    for i = 1:size(G.Edges, 1)
        if G.Nodes.Volt(G.Edges.EndNodes(i, 1)) ~= G.Nodes.Volt(G.Edges.EndNodes(i, 2))
            edge_trans = [edge_trans i];
        end
    end
    highlight(p, 'Edges', edge_trans, 'LineWidth', 1.5, 'LineStyle','-.');

    Gz = digraph;   % easy to read
    ax = gca;
    
    % 1 scale
    Z   = size(Gz.Nodes, 1);                          % number of nodes
    E   = size(Gz.Edges, 1);                          % number of edges
    L = numel(unique(Gz.Nodes.Volt));

    % 2 indegree/outdegree
    d_out = outdegree(Gz);
    d_in  = indegree(Gz);
    [~, k_out_deg] = max(d_out);
    [~, k_in_deg ] = max(d_in);
    
    
    % 3 text
    info_lines = {
        sprintf('Number of nodes:             %d', Z)
        sprintf('Number of edges:             %d', E)
        sprintf('Number of voltage levels:    %s', ternum(L))
        sprintf('Max out-degree:              node %d (%d edges)', k_out_deg, d_out(k_out_deg))
        sprintf('Max in-degree:               node %d (%d edges)', k_in_deg,  d_in(k_in_deg))
    };
    info_text = strjoin(info_lines, newline);
    
    % 4 output onto picture）
    text(ax, 0.02, 0.98, info_text, 'Units','normalized', ...
        'VerticalAlignment','top', 'FontName','Consolas', 'FontSize',10, ...
        'BackgroundColor',[1 1 1 0.75], 'Margin',6, 'EdgeColor',[0.2 0.2 0.2]);
    

    if_correct = 1;

end

function s = ternum(val)
    if isnan(val)
        s = '-';
    else
        s = num2str(val);
    end
end


to_picture_wheel.m:
function p = to_picture_wheel(G, TopPilotnode, cell_zone, CaseName)

    % to make wheel graph

    % A raw graph is made form the data after combination and before
    % aggregation 

    % Delete single edges between two zones to show wheels

define_constants;

% color
level_opvolt = unique(G.Nodes.OpVolt);
num_level_opvolt = length(level_opvolt);
color = color_colorbrewer(num_level_opvolt);

% width of each edge
table_zonescale = G.Edges;
P_edge_zonescale = table_zonescale.SendingMW';
P_log = log10(P_edge_zonescale + 0.1);
width_max = 3;
width_min = 1;
P_log_norm = (P_log - min(P_log))/(max(P_log) - min(P_log));
width_edge_zonescale = width_min + (width_max - width_min)*P_log_norm;

% color of each edge
color_edge_1 = [55,126,184]/255;
color_edge_2 = [228,26,28]/255;%from blue to red
color_edge = [];
for i = 1:size(table_zonescale,1)
    color_edge = [color_edge; color_edge_1 + (color_edge_2 - color_edge_1) * P_log_norm(i)];
end

%size of each supernode
% num_node_eachzone = G.Nodes.num_node;
% num_log = log10(num_node_eachzone + 1);
% marker_min = 5;
% marker_max = 20;
% num_log_norm = (num_log - min(num_log))/(max(num_log) - min(num_log));
% size_node_zonescale = marker_min + (marker_max - marker_min)*num_log_norm;

num_node_eachzone = G.Nodes.num_node;
marker_min = 7;
level_size_node = ceil(log10(num_node_eachzone + 1));
size_node_zonescale = marker_min + (level_size_node - 1)*5;

%color of supernode
% color_zone_1 = [204,235,197]/255;
% color_zone_2 = [77,175,74]/255;%from light green to dark green
% color_zone = [];
% for i = 1:size(G.Nodes,1)
%     color_zone = [color_zone; color_zone_1 + (color_zone_2 - color_zone_1) * num_log_norm(i)];
% end

% color_idx = floor(linspace(1,size(color, 1), num_level_opvolt));
% color = color(color_idx,:);
[~,idx_color_eachzone] = ismember(G.Nodes.OpVolt, level_opvolt);


%% output diagraph
figure
%p = plot(G, Layout="layered");
p = plot(G, Layout="force");
% p.LineWidth = width_edge_zonescale;
% linewidth
p.LineWidth = 1.0;
idx_nonwheel = find_single_edge(table2array(G.Edges), [1 2]);
if ~isempty(idx_nonwheel)
    highlight(p, 'Edges', idx_nonwheel, 'LineWidth', 0.001);
    highlight(p, 'Edges', idx_nonwheel, 'EdgeColor', [255 255 255]/255);
    %highlight(p, 'Edges', idx_nonwheel, 'EdgeAlpha', 0.12);
    highlight(p, 'Edges', idx_nonwheel, 'ArrowSize', 5);
end
%p.EdgeColor = 'black';
%p.MarkerSize = size_node_zonescale;
p.MarkerSize = size_node_zonescale;
p.NodeLabel = cellstr(string(G.Nodes.num_node));
%p.EdgeColor = color_edge;
color_used = [];
for thiszone = 1:size(G.Nodes,1)
    color_used = [color_used; color(idx_color_eachzone(thiszone), :)];
end
p.NodeColor = color_used;
% add a title
title(['Inter-Zonal Power Transfer - ', CaseName], 'Interpreter', 'none');

% add a legend
hold on;
h = gobjects(num_level_opvolt,1);
for k = 1:num_level_opvolt
    h(k) = scatter(NaN, NaN, 60, color(k,:), 'filled', ...
        'DisplayName', sprintf('%d kV', level_opvolt(k)));
end
hold off;
legend(h, 'Location', 'best');

Gz = G;   % easy to read
ax = gca;

% 1 scale
Z   = size(Gz.Nodes, 1);                          % num of zones
E   = size(Gz.Edges, 1);                          % num of zonal edges
if ismember('num_node', Gz.Nodes.Properties.VariableNames)
    Nbus = sum(Gz.Nodes.num_node);          % num of total nodes
else
    Nbus = NaN;
end
if ismember('OpVolt', Gz.Nodes.Properties.VariableNames)
    L = numel(unique(Gz.Nodes.OpVolt(~isnan(Gz.Nodes.OpVolt))));  % num of volt levels
else
    L = NaN;
end

% 2 indegree/outdegree
d_out = outdegree(Gz);
d_in  = indegree(Gz);
[~, k_out_deg] = max(d_out);
[~, k_in_deg ] = max(d_in);

% 3 Number of powerflow wheels
num_wheel = (size(G.Edges, 1) - length(idx_nonwheel))/2;


% 3 text
info_lines = {
    sprintf('Number of zones (Z):         %d', Z)
    sprintf('Number of wheels (Z):        %d', num_wheel)
    sprintf('Total number of nodes (N):   %s', ternum(Nbus))
    sprintf('Number of voltage levels:    %s', ternum(L))
    

};
info_text = strjoin(info_lines, newline);

% 4 output onto picture
text(ax, 0.02, 0.98, info_text, 'Units','normalized', ...
    'VerticalAlignment','top', 'FontName','Consolas', 'FontSize',10, ...
    'BackgroundColor',[1 1 1 0.75], 'Margin',6, 'EdgeColor',[0.2 0.2 0.2]);

% mark TopPilotnode 
top_zone = -1;
for i = 1:size(cell_zone, 1)
    if ismember(TopPilotnode, cell_zone{i, 1})
        top_zone = i;
        break;
    end
end
highlight(p, [top_zone], 'Marker','s');

end


function s = ternum(val)
    if isnan(val)
        s = '-';
    else
        s = num2str(val);
    end

end


waterfall.m:
function [TopPilotnode, case_WT] = waterfall(CaseName)

P_inject = 1;

define_constants;

%% ========================= 1. Data Preparation and Network Construction =========================
%CaseName = 'pglib_opf_case73_ieee_rts';

BaseCase = loadcase(CaseName);   

TestCase = ext2int(runpf(BaseCase)); % runpf: MATPOWER

OperatingVoltages = (TestCase.bus(TestCase.branch(:, F_BUS), BASE_KV) + TestCase.bus(TestCase.branch(:, T_BUS), BASE_KV))/2;

IsEdgeTransformer = (TestCase.bus(TestCase.branch(:, F_BUS), BASE_KV) ~= TestCase.bus(TestCase.branch(:, T_BUS), BASE_KV));

SysEdgeTableCapWeights = table([TestCase.branch(:, F_BUS), TestCase.branch(:, T_BUS)], ... % 1. 两端节点编号，N×2数组
    [1:size(TestCase.branch, 1)]', ... 
    TestCase.branch(:, RATE_A), ... 
    abs(TestCase.branch(:, BR_X)), ... 
    OperatingVoltages, ... 
    IsEdgeTransformer, ... 
    'VariableNames',["EndNodes", "EdgeIndex", "Cap", "Weight", "OpVolt", "IsTrafo"]);

NominalGraph = graph(SysEdgeTableCapWeights);

NominalGraph.Nodes.Voltages = TestCase.bus(:, BASE_KV);

NominalDiGraph = digraph(SysEdgeTableCapWeights);
NominalDiGraph.Nodes.Voltages = TestCase.bus(:, BASE_KV);

%% ========================= 2. Network Topological Analysis =========================

GeodesicMatrix = distances(NominalGraph,'Method','unweighted');

HighestVoltagesNodes = find(NominalGraph.Nodes.Voltages == max(NominalGraph.Nodes.Voltages));
LowestVoltagesNodes = find(NominalGraph.Nodes.Voltages == min(NominalGraph.Nodes.Voltages));

NodeRemoteness = sum(GeodesicMatrix,1);

RemotestNode = find(NodeRemoteness == max(NodeRemoteness));
RemotestNode = RemotestNode(1);

TopPilotnode = HighestVoltagesNodes(find(NodeRemoteness(HighestVoltagesNodes) == max(NodeRemoteness(HighestVoltagesNodes))));
TopPilotnode = TopPilotnode(1);

NodalDistanceFromPilot = GeodesicMatrix(:, TopPilotnode);

BottomPilotNode = LowestVoltagesNodes(find(NodalDistanceFromPilot(LowestVoltagesNodes) == min(NodalDistanceFromPilot(LowestVoltagesNodes))));
BottomPilotNode = BottomPilotNode(1);

%% ========================= 3. PTDF Matrix Analysis and Minimum Set Cover =========================

CurPTDFMatrix = makePTDF(TestCase.baseMVA, TestCase.bus, TestCase.branch, TopPilotnode); 

binaryMatrix = (abs(CurPTDFMatrix) > 0.001);


coveredRows = false(size(binaryMatrix, 1), 1);%Rows:branch    Clm:bus
selectedCols = [];   
counter = 0;
% Greedy method
while ~all(coveredRows)
    coverage = sum(binaryMatrix(~coveredRows, :), 1);    
    [cov_max, bestCol] = max(coverage); 
    if cov_max == 0
        warning('No further branches can be covered by any bus. Terminating selection.');
        break;
    end
    selectedCols = [selectedCols; bestCol];                
    coveredRows = coveredRows | binaryMatrix(:, bestCol); 
    counter = counter + 1;
end

%% ========================= 4. Synthetic Flow Scenario Construction and Edge Direction Processing =========================

FlowTestCase = BaseCase;                
FlowTestCase.bus(:, PD) = 0;            
FlowTestCase.gen(:, PG) = 0;

FlowTestCase.bus(TopPilotnode, PD) = -P_inject * size(selectedCols,1);  
FlowTestCase.bus(selectedCols, PD) = P_inject;                         

case_WT = ext2int(rundcpf(FlowTestCase)); % case of waterfall
end


waterfall_technique.m:
%% ========================= 0.Preparation =========================
P_inject = 1; 
%clear;                  
%define_constants;       

%rng(0);
%% ========================= 1. Data Preparation and Network Construction =========================
%CaseName = 'pglib_opf_case73_ieee_rts';

BaseCase = loadcase(CaseName);   

% pglib_opf_case14_ieeeprac 
% pglib_opf_case30_as
% pglib_opf_case30_ieee
% ...

TestCase = ext2int(runpf(BaseCase)); % runpf: MATPOWER

OperatingVoltages = (TestCase.bus(TestCase.branch(:, F_BUS), BASE_KV) + TestCase.bus(TestCase.branch(:, T_BUS), BASE_KV))/2;

IsEdgeTransformer = (TestCase.bus(TestCase.branch(:, F_BUS), BASE_KV) ~= TestCase.bus(TestCase.branch(:, T_BUS), BASE_KV));

SysEdgeTableCapWeights = table([TestCase.branch(:, F_BUS), TestCase.branch(:, T_BUS)], ... % 1. 两端节点编号，N×2数组
    [1:size(TestCase.branch, 1)]', ... 
    TestCase.branch(:, RATE_A), ... 
    abs(TestCase.branch(:, BR_X)), ... 
    OperatingVoltages, ... 
    IsEdgeTransformer, ... 
    'VariableNames',["EndNodes", "EdgeIndex", "Cap", "Weight", "OpVolt", "IsTrafo"]);

NominalGraph = graph(SysEdgeTableCapWeights);

NominalGraph.Nodes.Voltages = TestCase.bus(:, BASE_KV);

NominalDiGraph = digraph(SysEdgeTableCapWeights);
NominalDiGraph.Nodes.Voltages = TestCase.bus(:, BASE_KV);

%% ========================= 2. Network Topological Analysis =========================
tic
GeodesicMatrix = distances(NominalGraph,'Method','unweighted');
toc

HighestVoltagesNodes = find(NominalGraph.Nodes.Voltages == max(NominalGraph.Nodes.Voltages));
LowestVoltagesNodes = find(NominalGraph.Nodes.Voltages == min(NominalGraph.Nodes.Voltages));

NodeRemoteness = sum(GeodesicMatrix,1);

RemotestNode = find(NodeRemoteness == max(NodeRemoteness));
RemotestNode = RemotestNode(1);

TopPilotnode = HighestVoltagesNodes(find(NodeRemoteness(HighestVoltagesNodes) == max(NodeRemoteness(HighestVoltagesNodes))));
TopPilotnode = TopPilotnode(1);

NodalDistanceFromPilot = GeodesicMatrix(:, TopPilotnode);

BottomPilotNode = LowestVoltagesNodes(find(NodalDistanceFromPilot(LowestVoltagesNodes) == min(NodalDistanceFromPilot(LowestVoltagesNodes))));
BottomPilotNode = BottomPilotNode(1);

%% ========================= 3. PTDF Matrix Analysis and Minimum Set Cover =========================

CurPTDFMatrix = makePTDF(TestCase.baseMVA, TestCase.bus, TestCase.branch, TopPilotnode); 

binaryMatrix = (abs(CurPTDFMatrix) > 0.001);


coveredRows = false(size(binaryMatrix, 1), 1);%Rows:branch    Clm:bus
selectedCols = [];   
counter = 0;
% Greedy method
while ~all(coveredRows)
    coverage = sum(binaryMatrix(~coveredRows, :), 1);    
    [cov_max, bestCol] = max(coverage); 
    if cov_max == 0
        warning('No further branches can be covered by any bus. Terminating selection.');
        break;
    end
    selectedCols = [selectedCols; bestCol];                
    coveredRows = coveredRows | binaryMatrix(:, bestCol); 
    counter = counter + 1;
end

%% ========================= 4. Synthetic Flow Scenario Construction and Edge Direction Processing =========================

FlowTestCase = BaseCase;                
FlowTestCase.bus(:, PD) = 0;            
FlowTestCase.gen(:, PG) = 0;



FlowTestCase.bus(TopPilotnode, PD) = -P_inject * size(selectedCols,1);  
FlowTestCase.bus(selectedCols, PD) = P_inject;                         

%case_WT = ext2int(runpf(FlowTestCase)); % case of waterfall
case_WT = ext2int(rundcpf(FlowTestCase)); % case of waterfall

% num_step_max = 10;
% mpopt = mpoption('pf.nr.max_it', num_step_max); 
% [case_WT, success] = runpf(FlowTestCase, mpopt);
% while(~success)
%     num_step_max = num_step_max + 20;
%     mpopt = mpoption('pf.nr.max_it', num_step_max);
%     [case_WT, success] = runpf(FlowTestCase, mpopt);
% end
% case_WT = ext2int(case_WT);


zonal_aggregation.m:

function [digraph_zonescale_offset, digraph_zonescale_wheel, table_edge_withtrans_temp] = zonal_aggregation(table_islands, digraph_nominal, ...
    index_edge_withtrans)

num_zone = size(table_islands,1);
cell_zone = table_islands;

num_node_zone = [];
for thiszone = 1:num_zone
    num_node_zone = [num_node_zone length(cell_zone{thiszone,1})];
end

num_edge_zone = [];
for thiszone = 1:num_zone
    num_edge_zone = [num_edge_zone length(cell_zone{thiszone,2})];
end
num_edge_withtrans = length(index_edge_withtrans);

volt_zone = [];
for thiszone = 1:num_zone
    volt_zone = [volt_zone digraph_nominal.Nodes.Volt(cell_zone{thiszone,1}(1))];
end

table_node = sortrows(digraph_nominal.Nodes, "uni_index");  %sorted
table_edge = sortrows(digraph_nominal.Edges, "EdgeOrigIndex");  %sorted

table_edge_withtrans = table_edge(find(table_edge.IsTrafo == 1), :) ;
node_edgewithtrans_from = table_edge_withtrans.EndNodes(:,1);
node_edgewithtrans_to = table_edge_withtrans.EndNodes(:,2);
table_edge_withtrans.fromzone = table_node.zone(node_edgewithtrans_from);
table_edge_withtrans.tozone = table_node.zone(node_edgewithtrans_to);


table_edge_withtrans_temp = sort_table(table_edge_withtrans,'key', {'fromzone', 'tozone'}, 'directions', 'ascend');% multi-graph
edge_forzone_combined = edge_combination(table_edge_withtrans_temp, 'SendingMW', {'fromzone','tozone'});
edge_forzone_combined = sortrows(edge_forzone_combined, [1 2]);

edge_forzone_offsetted = edge_offset(edge_forzone_combined, 3, [1 2]);
edge_forzone_offsetted = sortrows(edge_forzone_offsetted, [1 2]);

edge_forzone = [edge_forzone_offsetted [1:size(edge_forzone_offsetted,1)]'];

table_zonescale = table([edge_forzone(:,1), edge_forzone(:,2)], [edge_forzone(:,3)], edge_forzone(:,4),...
    'VariableNames',["EndNodes", "SendingMW", "ori_index"]);
table_zonescale = sortrows(table_zonescale, "ori_index");
num_trans = size(table_edge_withtrans_temp,1);

digraph_zonescale_offset = digraph(table_zonescale);
digraph_zonescale_offset.Nodes.num_node = num_node_zone';
digraph_zonescale_offset.Nodes.OpVolt = volt_zone';

%wheel_1 after combination before aggregation
edge_forzone_wheel = [edge_forzone_combined [1:size(edge_forzone_combined,1)]'];
table_zonescale_wheel = table([edge_forzone_wheel(:,1), edge_forzone_wheel(:,2)], [edge_forzone_wheel(:,3)], edge_forzone_wheel(:,4),...
    'VariableNames', ["EndNodes", "SendingMW", "ori_index"]);
table_zonescale_wheel = sortrows(table_zonescale_wheel, "ori_index");

digraph_zonescale_wheel = digraph(table_zonescale_wheel);
digraph_zonescale_wheel.Nodes.num_node = num_node_zone';
digraph_zonescale_wheel.Nodes.OpVolt = volt_zone';


end


zonelabeller.m:
function [table_islands, graph_nominal,...
    digraph_nominal, index_edge_withtrans, index_nodes_oftrans] = zonelabeller(case_WT)

    define_constants;       

    [volt_branch_ope, bool_branch_trans, table_branch_ope] = powerflow_to_table(case_WT);

    graph_nominal = case_to_graph(case_WT, 'graph');

    digraph_nominal = case_to_graph(case_WT, 'digraph');
    digraph_nominal = my_flipedge(digraph_nominal, 'SendingMW');
    
    
    index_edge_withtrans = find(graph_nominal.Edges.IsTrafo == 1); %index
    index_nodes_oftrans = graph_nominal.Edges.EndNodes(index_edge_withtrans,:);
    
    graph_nominal_notrans = rmedge(graph_nominal, index_edge_withtrans);
    
    
    
    index_comp = conncomp(graph_nominal_notrans); 
    num_comp = max(index_comp); 
    
    table_islands = cell(num_comp,2);
    
    for k = 1:num_comp
        table_islands{k, 1} = [find(index_comp == k)]; %nodes, ori_index
        for thisedge = 1:size(graph_nominal_notrans.Edges,1)
            if ismember(graph_nominal_notrans.Edges.EndNodes(thisedge,1),table_islands{k, 1}) && ...
                    ismember(graph_nominal_notrans.Edges.EndNodes(thisedge,2),table_islands{k, 1})
                table_islands{k, 2} = [table_islands{k, 2} graph_nominal_notrans.Edges.EdgeOrigIndex(thisedge)];% edge,ori_index
            end
        end
    end
    
    
    num_island = max(index_comp);
    
    index_edge_island = cell(num_island,1);
    for k = 1:num_island
        index_edge_island{k} = find(ismember(digraph_nominal.Edges.EdgeOrigIndex, table_islands{k,2}));
    end
    
    orindex_edge_withtrans = graph_nominal.Edges.EdgeOrigIndex(index_edge_withtrans);
    index_edge_withtrans = find(ismember(digraph_nominal.Edges.EdgeOrigIndex, orindex_edge_withtrans));
    
    %% 以下代码标注每个节点和branch的zone。变压器的zone 为0 
    % The following code labels the zone of each node and branch. The zone of transformers is set to 0.
    
    % digraph_nominal.Nodes.matpower_index = case_WT.bus(:, BUS_I);
    % digraph_nominal.Nodes.angle = case_WT.bus(:, VA);
    
    % zone for each node 
    for thiszone = 1:num_island
        digraph_nominal.Nodes.zone(table_islands{thiszone,1}) = thiszone;
    end
    
    %zone for each edge 
    for thiszone = 1:num_island
        index_edge_thiszone = find(ismember(digraph_nominal.Edges.EdgeOrigIndex, table_islands{thiszone, 2})); %index of given zone in digraph 
        digraph_nominal.Edges.zone(index_edge_thiszone) = thiszone;
    end


end
