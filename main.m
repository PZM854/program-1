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