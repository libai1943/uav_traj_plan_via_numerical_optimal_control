function RunMe(caseName,showFigures,runMPC)
% RUNME One-click fixed-wing trajectory-planning demonstration.
% Requires MATLAB R2021b or newer and Optimization Toolbox.
% RunMe; RunMe('two_columns'); RunMe('all',false) for a headless run.
if nargin<1,caseName='all';end
if nargin<2,showFigures=true;end
if nargin<3,runMPC=true;end
root=fileparts(mfilename('fullpath'));oldPath=path;cleanup=onCleanup(@()path(oldPath)); %#ok<NASGU>
addpath(fullfile(root,'matlab'));
assert(exist('fmincon','file')==2 && license('test','Optimization_Toolbox'),'Optimization Toolbox with fmincon is required.');
cfg=jsondecode(fileread(fullfile(root,'scenes.json')));
assert(strcmp(caseName,'all')||any(strcmp(caseName,{cfg.cases.name})),'Unknown case name.');
out=fullfile(root,'outputs','matlab');if ~isfolder(out),mkdir(out);end
results=cell(0,1);
for i=1:numel(cfg.cases)
 c=cfg.cases(i);if ~strcmp(caseName,'all') && ~strcmp(caseName,c.name),continue;end
 fprintf('\n=== %s: 3-D A* -> cubic spline -> corridors -> NLP ===\n',c.name);
 result=plan_fixedwing(cfg,c);
 folder=fullfile(out,c.name);if ~isfolder(folder),mkdir(folder);end
 fid=fopen(fullfile(folder,'result.json'),'w');assert(fid>=0,'Cannot open result file.');
 fwrite(fid,jsonencode(result),'char');fclose(fid);save(fullfile(folder,'result.mat'),'result');
 plot_result(result,cfg,folder,showFigures);disp(result.checks);
 if strcmp(c.name,'two_columns') && runMPC
  disp('Tracking supplement: local linear model -> QP -> first input -> feedback');
  tracking=track_mpc(result,cfg);fid=fopen(fullfile(folder,'mpc_tracking.json'),'w');
  fwrite(fid,jsonencode(tracking),'char');fclose(fid);plot_tracking(tracking,folder,showFigures);
  fprintf('MPC/open-loop RMS position error: %.3f / %.3f m\n',tracking.rms_mpc_position_error_m,tracking.rms_open_loop_position_error_m);
 end
 results{end+1}=struct('name',result.name,'objective',result.objective,'iterations',result.iterations,...
  'seconds',result.seconds,'checks',result.checks); %#ok<AGROW>
end
report=struct('matlab',version,'config',cfg,'results',[results{:}]);
fid=fopen(fullfile(out,'summary.json'),'w');fwrite(fid,jsonencode(report),'char');fclose(fid);
fprintf('\nAll selected cases passed. Figures/data saved to %s\n',out);
end
