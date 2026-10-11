function result=plan_fixedwing(cfg,c)
% Scaled direct transcription; analytic derivatives; local SQP with fmincon.
N=cfg.intervals;M=N+1;dt=cfg.duration_s/N;obs=c.columns;
sx=[600;150;150;20;.3;.6;.6];su=[2;1;.4];scale=[repmat(sx,M,1);repmat(su,M,1)];
xs=cfg.start_state(:);xf=c.goal_state(:);
[raw,expanded]=astar3d(xs(1:3)',xf(1:3)',obs,cfg);
[seed,splinePath,knots,splineClearance]=spline_seed(raw,obs,cfg);
[boxes,owner,plo,phi]=build_corridors(seed,obs,cfg);
X=zeros(7,M);X(1:3,:)=seed';dp=zeros(3,M);
for d=1:3,dp(d,:)=gradient(X(d,:),dt);end
X(4,:)=vecnorm(dp);X(5,:)=atan2(dp(3,:),vecnorm(dp(1:2,:)));
X(6,:)=unwrap(atan2(dp(2,:),dp(1,:)));X(7,:)=atan(X(4,:)/9.81.*gradient(X(6,:),dt));
X(:,1)=xs;X(:,end)=xf;
U=[gradient(X(4,:),dt)+9.81*sin(X(5,:));cos(X(5,:))./cos(X(7,:));gradient(X(7,:),dt)];
lx=repmat(cfg.state_lower(:),1,M);hx=repmat(cfg.state_upper(:),1,M);
lx(1:3,:)=plo';hx(1:3,:)=phi';
assert(all(xs>=lx(:,1))&&all(xs<=hx(:,1))&&all(xf>=lx(:,end))&&all(xf<=hx(:,end)),'Endpoint outside its corridor/physical bounds.');
lx(:,1)=xs;hx(:,1)=xs;lx(:,end)=xf;hx(:,end)=xf;
lb=[lx(:);repmat(cfg.input_lower(:),M,1)]./scale;ub=[hx(:);repmat(cfg.input_upper(:),M,1)]./scale;
z0=min(max([X(:);U(:)]./scale,lb),ub);[initialX,initialU]=unpack(z0);
W=ones(1,M);W([1,end])=.5;history=[];
% Central-difference check is part of the one-click run, not a stored assertion.
[~,g0]=objective(z0);[~,~,G0,E0]=constraints(z0);derr=0;
for j=round(linspace(8,numel(z0)-3,18))
 dz=zeros(size(z0));dz(j)=1e-6;[cp,ep]=constraints(z0+dz);[cm,em]=constraints(z0-dz);
 derr=max([derr;abs((objective(z0+dz)-objective(z0-dz))/2e-6-g0(j));...
  norm((cp-cm)/2e-6-G0(j,:)',inf);norm((ep-em)/2e-6-E0(j,:)',inf)]);
end
assert(derr<1e-5,'Analytic derivative check failed.');
fprintf('A*: %d points; %d expanded nodes; derivative error %.3g\n',size(raw,1),expanded,derr);
options=optimoptions('fmincon','Algorithm','sqp','Display','final',...
 'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'MaxIterations',350,...
 'MaxFunctionEvaluations',150000,'ConstraintTolerance',2e-7,'OptimalityTolerance',1e-6,'OutputFcn',@record);
tic;[z,J,flag,stats]=fmincon(@objective,z0,[],[],[],[],lb,ub,@constraints,options);seconds=toc;
assert(flag>0,'fmincon did not converge; successful results will not be exported.');
[X,U]=unpack(z);[cc,ee]=constraints(z);residual=max([0;cc;abs(ee);lb-z;z-ub]);
% Independent interval RK4 checks and a separate whole-horizon open-loop rollout.
dense=zeros(N*21,8);rollout=zeros(N*20+1,8);rollout(1,:)=[0 X(:,1)'];
defect=0;violation=0;forward=X(:,1);
for k=1:N
 local=X(:,k);box=boxes(owner(k),:);
 for j=0:20
  dense((k-1)*21+j+1,:)=[(k-1)*dt+j*dt/20,local'];
  violation=max([violation;box(1:3)'-local(1:3);local(1:3)-box(4:6)']);
  if j==20,break;end
  u0=U(:,k)+j/20*(U(:,k+1)-U(:,k));u1=U(:,k)+(j+1)/20*(U(:,k+1)-U(:,k));
  local=rk4_step(local,u0,u1,dt/20);forward=rk4_step(forward,u0,u1,dt/20);
  rollout((k-1)*20+j+2,:)=[(k-1)*dt+(j+1)*dt/20,forward'];
 end
 defect=max(defect,norm(local(1:3)-X(1:3,k+1)));
end
clearance=inf;globalClear=inf;boxClear=inf;
for b=1:size(obs,1)
 clearance=min(clearance,min(vecnorm(dense(:,2:3)-obs(b,1:2),2,2)-obs(b,3)));
 globalClear=min(globalClear,min(vecnorm(rollout(:,2:3)-obs(b,1:2),2,2)-obs(b,3)));
end
for b=1:size(boxes,1),boxClear=min(boxClear,box_clearance(boxes(b,1:3),boxes(b,4:6),obs));end
checks=struct('max_scaled_constraint_violation',residual,'min_dense_clearance_m',clearance,...
 'min_whole_box_clearance_m',boxClear,'max_dense_corridor_violation_m',violation,...
 'max_interval_position_defect_m',defect,'spline_certified_clearance_m',splineClearance,...
 'open_loop_terminal_position_error_m',norm(forward(1:3)-X(1:3,end)),'open_loop_min_clearance_m',globalClear);
assert(residual<=1e-5 && violation<=1e-5 && clearance>=cfg.safety_margin_m && boxClear>cfg.safety_margin_m && defect<=.5,...
 'Post-solve validation failed: refine the mesh or initialization.');
result=struct('name',c.name,'X',X','U',U','initial_X',initialX','initial_U',initialU',...
 'astar_path',raw,'spline_path',splinePath,'spline_knots',knots,'boxes',boxes,...
 'interval_box',owner-1,'columns',obs,'time',(0:N)'*dt,'dense',dense,'open_loop_rollout',rollout,...
 'objective',J,'iterations',stats.iterations,'seconds',seconds,'solver_exitflag',flag,...
 'solver_message',stats.message,'astar_expanded_nodes',expanded,'history',history,...
 'checks',checks,'derivative_check_max_error',derr);

 function [xx,uu]=unpack(z)
  w=z.*scale;xx=reshape(w(1:7*M),7,M);uu=reshape(w(7*M+1:end),3,M);
 end
 function [J,dJ]=objective(z)
  [xx,uu]=unpack(z);du=uu-[0;1;0];w=[.15;2;.8];
  J=dt*sum(sum(w.*(du./su).^2.*W))+.02*dt*sum((xx(2,:)/100).^2.*W);
  gx=zeros(size(xx));gx(2,:)=.04*dt*xx(2,:)/1e4.*W;
  gu=2*dt*w.*du./su.^2.*W;dJ=[gx(:);gu(:)].*scale;
 end
 function [ine,eq,GI,GE]=constraints(z)
  [xx,uu]=unpack(z);F=zeros(7,M);AA=zeros(7,7,M);BB=zeros(7,3,M);
  for kk=1:M,[F(:,kk),AA(:,:,kk),BB(:,:,kk)]=fw_motion(xx(:,kk),uu(:,kk));end
  eq=reshape((xx(:,2:end)-xx(:,1:end-1)-dt/2*(F(:,1:end-1)+F(:,2:end)))./sx,[],1);
  ine=zeros(size(obs,1)*M,1);GI=zeros(numel(z),numel(ine));GE=zeros(numel(z),numel(eq));
  for bb=1:size(obs,1)
   P=xx(1:2,:)-obs(bb,1:2)';rr=obs(bb,3)+cfg.safety_margin_m;
   ine((bb-1)*M+(1:M))=1-sum(P.^2,1)'/rr^2;
   for kk=1:M,GI((kk-1)*7+(1:2),(bb-1)*M+kk)=-2*P(:,kk).*sx(1:2)/rr^2;end
  end
  for kk=1:N
   ix=(kk-1)*7+(1:7);iu=7*M+(kk-1)*3+(1:3);col=(kk-1)*7+(1:7);
   GE(ix,col)=(diag(1./sx)*(-eye(7)-dt/2*AA(:,:,kk))*diag(sx))';
   GE(ix+7,col)=(diag(1./sx)*(eye(7)-dt/2*AA(:,:,kk+1))*diag(sx))';
   GE(iu,col)=(diag(1./sx)*(-dt/2*BB(:,:,kk))*diag(su))';
   GE(iu+3,col)=(diag(1./sx)*(-dt/2*BB(:,:,kk+1))*diag(su))';
  end
 end
 function stop=record(~,v,state)
  stop=false;if strcmp(state,'iter'),history=[history;v.fval,v.constrviolation];end
 end
end
