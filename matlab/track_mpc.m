function output=track_mpc(S,cfg)
% Reference-linearized tracking MPC, input bounds only, exact-state feedback.
% Full nonlinear plant, initial displacement/heading error, no wind/estimator.
h=.2;M=15;tt=(0:round(cfg.duration_s/h))*h;
ref=interp1(S.time,S.X,tt,'pchip')';uu=interp1(S.time,S.U,tt,'linear')';
last=ref(:,end);
for j=1:M
 ref(:,end+1)=last+[20*j*h;0;0;0;0;0;0];uu(:,end+1)=[0;1;0];
end
actual=ref(:,1)+[0;5;3;0;0;pi/60;0];opened=actual;
rows=zeros(numel(tt),22);Uhist=zeros(numel(tt)-1,3);residuals=zeros(numel(tt)-1,1);guess=zeros(3*M,1);
AA=zeros(7,7,size(ref,2)-1);BB=zeros(7,3,size(ref,2)-1);QQ=zeros(7,7,size(ref,2)-1);defects=zeros(7,size(ref,2)-1);
for i=1:size(ref,2)-1
 [~,Ac,Bc]=fw_motion(ref(:,i),uu(:,i));aug=expm([Ac Bc;zeros(3,10)]*h);
 AA(:,:,i)=aug(1:7,1:7);BB(:,:,i)=aug(1:7,8:10);
 [~,Av,~]=fw_motion(ref(:,i+1),uu(:,i));C=[eye(3),zeros(3,4);Av(1:3,:)];
 QQ(:,:,i)=C'*diag([ones(1,3)/25 ones(1,3)/4])*C+diag([0 0 0 .01 .03 .03 .03]);
 defects(:,i)=ref(:,i+1)-rk4_step(ref(:,i),uu(:,i),uu(:,i),h);
end
for k=1:numel(tt)-1
 H=zeros(3*M);f=zeros(3*M,1);P=eye(7);G=zeros(7,3*M);offset=zeros(7,1);e=ref(:,k)-actual;
 for j=1:M
  idx=k+j-1;A=AA(:,:,idx);B=BB(:,:,idx);Q=QQ(:,:,idx);
  P=A*P;G=A*G;G(:,(j-1)*3+(1:3))=G(:,(j-1)*3+(1:3))-B;
  offset=A*offset+defects(:,idx);H=H+2*G'*Q*G;f=f+2*G'*Q*(P*e+offset);
 end
 H=H+2*kron(eye(M),diag([.02 .15 .05]))+1e-9*eye(3*M);H=(H+H')/2;
 lo=reshape(cfg.input_lower(:)-uu(:,k:k+M-1),[],1);hi=reshape(cfg.input_upper(:)-uu(:,k:k+M-1),[],1);
 [du,res]=box_qp(H,f,lo,hi,guess);residuals(k)=res;
 guess=[du(4:end);zeros(3,1)]; % Shift the previous solution: warm start.
 u=uu(:,k)+du(1:3);Uhist(k,:)=u';rows(k,:)=[tt(k),actual',opened',ref(:,k)'];
 for q=1:5,actual=rk4_step(actual,u,u,h/5);opened=rk4_step(opened,uu(:,k),uu(:,k),h/5);end
end
rows(end,:)=[tt(end),actual',opened',ref(:,numel(tt))'];
eM=vecnorm(rows(:,2:4)-rows(:,16:18),2,2);eO=vecnorm(rows(:,9:11)-rows(:,16:18),2,2);clearance=inf;
for b=1:size(S.columns,1),o=S.columns(b,:);clearance=min(clearance,min(vecnorm(rows(:,2:3)-o(1:2),2,2)-o(3)));end
output=struct('time',tt','rows',rows,'inputs',Uhist,'qp_residuals',residuals,'sample_time_s',h,...
 'prediction_steps',M,'mpc_position_error',eM,'open_loop_position_error',eO,...
 'min_sampled_mpc_clearance_m',clearance,'rms_mpc_position_error_m',sqrt(mean(eM.^2)),...
 'rms_open_loop_position_error_m',sqrt(mean(eO.^2)));
end

function [x,residual]=box_qp(H,f,lo,hi,guess)
% Strictly convex box QP: log-barrier damped Newton iterations.
x=min(max(guess,lo+.001*(hi-lo)),hi-.001*(hi-lo));
for mu=[1 .1 .01 .001 1e-4 1e-5 1e-6 1e-7]
 for it=1:70
  dl=x-lo;du=hi-x;grad=H*x+f-mu./dl+mu./du;Hess=H+diag(mu./dl.^2+mu./du.^2);p=-Hess\grad;
  if -grad'*p<1e-12,break;end
  alpha=1;neg=p<0;pos=p>0;
  if any(neg),alpha=min(alpha,.99*min(-dl(neg)./p(neg)));end
  if any(pos),alpha=min(alpha,.99*min(du(pos)./p(pos)));end
  old=.5*x'*H*x+f'*x-mu*sum(log(dl)+log(du));
  while alpha>1e-12
   xx=x+alpha*p;val=.5*xx'*H*xx+f'*xx-mu*sum(log(xx-lo)+log(hi-xx));
   if val<=old+1e-4*alpha*grad'*p,break;end
   alpha=alpha/2;
  end
  x=x+alpha*p;
 end
end
residual=norm(x-min(max(x-(H*x+f),lo),hi),inf);assert(residual<2e-3,'MPC QP projected-gradient check failed.');
end
