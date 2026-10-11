function [seed,dense,knots,clearance]=spline_seed(raw,obs,cfg)
% Visibility shortcut -> clamped cubic spline -> arc-length resampling.
% The cubic second derivative bounds the deviation from each checked chord.
margin=cfg.safety_margin_m;knots=raw(1,:);i=1;
while i<size(raw,1)
 j=size(raw,1);
 while j>i+1 && segment_clearance(raw(i,:),raw(j,:),obs)<=margin+cfg.astar_extra_clearance_m,j=j-1;end
 knots(end+1,:)=raw(j,:);i=j; %#ok<AGROW>
end
alt=raw(1:2:end,:);if any(alt(end,:)~=raw(end,:)),alt(end+1,:)=raw(end,:);end
candidates={knots,alt,raw};
for ci=1:3
 knots=candidates{ci};arc=[0;cumsum(sqrt(sum(diff(knots).^2,2)))];
 pp=spline(arc',[[1;0;0],knots',[1;0;0]]);
 [breaks,coefs,~,~,dim]=unmkpp(pp);pp2=mkpp(breaks,[6*coefs(:,1),2*coefs(:,2)],dim);
 q=unique([linspace(0,arc(end),1601),arc']);dense=ppval(pp,q)';
 s2=ppval(pp2,q)';curvature=max(vecnorm(s2(1:end-1,:),2,2),vecnorm(s2(2:end,:),2,2));
 errs=curvature.*diff(q').^2/8;clearance=inf;
 for k=1:numel(q)-1,clearance=min(clearance,segment_clearance(dense(k,:),dense(k+1,:),obs)-errs(k));end
 lo=cfg.state_lower(1:3)';hi=cfg.state_upper(1:3)';
 if clearance>margin+2 && all(all(dense>=lo+max(errs))) && all(all(dense<=hi-max(errs)))
  distance=[0;cumsum(vecnorm(diff(dense),2,2))];
  query=interp1(distance,q',linspace(0,distance(end),cfg.intervals+1));seed=ppval(pp,query)';return;
 end
end
error('Cubic spline is unsafe. Increase A* clearance or refine the grid.');
end
