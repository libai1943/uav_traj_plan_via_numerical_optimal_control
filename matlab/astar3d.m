function [path,expanded] = astar3d(start,goal,obs,cfg)
% ASTAR3D 26-neighbor XYZ search; exact XY segment-column collision checking.
step=cfg.grid_step_m(:)';lo=cfg.search_lower_m(:)';hi=cfg.search_upper_m(:)';
dims=round((hi-lo)./step)+1;s=round((start-lo)./step)+1;g=round((goal-lo)./step)+1;
assert(all(abs(lo+(s-1).*step-start)<1e-9) && all(abs(lo+(g-1).*step-goal)<1e-9),'Endpoints must lie on the 3-D grid.');
margin=cfg.safety_margin_m+cfg.astar_extra_clearance_m;
assert(segment_clearance(start,start,obs)>margin && segment_clearance(goal,goal,obs)>margin,'A* endpoint clearance failed.');
[ix,iy,iz]=ndgrid(1:dims(1),1:dims(2),1:dims(3));coords=[ix(:),iy(:),iz(:)];points=lo+(coords-1).*step;
n=size(coords,1);cost=inf(n,1);f=inf(n,1);order=inf(n,1);parent=zeros(n,1);closed=false(n,1);
si=sub2ind(dims,s(1),s(2),s(3));gi=sub2ind(dims,g(1),g(2),g(3));
cost(si)=0;f(si)=norm(goal-start);order(si)=0;serial=0;expanded=0;
offsets=[];for a=-1:1,for b=-1:1,for c=-1:1,if any([a b c]),offsets(end+1,:)=[a b c];end,end,end,end %#ok<AGROW>
while true
 best=min(f);assert(isfinite(best),'3-D A* found no collision-free grid path.');
 ids=find(f==best);[~,j]=min(order(ids));id=ids(j);
 if id==gi
  ids=gi;while ids(1)~=si,ids=[parent(ids(1));ids];end %#ok<AGROW>
  path=points(ids,:);return;
 end
 f(id)=inf;closed(id)=true;expanded=expanded+1;a=points(id,:);
 for k=1:26
  next=coords(id,:)+offsets(k,:);if any(next<1)||any(next>dims),continue;end
  ni=sub2ind(dims,next(1),next(2),next(3));if closed(ni),continue;end
  b=points(ni,:);if segment_clearance(a,b,obs)<=margin,continue;end
  alt=cost(id)+norm(b-a);
  if alt<cost(ni)-1e-10
   cost(ni)=alt;parent(ni)=id;serial=serial+1;order(ni)=serial;f(ni)=alt+norm(goal-b);
  end
 end
end
end
