function [boxes,owner,nodeLo,nodeHi]=build_corridors(seed,obs,cfg)
% Two intervals per box. Greedily expand six faces while the whole box is safe.
N=cfg.intervals;assert(mod(N,2)==0,'intervals must be even.');
lower=cfg.state_lower(1:3)';upper=cfg.state_upper(1:3)';margin=cfg.safety_margin_m;
nodeLo=repmat(lower,N+1,1);nodeHi=repmat(upper,N+1,1);boxes=zeros(N/2,6);owner=zeros(N,1);
for ib=1:N/2
 ids=(2*ib-1):(2*ib+1);owner(ids(1:end-1))=ib;
 lo=max(min(seed(ids,:),[],1)-[2 2 8],lower);hi=min(max(seed(ids,:),[],1)+[2 2 8],upper);
 assert(box_clearance(lo,hi,obs)>margin,'Seed box collides: refine intervals or increase seed clearance.');
 baseLo=lo;baseHi=hi;cap=[34 40 18];
 for grow=1:26
  for face=1:6
   d=ceil(face/2);cl=lo;ch=hi;if mod(face,2),cl(d)=lo(d)-2;else,ch(d)=hi(d)+2;end
   if cl(d)<baseLo(d)-cap(d)||ch(d)>baseHi(d)+cap(d),continue;end
   cl=max(cl,lower);ch=min(ch,upper);
   if box_clearance(cl,ch,obs)>margin+.05,lo=cl;hi=ch;end
  end
 end
 boxes(ib,:)=[lo hi];nodeLo(ids,:)=max(nodeLo(ids,:),lo+1);nodeHi(ids,:)=min(nodeHi(ids,:),hi-1);
end
assert(all(nodeLo(:)<=nodeHi(:)),'Consecutive corridors have empty node intersections.');
end
