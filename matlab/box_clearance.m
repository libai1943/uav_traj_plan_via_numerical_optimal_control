function d=box_clearance(lo,hi,obs)
% Exact whole-box clearance, since columns span all permitted altitudes.
dxy=max(max(lo(1:2)-obs(:,1:2),obs(:,1:2)-hi(1:2)),0);
d=min(sqrt(sum(dxy.^2,2))-obs(:,3));
end
