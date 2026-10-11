function d=segment_clearance(a,b,obs)
% Minimum distance from the whole segment to physical vertical columns.
v=b(1:2)-a(1:2);w=obs(:,1:2)-a(1:2);t=max(0,min(1,(w*v')/max(v*v',1e-30)));
d=min(sqrt(sum((w-t.*v).^2,2))-obs(:,3));
end
