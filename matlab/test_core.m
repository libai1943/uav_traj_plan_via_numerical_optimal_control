function test_core
% Geometry failures that a node-only collision check would miss; model checks.
obs=[0 0 1];assert(segment_clearance([-2 0 100],[-2 0 100],obs)>0);
assert(segment_clearance([2 0 100],[2 0 100],obs)>0);
assert(segment_clearance([-2 0 100],[2 0 100],obs)<0);
assert(abs(box_clearance([-2 -2 80],[2 2 180],obs)+1)<1e-12);
assert(abs(box_clearance([2 -1 80],[3 1 180],obs)-1)<1e-12);
assert(norm(fw_motion([0;0;100;20;0;0;0],[0;1;0])-[20;0;0;0;0;0;0])<1e-12);
x=[12;-8;115;19;.12;-.18;.25];u=[.2;1.15;-.07];[~,A,B]=fw_motion(x,u);
for j=1:7,dx=zeros(7,1);dx(j)=1e-6;assert(norm((fw_motion(x+dx,u)-fw_motion(x-dx,u))/2e-6-A(:,j),inf)<2e-8);end
for j=1:3,du=zeros(3,1);du(j)=1e-6;assert(norm((fw_motion(x,u+du)-fw_motion(x,u-du))/2e-6-B(:,j),inf)<2e-8);end
disp('Geometry and dynamics regression tests passed.');
end
