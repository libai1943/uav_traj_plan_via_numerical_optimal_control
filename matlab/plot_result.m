function plot_result(r,cfg,folder,showFigures)
% Equal-metric 3-D axes, translucent columns, and all corridor fragments.
visible='off';if showFigures,visible='on';end
blue=[0 .4 .68];orange=[.67 .41 .12];teal=[.28 .56 .60];
X=r.X;U=r.U;raw=r.astar_path;seed=r.spline_path;
fig=figure('Color','w','Visible',visible,'Position',[30 30 1300 720]);ax=axes(fig);hold(ax,'on');
set(ax,'FontName','Arial','FontSize',12,'Projection','orthographic','GridAlpha',.15);
[cx,cy,cz]=cylinder(1,80);
for i=1:size(r.columns,1)
 o=r.columns(i,:);surf(ax,o(1)+o(3)*cx,o(2)+o(3)*cy,200*cz,'FaceColor',[.5 .5 .5],'FaceAlpha',.22,'EdgeColor','none');
 for h=[0 200],plot3(ax,o(1)+o(3)*cx(1,:),o(2)+o(3)*cy(1,:),h+0*cx(1,:),'Color',[.55 .55 .55],'LineWidth',.7);end
end
for i=1:size(r.boxes,1)
 b=r.boxes(i,:);v=zeros(8,3);for q=0:7,for d=1:3,if bitget(q,d),v(q+1,d)=b(d+3);else,v(q+1,d)=b(d);end,end,end
 faces=[1 2 4 3;5 6 8 7;1 2 6 5;3 4 8 7;1 3 7 5;2 4 8 6];
 patch(ax,'Vertices',v,'Faces',faces,'FaceColor',teal,'FaceAlpha',.035,'EdgeColor',teal,'EdgeAlpha',.45,'LineWidth',.55);
end
h1=plot3(ax,raw(:,1),raw(:,2),raw(:,3),'Color',[.55 .55 .55],'LineWidth',1.2);
h2=plot3(ax,seed(:,1),seed(:,2),seed(:,3),'--','Color',orange,'LineWidth',1.8);
h3=plot3(ax,X(:,1),X(:,2),X(:,3),'Color',blue,'LineWidth',3);
plot3(ax,X(1,1),X(1,2),X(1,3),'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue);
plot3(ax,X(end,1),X(end,2),X(end,3),'p','MarkerFaceColor',blue,'MarkerEdgeColor',blue,'MarkerSize',12);
h4=patch(ax,nan,nan,teal,'FaceAlpha',.2,'EdgeColor',teal);h5=patch(ax,nan,nan,[.8 .8 .8],'EdgeColor',[.5 .5 .5]);
legend(ax,[h1 h2 h3 h4 h5],{'3-D A* path','Cubic-spline initial guess','Optimized trajectory','Safe corridor fragments','Column obstacles'},...
 'Location','southoutside','NumColumns',3,'Box','off','FontSize',12);
xlabel(ax,'X (m)');ylabel(ax,'Y (m)');zlabel(ax,'Z (m)');axis(ax,[-20 640 -220 220 0 205]);daspect(ax,[1 1 1]);view(ax,-20,23);grid(ax,'on');
title(ax,[strrep(r.name,'_',' '),' | fixed-wing trajectory and safe corridors']);export(fig,'trajectory_3d',13.5,7.2);

fig=figure('Color','w','Visible',visible,'Position',[40 40 1200 660]);ax=axes(fig);hold(ax,'on');
for i=1:size(r.boxes,1)
 b=r.boxes(i,:);patch(ax,[b(1) b(4) b(4) b(1)],[b(2) b(2) b(5) b(5)],teal,'FaceAlpha',.045,'EdgeColor',teal,'EdgeAlpha',.4,'LineWidth',.55);
end
th=linspace(0,2*pi,160);
for i=1:size(r.columns,1)
 o=r.columns(i,:);patch(ax,o(1)+o(3)*cos(th),o(2)+o(3)*sin(th),[.82 .82 .82],'EdgeColor',[.4 .4 .4]);
 rr=o(3)+cfg.safety_margin_m;plot(ax,o(1)+rr*cos(th),o(2)+rr*sin(th),':','Color',[.4 .4 .4]);
end
h1=plot(ax,raw(:,1),raw(:,2),'Color',[.55 .55 .55],'LineWidth',1.3);
h2=plot(ax,seed(:,1),seed(:,2),'--','Color',orange,'LineWidth',1.7);h3=plot(ax,X(:,1),X(:,2),'Color',blue,'LineWidth',2.7);
plot(ax,X(1,1),X(1,2),'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue);
plot(ax,X(end,1),X(end,2),'p','MarkerFaceColor',blue,'MarkerEdgeColor',blue,'MarkerSize',12);
legend(ax,[h1 h2 h3],{'3-D A* path','Cubic-spline initial guess','Optimized trajectory'},'Location','southoutside','Orientation','horizontal','Box','off');
axis(ax,[-20 640 -220 220]);axis(ax,'equal');xlabel(ax,'X (m)');ylabel(ax,'Y (m)');grid(ax,'on');
set(ax,'GridAlpha',.15);title(ax,'Plan view | dotted circles include the 8 m safety margin');export(fig,'trajectory_top',12,6.6);

fig=figure('Color','w','Visible',visible,'Position',[50 50 1300 700]);
data=[X(:,4),rad2deg(X(:,5)),rad2deg(X(:,7)),U(:,1),U(:,2),rad2deg(U(:,3))];
labels={'Airspeed V (m/s)','Flight-path angle gamma (deg)','Bank angle mu (deg)',...
 'Tangential input a_t (m/s^2)','Load factor n','Bank rate omega_mu (deg/s)'};
bounds=[cfg.state_lower(4),cfg.state_upper(4);rad2deg([cfg.state_lower(5),cfg.state_upper(5)]);...
 rad2deg([cfg.state_lower(7),cfg.state_upper(7)]);cfg.input_lower(1),cfg.input_upper(1);...
 cfg.input_lower(2),cfg.input_upper(2);rad2deg([cfg.input_lower(3),cfg.input_upper(3)])];
for j=1:6
 ax=subplot(2,3,j,'Parent',fig);plot(ax,r.time,data(:,j),'Color',blue,'LineWidth',1.8);hold(ax,'on');
 yline(ax,bounds(j,1),':','Color',[.5 .5 .5]);yline(ax,bounds(j,2),':','Color',[.5 .5 .5]);
 xlabel(ax,'Time (s)');ylabel(ax,labels{j});grid(ax,'on');set(ax,'GridAlpha',.15);
end
sgtitle(fig,'Optimized states and inputs | dotted lines: prescribed bounds');export(fig,'states_inputs',13,7);

fig=figure('Color','w','Visible',visible,'Position',[60 60 1100 400]);
ax=subplot(1,2,1,'Parent',fig);plot(ax,1:size(r.history,1),r.history(:,1),'Color',blue);xlabel(ax,'SQP iteration');ylabel(ax,'Objective J');grid(ax,'on');
ax=subplot(1,2,2,'Parent',fig);semilogy(ax,1:size(r.history,1),max(r.history(:,2),1e-16),'Color',blue);
xlabel(ax,'SQP iteration');ylabel(ax,'Maximum scaled constraint violation');grid(ax,'on');export(fig,'convergence',11,4);

function export(f,name,w,h)
 set(f,'PaperUnits','inches','PaperPosition',[0 0 w h],'PaperSize',[w h],'InvertHardcopy','off');
 print(f,fullfile(folder,[name '.png']),'-dpng','-r240');if ~showFigures,close(f);end
end
end
