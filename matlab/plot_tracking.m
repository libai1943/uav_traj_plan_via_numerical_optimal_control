function plot_tracking(r,folder,showFigures)
visible='off';if showFigures,visible='on';end
fig=figure('Color','w','Visible',visible,'Position',[70 70 1200 460]);blue=[0 .4 .68];orange=[.67 .41 .12];
ax=subplot(1,2,1,'Parent',fig);hold(ax,'on');
plot(ax,r.rows(:,16),r.rows(:,17),'--','Color',[.45 .45 .45],'LineWidth',1.5);
plot(ax,r.rows(:,9),r.rows(:,10),'Color',orange,'LineWidth',1.5);
plot(ax,r.rows(:,2),r.rows(:,3),'Color',blue,'LineWidth',1.8);axis(ax,'equal');
xlabel(ax,'X (m)');ylabel(ax,'Y (m)');title(ax,'Same initial position and heading error');
legend(ax,{'Planned reference','Open loop','MPC feedback'},'Location','best','Box','off');grid(ax,'on');
ax=subplot(1,2,2,'Parent',fig);hold(ax,'on');
plot(ax,r.time,r.open_loop_position_error,'Color',orange,'LineWidth',1.5);
plot(ax,r.time,r.mpc_position_error,'Color',blue,'LineWidth',1.8);
xlabel(ax,'Time (s)');ylabel(ax,'Position error norm (m)');title(ax,'Position / velocity tracking MPC');
legend(ax,{'Open loop','MPC feedback'},'Location','best','Box','off');grid(ax,'on');
set(fig,'PaperUnits','inches','PaperPosition',[0 0 12 4.6],'PaperSize',[12 4.6],'InvertHardcopy','off');
print(fig,fullfile(folder,'mpc_tracking.png'),'-dpng','-r240');if ~showFigures,close(fig);end
end
