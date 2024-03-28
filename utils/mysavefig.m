function mysavefig(h, filename, outdir, fontsize, aspect, ticklength)
if ~exist('ticklength','var')
    ticklength = 0;
end
set(gca,'FontSize',fontsize);
set(gca,'FontName','Arial')
set(gca,'TickDir','in');
ax = gca;
ax.TickLength = ax.TickLength*ticklength;
%set(gca,'LineWidth',1); 
xl = get(gca,'XLabel');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', fontsize)
set(xl, 'FontSize', fontsize);
yl = get(gca,'YLabel');
yAX = get(gca,'YAxis');
set(yAX,'FontSize', fontsize)
set(yl, 'FontSize', fontsize);
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 aspect];
saveas(h,fullfile(outdir,sprintf('%s.pdf',filename)),'pdf');