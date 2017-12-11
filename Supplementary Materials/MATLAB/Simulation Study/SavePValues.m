function []=SavePValues(test,dist,n,n1,n2,r,res,Kg,PValuesR,PValuesC)  %#ok<INUSL>
PValuesCorrected=sort(PValuesC);
PValuesRaw=sort(PValuesR);
dir='PValues\';

if ~exist(dir,'dir'),  mkdir(dir); end;
if ~exist([dir test],'dir'),  mkdir([dir test]); end;
if ~exist([dir test '\p'],'dir'),    mkdir([dir test '\p']); end;
if ~exist([dir test '\eps'],'dir'),  mkdir([dir test '\eps']); end;
if ~exist([dir test '\fig'],'dir'),  mkdir([dir test '\fig']); end;
if ~exist([dir test '\png'],'dir'),  mkdir([dir test '\png']); end;
if ~exist([dir test '\png\Small'],'dir'),   mkdir([dir test '\png\Small']); end;
if ~exist([dir test '\png\Normal'],'dir'),  mkdir([dir test '\png\Normal']); end;
if ~exist([dir test '\png\Large'],'dir'),   mkdir([dir test '\png\Large']); end;


if strcmp(test,'OST')
    fbasename=[test num2str(dist) 'N' num2str(n) 'R' num2str(r)];
else
    fbasename=[test num2str(dist) 'N1' num2str(n1) 'N2' num2str(n2) 'R' num2str(r)];
end
fname_raw      =[fbasename 'R.txt'];
fname_corrected=[fbasename 'C.txt'];

save([dir test '\p\' fname_raw],'PValuesRaw','-ascii','-double');
save([dir test '\p\' fname_corrected],'PValuesCorrected','-ascii','-double');
save([dir test '\Kg' strtok(fbasename,'R') '.txt'],'Kg','-ascii','-double');

N=r*res;
CDF_PValuesRaw=(1:length(PValuesRaw))/N;
CDF_PValuesCorrected=(1:length(PValuesCorrected))/N;

fig=figure('Visible','Off','name',['Correction of p-values - ' test],...
       'units',' centimeters',...
       'position',[5   5   15   15],...
       'Color','white');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [15 15]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 15 15]);   

plot([0 1], [0 1],'b');
hold on;
axis square;
set(gca,'PlotBoxAspectRatio',[1 1 1]);
set(gca,'XLim',[0 1]);
set(gca,'YLim',[0 1]);
set(gca,'XTick',0:0.2:1);
set(gca,'YTick',0:0.2:1);
set(gca,'units','centimeters');
set(gca,'FontWeight','light');
set(gca,'FontName','Courier New');
set(gca,'Position',[0.8 0.8 13.4 13.4]);

plot(r*PValuesRaw,r*CDF_PValuesRaw,'k');
plot(r*PValuesCorrected,r*CDF_PValuesCorrected,'r');


print(gcf,'-r0','-dpng','-painters',[dir test '\png\Normal\' fbasename 'N' '.png']);
saveas(gcf,[dir test '\fig\' fbasename '.fig']);

% Save .eps
set(gca,'FontSize',25);
set(gca,'Position',[1.75 1 12 12]);
print(fig,'-depsc2',[dir test '\eps\' fbasename '.eps']);
set(gca,'Position',[1.5 1.5 12 12]);
set(gca,'FontSize',20);
print(gcf,'-r50','-dpng', '-painters',[dir test '\png\Small\'  fbasename '.png']);
set(gca,'FontSize',3);
print(gcf,'-r300','-dpng', '-painters',[dir test '\png\Large\'  fbasename '.png']);
close(fig);