function res = plot_error_ratio(Image,Title)

[MaxDur,MaxShift]=size(Image);
Legend = cell(MaxDur,1);
l = 1;
for RD = 1:10:MaxDur
         hold on
         plot(1:2:25,Image(RD,1:2:MaxShift),'-o','LineWidth',2)
         set(gca,'fontsize',18)

         Legend{l} = strcat('RD=',num2str(RD+1));
         l=l+1;

         ax=gca;
         title(Title);
         xlabel('Shift');
         ylabel('RED');

         ax.XTick=0:3:MaxShift;
         ax.YDir = 'normal';         
end
hold on
plot(0:26,zeros(27,1),'k--','LineWidth',1.5)
h = legend(Legend(1:l-1));
set(h,'FontSize',10)

res = 1;
end

