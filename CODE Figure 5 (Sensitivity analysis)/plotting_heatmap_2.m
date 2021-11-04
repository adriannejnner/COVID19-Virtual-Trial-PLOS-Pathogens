function plotting_heatmap(mega_metric_matrix,amount1,amount2,j,col,comp_mat)

for i = amount1:amount2
   grid_metric_1 = linspace(min(comp_mat(:,i)),max(comp_mat(:,i)),11);
   if isnan(mega_metric_matrix(j,i))==1
       fill([i-1 i-1 i i],[j-1 j j j-1],'k','EdgeColor','none') 
   else
       val = find(grid_metric_1>=mega_metric_matrix(j,i),1);
       fill([i-1 i-1 i i],[j-1 j j j-1],col(val,:),'EdgeColor','none') 
   end
end

end