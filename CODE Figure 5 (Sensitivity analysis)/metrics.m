    function metric_vector = metrics(time,sol,timeORIG,solORIG)

locs = find(time>1);
locsORIG = find(timeORIG>1);
metric_vector(1) = max(sol(1,locs))-max(solORIG(1,locsORIG)); %max V
metric_vector(2) = max(sol(5,:))-max(solORIG(5,:)); %D max
metric_vector(3) = min(sol(2,:)+sol(4,:))-min(solORIG(2,:)+solORIG(4,:)); %min S+R
%metric_vector(4) = min(sol(2,:))-min(solORIG(2,:));% min S
metric_vector(4) = max(sol(7,:))-max(solORIG(7,:)); %max MphiI
metric_vector(5) = max(sol(10,:))-max(solORIG(10,:)); %max T
metric_vector(6) = max(sol(11,:))-max(solORIG(11,:)); %max IL-6
metric_vector(7) = max(sol(17,:))-max(solORIG(17,:)); %max IFN


for i = 1:length(sol(17,:))-1
   tgrid_vec(i) = time(i+1)-time(i); 
end
metric_vector(8) = sum(sol(17,1:end-1).*tgrid_vec);


tunder =  find([sol(2,:)+sol(4,:)]<=0.16*0.3);
tabove =  find(sol(2,tunder:end)+sol(4,tunder:end)>=0.16*0.3);
if isempty(tabove)==1 & isempty(tunder)==1
    metric_vector(9) = 0-1.758920349788992;
else
    if isempty(tabove)==1
        tabove = 30;
    end
    metric_vector(9) = time(tabove(1)+tunder(1))-time(tunder(1))-1.758920349788992;
end

peak_loc = find(sol(17,:)==max(sol(17,:)));
metric_vector(10) = time(peak_loc)-1.797097448511531; %time of IFN peak


end