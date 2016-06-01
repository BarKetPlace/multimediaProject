function plotAdditivity( ifig, Ex,Ey,En_model )
% function plotAdditivity( ifig, Ex,Ey,En_model )


[M, nbframe]=size(Ex);
figure(ifig), clf; 

stem(M*[1:nbframe],max(Ey(:))*ones(1,nbframe),'--k'); hold on;
plot(Ey(:));hold on;
plot(Ex(:)+En_model(:));
legend('Frames','Ey', 'Ex + En_{model}');
set(gca,'Xtick',[round(M/2):M:length(Ex(:))]);
set(gca,'XtickLabel',[1:nbframe]);
xlabel('Frame number');
title('Additivity check');

end

