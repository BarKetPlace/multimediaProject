function plotMFCC(ifig, Ex,Ey,Exhat)
% function plotMFCC(Ex,Ey,Exhat)
[M, nbframe]=size(Ex);

figure(ifig), clf;
stem(M*[1:nbframe], max(Ex(:))*ones(1,nbframe),'--k'); hold on;
plot(Ey(:),'LineWidth',2); hold on;
plot(Ex(:),'LineWidth',2); hold on;
plot(Exhat(:));hold on;
legend('Frames','Ey','Ex','Exhat');
xlabel('Mel space coefficent');
ylabel('Value f coefficient');
set(gca,'Xtick',[round(M/2):M:length(Ex(:))]);
set(gca,'XtickLabel',[1:nbframe]);
title('Denoising results');
end

