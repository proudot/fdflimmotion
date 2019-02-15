sampleNb=50;
E=load('build/energy.dlm');E1=reshape(E,[sampleNb,sampleNb,sampleNb]);
tickNb=5;
phaseTicks=linspace(0,2*pi,tickNb)
ampTicks=linspace(0,60000,tickNb)
CTicks=linspace(1000,19000,tickNb)

init=[10408.9  1.14541 2.18401];

tickPos=linspace(1,sampleNb,tickNb)

figure();
subplot(1,2,1);
imshow((mat2gray(log(squeeze([min(E1,[],3)]))))); 
axis on
hold on
plot(init(1)*init(2)*sampleNb/max(ampTicks),init(3)*sampleNb/max(phaseTicks),'g.', 'MarkerSize',20);
ampPhaseE=mat2gray((squeeze([min(E1,[],3)])));

lm=locmax2d(-ampPhaseE,1);lm(lm~=0)=255;
[x,y]=find(255-lm);
plot(y,x,'r.','MarkerSize',20);

set(gca,'XTick',tickPos)
set(gca,'YTick',tickPos)
set(gca,'XTickLabel',arrayfun(@(x) num2str(x),ampTicks,'unif',0));
set(gca,'YTickLabel',arrayfun(@(x) num2str(x),phaseTicks,'unif',0));

xlabel('Amplitude (a.i.u.)')
ylabel('Phase (rad)')


subplot(1,2,2);
imshow((mat2gray(log(squeeze([min(E1,[],2)])))));
axis on
hold on
CPhaseE=mat2gray((squeeze([min(E1,[],2)])));
lm=locmax2d(-CPhaseE,1);lm(lm~=0)=255;
[x,y]=find(255-lm);
plot(y,x,'r.', 'MarkerSize',20);
plot(init(1)*sampleNb/max(CTicks),init(3)*sampleNb/max(phaseTicks),'g.','MarkerSize',20);
set(gca,'XTick',tickPos)
set(gca,'YTick',tickPos)
set(gca,'XTickLabel',arrayfun(@(x) num2str(x),CTicks,'unif',0));
set(gca,'YTickLabel',arrayfun(@(x) num2str(x),phaseTicks,'unif',0));

xlabel('Constant component (a.i.u.)')
ylabel('Phase (rad)')
