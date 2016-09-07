function scatter_bougnoux()
global fgnb
n=100;
fgnb=getFgnb(n);
colors=repmat([1 0 0],n,1);
fgnb1=num2cell(fgnb,2);
choice=cellfun(@isreal,fgnb1);
chosen=colors(choice,:);
colors(choice,:)=repmat([0 0 1],size(chosen,1),1);
scatter([abs(fgnb(:,1)); 900],[abs(fgnb(:,2)); 1100],20,[colors; [0 1 0]],'+')
hold on
plot([0 7200],[0 8800],'y')
title('Scatter plot of bougnoux results. b=real, r=imag, g=without noise');
xlabel('abs(f2)');
ylabel('abs(f1)');
hold off
end

function a=realify(x)

if isreal(x)
    a = x
end
    a = abs(x)
end

function fgnb=getFgnb(n)
fgnb=zeros(n,2);
tic();
points=rand(3,60);
for i=1:n
    temp=FDemo(900,1100,points);
    fgnb(i,:)=temp;
end
toc();
end