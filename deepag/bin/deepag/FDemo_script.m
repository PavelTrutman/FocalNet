n=10
fgnb=zeros(10,2)
for i=1:n
    temp=FDemo();
    fgnb(i,:)=temp;
end
%plot(fgnb)
fgnb(:,1)./fgnb(:,2)

