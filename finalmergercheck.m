function [histodiff,intensitydiff,textdiff]=finalmergercheck(x1,x2,imgorig,wat31,wat41,gaboutt)
im=imgorig;
II1=double(im2uint8(imadjust(im2double(rgb2gray(imgorig)))));

% [Ir,map]=gray2ind(im2double(im(:,:,1)),16);
% [Ig,map]=gray2ind(im2double(im(:,:,2)),16);
% [Ib,map]=gray2ind(im2double(im(:,:,3)),16);
% LL=double(Ir)+16*double(Ig)+256*double(Ib);
oo1=II1(wat31==x1);
oo2=II1(wat31==x2);
cov1=std(oo1);
cov2=std(oo2);
% size(oo1)
% size(oo2)
mew1=mean(oo1);
mew2=mean(oo2);
II1=double((((rgb2gray(imgorig)))));
oo1=double(II1(wat31==x1));
%oo1=fix((oo1-min(min(oo1)))/(max(max(oo1))-min(min(oo1)))*255);
count1ww=hist(oo1,-50:400);
count1ww=count1ww/sum(count1ww);

% count1ww=(count1ww-1)/(4096-1);
oo2=double(II1(wat31==x2))+double(fix(mew1-mew2));
%oo2=fix((oo2-min(min(oo2)))/(max(max(oo2))-min(min(oo2)))*255);
count2ww=hist(oo2,-50:400);
count2ww=count2ww/sum(count2ww);

%bhattdww=sum(sqrt(count1ww.*count2ww))
diffbhatt=(((count1ww-count2ww).^2))./(count1ww+count2ww);
% abs(mew1-mew2)

for i=1:size(diffbhatt,2)
    if isnan(diffbhatt(i))
        diffbhatt(i)=0;
    end
end
% sum(diffbhatt)
%  (mew1-mew2)
%   sum(diffbhatt)
histodiff=sum(diffbhatt);
%var(cat(1,oo1/255,oo2/255))
mew1=mew1/255;
mew2=mew2/255;
cov1=cov1/(255);
cov2=cov2/(255);

intensitydiff=(mew1-mew2)*255;
ooo1=gaboutt(wat31==x1);
ooo2=gaboutt(wat31==x2);
textdiff=abs(mean(ooo1)-mean(ooo2));
intensitydiff=(mew1-mew2);%*(1-((std(cat(1,oo1/255,oo2/255)))^(mew1+mew2)/((cov1^mew1)*(cov2^mew2))));
% figure,plot(-50:400,count1ww);
% axis([-50 400 0 0.05]);
%figure,imshow(double(edge(double(wat41==x1)))+im2double(II1/255))
% figure,plot(-50:400,count2ww);
% axis([-50 400 0 0.05]);
%figure,imshow(double(edge(double(wat41==x2)))+im2double(II1/255))
% figure,plot(-50:400,diffbhatt);
end