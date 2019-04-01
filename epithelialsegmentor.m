function finalclassesimage=epithelialsegmentor(inputimage)
% imagename : name of the input oral mucosa image
% finalclassesimage : final image containing the various regions in different classes
% usage : finalclassesimage=epithelialsegmentor('N1-09-6.bmp');
% This file requires the matitk library : copy the 3 files in matitk.rar within the folder where you run this code
% This file also require imRAG.m(for finding region adjacency matrix) & finalmergercheck.m(for final stage of region merging)
clc;
imagename=inputimage(1:size(inputimage)-4);
imgorig=uint8(imread(inputimage));
cc=colorspace('Lab<-RGB',(imgorig));
L = cc(:,:,1)/100; 
cc(:,:,1) = imadjust(L)*100;
imgorig=colorspace('RGB<-Lab',cc);
img=wiener2((rgb2gray(imgorig)),[5 5]);
counter=0;
I=im2double(img);
i=0;
nooffrequencies=4;
noofangles=6;
gabout=single(zeros(size(img,1),size(img,2),nooffrequencies*noofangles));
for countf=6:9
    f=2^(countf-1)*sqrt(2)/1024;
    b=1;
    Sx=1/f*(1/pi)*(sqrt(log(2)/2)*(2^b+1)/(2^b-1))
    Sy=(sqrt(log(2)))/(sqrt(2)*pi*f*tan(pi/(2*noofangles)))
    windowsize=fix(2*(Sx+2));
    for countt=0:noofangles-1
        counter=counter+1;
        theta=pi/noofangles*countt+pi/2;
        if isa(I,'double')~=1
            I = double(I);
        end
        j = sqrt(-1);
        clear G
        for x = -windowsize:windowsize
            for y = -windowsize:windowsize
                xPrime = x * cos(theta) + y * sin(theta);
                yPrime = y * cos(theta) - x * sin(theta);
                G(windowsize+x+1,windowsize+y+1) = (exp(-.5*((xPrime/Sx)^2+(yPrime/Sy)^2))*exp(2*i*pi*f*xPrime));
                Ggauss(windowsize+x+1,windowsize+y+1)=(exp(-.5*((xPrime/Sx)^2+(yPrime/Sy)^2)));
            end
        end
        realG=real(G)-mean2(real(G));
        imagG=imag(G)-mean2(imag(G));
		Imgabout = conv2(I,double((imagG)),'same');
        Regabout = conv2(I,double((realG)),'same');
		gabout(:,:,counter) = sqrt(Regabout.^2+Imgabout.^2);
		gabout(:,:,counter)=(gabout(:,:,counter)-min(min(gabout(:,:,counter))))/(max(max(gabout(:,:,counter)))-min(min(gabout(:,:,counter))));
    end
end
v=1;
save (strcat(imagename,'gabor'),'gabout');
for v=1:nooffrequencies*noofangles
    gabout2(:,:,v)=tanh(0.25*gabout(:,:,v));
end
clear Regabout Imgabout img1 gabout
v=1;
for countf=6:9
f=(2^(countf-1))*sqrt(2)/1024;
  for countt=0:noofangles-1
    Sx=0.5*1/(f)
    gabout2(:,:,v)=imfilter(gabout2(:,:,v),fspecial('gaussian',[fix(2*Sx) fix(2*Sx)],Sx));
    gabout2(:,:,v)=(gabout2(:,:,v)-min(min(gabout2(:,:,v))))/(max(max(gabout2(:,:,v)))-min(min(gabout2(:,:,v))));
    gabout3=single(gabout2(:,:,v));
    [gaboutx(:,:,v),gabouty(:,:,v)]=gradient(imadjust(gabout3));
    v=v+1;
  end
end

counter=v-1;
for v=1:counter
gaboutx(:,:,v)=gaboutx(:,:,v)/max(max(gaboutx(:,:,v)));
end
for v=1:counter
gabouty(:,:,v)=gabouty(:,:,v)/max(max(gabouty(:,:,v)));
end
clear gabout3
for i=1:size(gaboutx,1)
    for j=1:size(gabouty,2)
        e=squeeze(gaboutx(i,j,:));
        g=squeeze(gabouty(i,j,:));
        E=e'*e;
        F=e'*g;
        G=g'*g;
        gradienttexture(i,j)=sqrt((E-G)^2+4*(F^2));
    end
end
hy = fspecial('sobel');
hx = hy';
% 
for i=1:3
    orig(:,:,i)=wiener2(imgorig(:,:,i),[5 5]);
end
cc=colorspace('Lab<-RGB',orig);
for i=1:3
Iy(:,:,i) = imfilter(double(cc(:,:,i)), hy, 'replicate');
Ix(:,:,i) = imfilter(double(cc(:,:,i)), hx, 'replicate');
end
gradmag = sqrt(sum(Ix.^2,3) + sum(Iy.^2,3));
for i=1:4
    gradimage(:,:,i)=(gradmag(1:size(img,1),1:size(img,2))/40)+imadjust((gradienttexture),[0 0.5],[0 1]);
    gradimage(:,:,i)=gradimage(:,:,i)/2;
end
save (imagename,'gradimage');
save (strcat(imagename,'gradtex'),'gradienttexture');
figure,imshow(gradimage(:,:,1));
grim1=gradimage(:,:,1)+20/255; grim=gradimage(:,:,1);
se=strel('disk',1);
for i=1:15
    grim1=imerode(grim1,se);
    for j=1:size(grim,1)
        for k=1:size(grim,2)
            grim1(j,k)=max([grim(j,k) grim1(j,k)]);
        end
    end
end
figure,imshow(grim1)
for i=1:4
    gradim(:,:,i)=(grim1);
end
Water=matitk('SWS',[0.1 0.0],gradim);  
figure,imagesc(Water(:,:,4));
labimg=0*Water(:,:,4);
wat=Water(:,:,4);
kk=unique(Water(:));
counter=0;
warning off all
for i=1:size(kk,1)
    counter=counter+1;
    wat(wat==kk(i))=counter;
end
% Region merging starts here %
%1st part of region merging based on RGB color space

while(1)
wat1=wat;
clear watx waty
watpad=padarray(wat,[1 1],0);
watx=abs(diff(watpad,1,1));
waty=abs(diff(watpad,1,2));
watx1=watx(1:size(watx,1)-1,2:size(watx,2)-1);
waty1=waty(2:size(waty,1)-1,1:size(waty,2)-1);
wat1((sign(watx1+waty1)>0))=0;
IIdash=rgb2gray(im2double(imgorig));
II=(IIdash);
colorcount=0;
adj=imRAG(wat1);
clear I1 I2
	for k=1:size(adj,1)
		for i=1:3
        II=imgorig(:,:,i);
        I1(:,i)=double(II(wat==adj(k,1)));
        I2(:,i)=double(II(wat==adj(k,2)));
        end
        n1=size(I1,1);
        n2=size(I2,1);
        if (n1<100000 && n2<100000)
            mew1=mean(I1,1);
            mew2=mean(I2,1);
            cov1=cov(I1);
            cov2=cov(I2);
            n1=size(I1,1);
            n2=size(I2,1);
            S=(n1*cov1+n2*cov2)/(n1+n2);
            T=(mew1-mew2)*inv((1/n1+1/n2)*S)*(mew1-mew2)';
            if (T<=300)
                wat(wat==adj(k,2))=adj(k,1);
                adj(adj==adj(k,2))=adj(k,1);
                colorcount=colorcount+1;
            end
        end
        clear I1 I2
    end
	if colorcount==0
		break;
	end
end
save(strcat(imagename,'firstregionmerging'),'wat');
%Calculation of combined texture feature for region merging
%clear various variables
clear Imgabout Regabout gabout x y xPrime yPrime f gsize imgorig img img1 I L cc
load([imagename 'gabor']);
i=uint8(0);
p=uint8(0);
for i=1:24
    padgabout=gray2ind(imadjust(gabout(:,:,i)),32);
    for p=1:32
        gaboutcnt(:,:,p)=single(conv2(single(padgabout==(p-1)),ones(21,21),'same'));
    end
    phor=[0 0 0;0 1 -1;0 0 0];
    nhor=[0 0 0;-1 1 0;0 0 0];
    nver=[0 0 0;0 1 0;0 -1 0];
    pver=[0 -1 0;0 1 0;0 0 0];
    for g=1:32
        gaboutdf(:,:,g)=single(conv2(gaboutcnt(:,:,g),phor,'same'));
    end
    gaboutd(:,:,1)=sum(abs(gaboutdf),3);
    for g=1:32
        gaboutdf(:,:,g)=single(conv2(gaboutcnt(:,:,g),nhor,'same'));
    end
    gaboutd(:,:,2)=sum(abs(gaboutdf),3);
    for g=1:32
        gaboutdf(:,:,g)=single(conv2(gaboutcnt(:,:,g),pver,'same'));
    end
    gaboutd(:,:,3)=sum(abs(gaboutdf),3);
    for g=1:32
        gaboutdf(:,:,g)=single(conv2(gaboutcnt(:,:,g),nver,'same'));
    end
    gaboutd(:,:,4)=sum(abs(gaboutdf),3);
    gabouttex=sum(gaboutd,3)/441;
    save([imagename 'gabouttex' num2str(i)],'gabouttex');
end
clear gaboutdf
clear gaboutcnt
for i=1:24
load([imagename 'gabouttex' num2str(i)]);
gaboutt(:,:,i)=gabouttex;
end
load([imagename 'firstregionmerging']);
II=(im2double(rgb2gray(imgorig)));

%2nd part of region merging
[Ir,map]=gray2ind(imgorig(:,:,1),16);
[Ig,map]=gray2ind(imgorig(:,:,2),16);
[Ib,map]=gray2ind(imgorig(:,:,3),16);
LL=double(Ir)+16*double(Ig)+256*double(Ib);
wat4=wat;

for m=1:6
    mergecount=0;
    wat3=wat4;
	watpad3=padarray(wat3,[1 1],0);
	wat3x=abs(diff(watpad3,1,1));
	wat3y=abs(diff(watpad3,1,2));
	wat3x1=wat3x(1:size(wat3x,1)-1,2:size(wat3x,2)-1);
	wat3y1=wat3y(2:size(wat3y,1)-1,1:size(wat3y,2)-1);
	wat3((sign(wat3x1+wat3y1)>0))=0;
	merged=zeros(1000,1);
	adjwat3=imRAG(wat3);
	for k=1:size(adjwat3,1)
       if size(wat4(wat4==adjwat3(k,1)),1)>1000 && size(wat4(wat4==adjwat3(k,2)),1)>1000
       [histdiffo,intensitydiff,texdiff]=finalmergercheck(adjwat3(k,1),adjwat3(k,2),imgorig,wat4,wat4,sqrt(mean(gaboutt.^2,3)));
        if (histdiffo<=0.2 && abs(intensitydiff)<=0.1 && texdiff<0.04)
			mergecount=mergecount+1;
			wat4(wat4==adjwat3(k,2))=adjwat3(k,1);
			adjwat3(adjwat3==adjwat3(k,2))=adjwat3(k,1);
		end
	   end
 	end
if mergecount==0
    break;
end
end
save([imagename '_finalclassesimage'],'wat4');
finalclassesimage=wat4;
% g=wat4;  
% cd ..
% imgorig=(imread('1/N1-09-6.bmp'));
% cd sem9
% img1=imgorig;
% img11=im2bw(img1,0.5);
% for i=1:size(img1,1)
% for j=1:size(img1,2)
% if g(i,j)==16
% img2(i,j)=1;
% else 
%     img2(i,j)=0;
% end
% end
% end
% img4=(im2double(logical((img2))));
% img5=imdilate(img4,strel('diamond',0));
% %img5=imopen(img5,strel('disk',10));
% figure,imshow(img5);
% img6=bwareaopen(img5,20000);
% 
% 
% % figure,imshow(img6);
% img7=imerode(img6,strel('diamond',0));
% %figure,imshow(img7/6+img);
% img8=bwareaopen(img7,0);
% %figure,imshow(img7/6+img+edge(im2double(img7)));
% img8=imclose(img8,strel('disk',8));
% 
% roi(:,:,1)=img8;
% roi(:,:,2)=img8;
% roi(:,:,3)=img8;
% roiedge(:,:,1)=edge(im2double(img8));
% roiedge(:,:,2)=edge(im2double(img8));
% roiedge(:,:,3)=edge(im2double(img8));
% % img1=((uint8((imread('testimage1.bmp')))));
% imgg=img1(1:size(img1,1),1:size(img1,2),:);
% figure,imshow(im2double(imgg)+im2double(roiedge(:,:,:))+im2double(roi(:,:,:))/10)
% for i=1:3
% imgmid=imgg(:,:,i);
% imgmid(img8==0)=0;
% imgg(:,:,i)=imgmid;
% end
% figure,imshow(imgg);
% imgg=img1(1:size(img1,1),1:size(img1,2),:);
% cd ..
% imggg=otsu(rgb2gray(imgg),3);
% for i=1:3
% imgmid=imgg(:,:,i);
% imgmid(imggg>0)=0;
% imgg(:,:,i)=imgmid;
% end
% figure,imshow(imgg);
% % figure,imshow(im2double(img1)+im2double(roiedge)+im2double(roi)/30);
