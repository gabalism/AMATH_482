clear,clc
%Michael Gabalis
% PCA for clustering MNIST Images
%% Loading in the data
[testimages,testlabels] = MNISTload('t10k-images.idx3-ubyte','t10k-labels.idx1-ubyte');
[images,labels] = MNISTload('train-images.idx3-ubyte','train-labels.idx1-ubyte');
%% Extracting different digits
a = 5;
b = 8;
c = 5; 
images = images(:,labels == a |labels ==b |labels==c);
%images = images(:,labels == a |labels ==b);
labels = labels(labels == a | labels ==b |labels ==c,:); 
%labels = labels(labels == a | labels ==b,:); 
testimages = testimages(:, testlabels==a |testlabels ==b |testlabels ==c);
testlabels = testlabels(testlabels==a | testlabels == b |testlabels ==c,:);
%testimages = testimages(:, testlabels==a |testlabels ==b);
%testlabels = testlabels(testlabels==a | testlabels == b,:);

%% Mean centering and pixe l covariance
[numpixels,numimages] = size(images);
%mean of each pixel
mu = mean(images,2);
%Normalization
X = images - repmat(mu,1,numimages);
images_centered = X/sqrt(numimages-1);


%% SVD and reconstructing the covariance
[U,S,V] = svd(images_centered,'econ');
eigs = diag(S*S'/(numimages-1))';

%% Variance explained per Principal Component

var = 100* cumsum(eigs)./sum(eigs);
figure(1)
subplot(1,2,1)
plot(diag(S)./ max(diag(S)),'r.','markersize',15);
title('Normalized Singular Values')
xlabel('index ')
ylabel('\sigma')

subplot(1,2,2)
plot(1:numpixels, var);
title('Cumulative variance explained')
xlabel('Eigenvalue #')
ylabel('Cumulative % Variance')

%% Finding the energy of the system
energy = 0;
total = sum(diag(S));
% how much energy we want our modes to capture.
% 75% does alright; 90% does very good.
threshold = 0.90; 
r = 0;
while energy < threshold
    r = r + 1;
    energy = energy + S(r,r)/total;
end
rank = r;
Proj = U(:,1:rank)'*images_centered;
%% Vizualizing the eigenvectors as images

figure(3)
for k = 1:9
    subplot(3,3,k)
    imagesc(reshape(U(:,k),28,28));
    colormap gray;
end

%% Rotate mean subtracted data and vizualize projections

figure

scatter3(Proj(1,labels==0), Proj(2,labels==0),Proj(3,labels==0), 'b','filled')
hold on

scatter3(Proj(1,labels==1), Proj(2,labels==1),Proj(3,labels==1), 'r','filled')
hold on

scatter3(Proj(1,labels==2), Proj(2,labels==2),Proj(3,labels==2), 'k','filled')
legend('2s', '3s', '5s')
xlabel('PC1', 'fontsize', 20)
ylabel('PC2', 'fontsize', 20)
zlabel('PC3', 'fontsize', 20)
hold off
box off

%% Use LDA on only two digit numbers
first = Proj(:,labels==1);
second = Proj(:,labels==0);
[U2,S2,V2,threshold,w,sortones,sorttwos]=LDA_trainer(first,second,rank);

Testnum = size(testimages,2);
Proj3 = U2(:,1:rank)'*(testimages(1:rank,:));  
pval = w'*Proj3;

Res = (pval>threshold);
err = abs(Res'-testlabels);
errNum = sum(err);
success = 1-errNum/Testnum;
%% Use LDA on two or more digits, plot error
Proj2 = U(:,1:rank)'*(testimages);
Mdl = fitcdiscr(Proj',labels, 'DiscrimType','pseudolinear');
labeled = predict(Mdl,Proj2');
confusionchart(testlabels,labeled)
title('LDA Classification')
set(gca,'fontsize',20)
%% Perform SVM
SVM = fitcecoc(images', labels);
test_labels = predict(SVM,testimages');

confusionchart(testlabels,test_labels)
title('SVM Classification')
set(gca,'fontsize',20)
%% Perform Decision Tree
Tree = fitctree(Proj', labels,'MaxNumSplits',3,'CrossVal','on');
%test_labels = predict(Tree,Proj2');
%confusionchart(testlabels,test_labels)
classError = kfoldLoss(Tree)

