clc,clear,close all
addpath('faces/')
faces_matrix = [];
imageSize = 40;
for i = 1:9
    im = imread(strcat('train',int2str(i),'.jpg'));
    im = rgb2gray(im);
    im = imresize(im,[imageSize,imageSize]);
    im_vec = double(im(:));
    faces_matrix = [faces_matrix,im_vec];
end
test_im = rgb2gray(imread('test9.jpg'));
test_im = imresize(test_im,[imageSize,imageSize]);
faces_matrix = faces_matrix- mean(faces_matrix,2);
test_vec = double(test_im(:)) - mean(faces_matrix,2);
[U,D,V] = svd(faces_matrix);
Ur = U(:,1:9);
features = D*V';
features = features(1:9,1:9);
new_feature = (test_vec'*Ur)';
errs = [];
for i = 1:9
    face_feature = features(:,i);
    error = norm(face_feature - new_feature);
    errs = [errs,error];
end
plot(errs,'*');
xlabel('face # in training set')
ylabel('Sqaured Error');