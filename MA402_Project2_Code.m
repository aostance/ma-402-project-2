% Andrea Stancescu, Jenna Varnell, Omry Brewster
% MA 402 Project 2 
% 5-2-2022
% MATLAB .mlx file with integrated text and code

load jokesSp22.mat
%%

% a) Find the jokes with the highest and lowest ratings over JTrain.
% aka find the mean of each column (joke) and find smallest and largest average / mean 
 avg = mean(JTrain,1);
% return 42 length vector with the mean score of each columns / joke 
[score_f, funniest] = max(avg);
fprintf("The funniest joke is the %s joke which says ''%s'' \n", char(jokes(funniest)), char(jokeTexts(funniest)))

[score_l, nonfunny] = min(avg);
fprintf("The least funny joke is the %s joke which says ''%s'' \n", char(jokes(nonfunny)), char(jokeTexts(nonfunny)))

person = mean(JTrain, 2);
% return a 63 x 1 length vector with mean score of each individual
[score_p, highestperson] = max(person);
fprintf("The individual who gave the highest average joke rating is the %d person with an avergae joke rating of %f \n", highestperson , score_p)

% find which joke the 11th person liked the least
person_data = JTrain(11,:);
overallworst_score = min(person_data);
[worst_scores, worst_jokes] = find(person_data == overallworst_score);
fprintf("The individual who gave the highest average joke rating gave two low ratings of %d. \nThe jokes that recieved these ratings were: \nthe %s joke which says ''%s'' \nand \nthe %s joke which says ''%s'' \n ", overallworst_score, char(jokes(worst_jokes(1))), char(jokeTexts(worst_jokes(1))), char(jokes(worst_jokes(2))), char(jokeTexts(worst_jokes(2))))  



figure(1)
pcolor(JTrain)
colorbar

figure(2)
individ30 = JTrain(30,:);
emptyrow = zeros(1,42);
pcolor([individ30 ; emptyrow])
colorbar


figure(3)
individ11 = JTrain(11, :);
emptyrow = zeros(1, 42);
pcolor([individ11 ; emptyrow])
colorbar

rho = corrcoef(JTrain);
% maximum possible correlation is 1
for i = 1:42
    rho(i, i) = 0;
end 
colmax = max(rho);
largesttwo = maxk(colmax, 4);



[r1,c1] = find(rho == largesttwo(1));
% joke 20 and 22 have the highest correlation 
fprintf("The 2 jokes with the highest correlation are jokes %d and %d with a correlation of %f \nThis is the %s joke which says ''%s'' \nand the %s joke which says ''%s''. \n ", r1(1), r1(2), rho(r1(1),r1(2)), char(jokes(r1(1))), char(jokeTexts(r1(1))), char(jokes(r1(2))), char(jokeTexts(r1(2))) )


fprintf("\n")
[r2, c2] = find(rho == largesttwo(3));
% joke 11 and 12 have the second highest correlation 
fprintf("The 2 jokes with the second highest correlation are jokes %d and %d with a correlation of %f \nThis is the %s joke which says ''%s'' \nand the %s joke which says ''%s''. \n ", r2(1), r2(2), rho(r2(1),r2(2)), char(jokes(r2(1))), char(jokeTexts(r2(1))), char(jokes(r2(2))), char(jokeTexts(r2(2))) )

rho2 = corrcoef(JTrain');
for i = 1:63
    rho2(i, i) = 0;
end 
colmax = max(rho2);
largesttwo = maxk(colmax, 4);



[r1,c1] = find(rho2 == largesttwo(1));
% person 5 and 32 have the highest correlation 
fprintf("The 2 individuals with the highest correlated responses are person %d and %d with a correlation of %f  \n ", r1(1), r1(2), rho2(r1(1),r1(2)) )


fprintf("\n")
[r2, c2] = find(rho2 == largesttwo(3));
% person 50 and 53 have the second highest correlation 
fprintf("The 2 individuals with the second highest correlated responses are person %d and %d with a correlation of %f  \n ", r2(1), r2(2), rho2(r2(1),r2(2)) )

% scatterplot of individuals' 5 and 32 ratings
figure(4)
hold on 
scatter(JTrain(5,:), JTrain(32,:))
xlabel("Individual's 5 Ratings")
ylabel("Individual's 32 Ratings")
title("Correlation between the individual 5 and 32")
hold off

%scatterplot of individuals' 50 and 53 ratings 
figure(5)
hold on
scatter(JTrain(50,:), JTrain(53,:))
xlabel("Individual's 50 Ratings")
ylabel("Individual's 53 Ratings")
title("Correlation between the individual 50 and 53")
hold off







%%

% k means clustering method

% a)Did you reduce the dimensionality of the data set first? If so, how many dimensions did you retain?
% below I am standardizing the data to have a mean of zero and standard deviation of 1
JTrainreduced = normalize(JTrain');
M = 2; % principal number of components to retain
[~,JTrainreduced] = pca(JTrainreduced,'NumComponents',M); 

% b) How many clusters did you end up using? How did you decide on this number?
% use elbow plot to choose k (optimal amount of clusters)

for k = 1:42
    %fprintf("k = %d\n", k); 
    [~,~,sumd] = kmeans(JTrainreduced,k);
    errs(k) = sum(sumd);
end 
%[idx3,C,sumdist3] = kmeans(JTrain,3)

figure(6)
plot(errs)
xlabel("Number of Clusters", 'FontSize',14); 
ylabel("Sum of Distances", 'FontSize',14); 
title("JTrain Data",'FontSize',14); 


% optimal k = 7, 19
% optimal k is observed where there is a bend / elbow in the plot

% c) Do the jokes that were clustered together seem related to you?

figure(7)
[idx,C] = kmeans(JTrainreduced, 7)

f(7) = makekmeansplot(JTrainreduced, 7);
title("Kmeans with k = 7")


figure(8)
[idx,C] = kmeans(JTrainreduced, 19)

f(19) = makekmeansplot(JTrainreduced, 19);
title("Kmeans with k = 19")




%%

X = JTrain(1:63, 1:20);
Y = JTrain(1:63, 21:42);
B = X\Y; %20 x 22 - each column gives the parameters for a line of best fit for one of last 22 jokes


predictions = JTest * B;
for i = 1: (36*22)
    if predictions(i) < 0
        predictions(i) = 0;
    elseif predictions(i) > 10
        predictions(i) = 10;
    end
end

% calculate error - compare to the average of the the training set for each joke is a 1x22 vector
avglastjokes = mean(Y, 1);
avgpredictions = mean(predictions, 1);
% error between average predictions - measured in frobenius norm
avgerr = norm(avglastjokes - avgpredictions, 'F')




figure(9)
scatter(avglastjokes, avgpredictions)
lsline
xlabel("Average score of the training jokes (with labeled measurements)")
ylabel("Average score of the testing jokes (unlabeled data)")
title("Comparison of predicted mean scores versus measured mean scores")

% the plot shows high correlation and frobenius error shows small error values suggesting that the scores were predicted well (in comparison with previous years)


train = JTrain(:, 21:42);
test = predictions;

figure(10)
sgtitle("Scatterplots of all the measured scores (green) and predicted scores (red)")
for i = 1:22
    subplot(4,6,i)
    hold on
    scatter(1:63, train(:, i), 'g')
    scatter(1:36, test(:, i), 'r')
    
    hold off
end 





trainx = JTrain(1:63, 1:20);
trainy = JTrain(1:63, 21:42);

testx = JTest;
% reducing the dimension via pca reduces the noise in the predictions
[COEFF,SCORE] = pca(trainx,'numComponents',2);

SCORE_test = (testx - mean(trainx))*COEFF; 
A = [ones(size(SCORE,1),1) , SCORE];
B = trainy;
C = A\B; % least squares matrix of coefficients -- size 3 x 22
testy = C(1,:) + (testx - mean(trainx))*COEFF * C(2:end, :);  


for i = 1:(36*22)
    if testy(i) < 0
        testy(i) = 0;
    elseif testy(i) > 10
        testy(i) = 10;
    end
end
% 
% calculate error - compare to the average of the the training set for each joke is a 1x22 vector
avgtraining = mean(trainy, 1);
avgtesting = mean(testy, 1);
% error between average predictions - measured in frobenius norm
avgerr = norm(avgtraining - avgtesting, 'F')


figure(11)
scatter(avgtraining, avgtesting)
lsline
xlabel("Average score of the reduced training jokes (with labeled measurements)")
ylabel("Average score of the testing jokes (unlabeled data)")
title("Comparison of predicted mean scores versus measured mean scores after dimensionality reduction")












