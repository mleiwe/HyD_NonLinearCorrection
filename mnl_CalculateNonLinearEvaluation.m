function [modelfunc,beta,OverallRMSE]=mnl_CalculateNonLinearEvaluation
%This function reads in a multi-stack tiff which contains the same image at
%increasing laser powers. It assumes the first 5 time points are linear and
%will produce a function and covariates which can be used in
%mnl_CorrectImageForNonLinearity to produce a corrected image.
%% Input: mNG-ch Tiff 512*512*15
[ImageID1,ImagePath1]=uigetfile('*.tif','Select the serial tiff file');
FileTif1=num2str([ImagePath1,ImageID1]);
InfoImage=imfinfo(FileTif1);
Data=tiffreadVolume(FileTif1);
LaserXs=linspace(0.2,3,15);%THe laser powers used

%% Reshape the data
nSteps=size(Data,3); %Added in nSteps instead of the hard coded value 15 for future use
Data = reshape(Data,[],nSteps);
idx = Data(:,15) <1000;
obs = Data(~idx,:);
obs=double(obs); % I added this in to make sure obs is a double and operations can be performed on it
Data = double(Data);
% View the plots
figure('Name','Intenisity changes with laser power')
for i=1:length(obs)
    if obs(i,1)>500 %only if the signl is strong
        plot(LaserXs,obs(i,:))
        hold on
    end
end
xlabel('Laser Power (%)')
ylabel('Intensity (AU)')
xlim([min(LaserXs) max(LaserXs)])

%% Calculate/Predict what the values should be if everything is linear
PredLinearObs=zeros(size(obs)); %Pre-allocate values
MaxLinearRange=5;
LinearLasers=LaserXs(1:MaxLinearRange);
k=nan(length(obs),1);

for i=1:length(obs)
    tempY=obs(i,1:MaxLinearRange); 
    k(i)=tempY/LinearLasers; % k=A\B solves the system of linear equations A*k = B.
    %Now predict the values if it is completely linear
    predY=LaserXs*k(i);
    PredLinearObs(i,:)=predY;
end
%Alternative method k for pixels containing signal only - to write
SignalIdx=obs(:,1)>=1000;
SignalK=k(SignalIdx);
    
%Visualisations
figure('Name','Evaluation of Differences to Linear')
subplot(2,6,3)
histogram(k,200,'FaceColor',[0,0,1],'EdgeAlpha',0)
title('Variation in  k')
xlabel('"k" value')

subplot(2,6,4)
histogram(SignalK,50,'FaceColor',[1,0,0],'EdgeAlpha',0)
title('Variation in  k in different pixels above 1000')
xlabel('"k" value')

subplot(2,3,5)
tempKs=nan(length(obs),2);
tempKs(:,1)=k;
tempKs(1:length(SignalK),2)=SignalK;
mnl_CumulativePlot4(tempKs)
legend('All Pixel Ks','Pixels Above Threshold')
xlabel('"k" value')

subplot(1,3,1)
for i=1:length(obs)
    if obs(i,1)>1000
        scatter(LaserXs,obs(i,:),'.b')
        hold on
    end
end
xlabel('Laser Power (%)')
ylabel('Intensity (AU)')
xlim([min(LaserXs) max(LaserXs)])
title('Intensity vs laser power - Raw')

subplot(1,3,3)
for i=1:length(obs)
    if obs(i,1)>=1000
        scatter(obs(i,:),PredLinearObs(i,:),'.b')
        hold on
    end
end
%Now plot a linear line between them
plot([min(obs(:)) max(obs(:))],[min(obs(:)) max(obs(:))],'--k')
xlabel('Original Intensity Value')
ylabel('Predicted Linear Intensity Value')
title({'Comparing Actual to predicted values', 'if truly linear'})
%% Reshape the values
TrueYs = reshape(obs,[],1);
PredYs = reshape(PredLinearObs,[],1);
%% Randomly sample evenly across the range
%Split into Sections of 10 percentiles to the nearest 1,000
maxValueY=floor(max(TrueYs(:)/1000))*1000;
FullMaxValue=ceil(max(TrueYs(:)/500))*500;
Divisions=linspace(0,maxValueY,11); 
Ystruct=struct('MinVal',[],'MaxVal',[],'nVals',[],'TrueValues_All',[],'PredValues_All',[],'TrueValues_Sample',[],'PredValues_Sample',[]);
for i=1:10
    LowLim=Divisions(i);
    HighLim=Divisions(i+1);
    t_idx= TrueYs>=LowLim & TrueYs<HighLim;
    %Now add to the structure
    Ystruct(i).MinVal=LowLim;
    Ystruct(i).MaxVal=HighLim;
    Ystruct(i).nVals=sum(t_idx);
    Ystruct(i).TrueValues_All=TrueYs(t_idx);
    Ystruct(i).PredValues_All=PredYs(t_idx);
    All_nVals(i)=sum(t_idx);
end
%Now resample evenly across all percentiles
N=min(All_nVals);
TrueYs_Resample=nan(10*N,1);
PredYs_Resample=nan(10*N,1);
for i=1:10
    t_idx=randsample(Ystruct(i).nVals,N,false);
    Ystruct(i).TrueValues_Sample=Ystruct(i).TrueValues_All(t_idx);
    Ystruct(i).PredValues_Sample=Ystruct(i).PredValues_All(t_idx);
    TrueYs_Resample(((i-1)*N+1):(i*N))=Ystruct(i).TrueValues_Sample;
    PredYs_Resample(((i-1)*N+1):(i*N))=Ystruct(i).PredValues_Sample;
end

%Visualise this
figure('Name','Randomly Sampled')
scatter(TrueYs_Resample,PredYs_Resample,'.b')
hold on
%Now plot a linear line between them
plot([min(obs(:)) max(obs(:))],[min(obs(:)) max(obs(:))],'--k')
xlabel('Original Intensity Value')
ylabel('Predicted Linear Intensity Value')
title({'Comparing Actual to predicted values', 'if truly linear'})
%% Fit the model
plotXs=0:FullMaxValue;
figure('Name','Potential fit')
%((b1).*exp(b(2).*b(3))
modelfunc= @(b,x) (b(1).*x)+(b(2).*(exp(b(3).*x)))-b(2);
beta0= [1,1,0.001];
beta = nlinfit(TrueYs_Resample, PredYs_Resample, modelfunc, beta0);
CorrYs=modelfunc(beta,plotXs);
%Plot the corrected function
subplot(1,2,1)
scatter(TrueYs_Resample,PredYs_Resample,'.b')
hold on
plot(plotXs,CorrYs,'r','LineWidth',2)
title('((b1).*exp(b(2).*b(3)')
xlabel('Original Intensity Value')
ylabel('Predicted Intensity Value')
%Calculate the Residuals
FitYs=modelfunc(beta,TrueYs);
ResidualYs=PredYs-FitYs;
subplot(1,2,2)
scatter(TrueYs,ResidualYs,'.k')
title('Residuals')
xlabel('Original Intensity Value')
ylabel('Residual Error')
%RMSE value
OverallRMSE=mean(sqrt(sum(ResidualYs.^2)));
% Split into 500 zones
Zones=0:500:FullMaxValue;
nZones=length(Zones);
ScaledResidualYs=ResidualYs./TrueYs;
AllResiduals=nan(100,nZones-1);
ScaledAllResiduals=nan(100,nZones-1);
for i=1:nZones-1
    LowerLim=Zones(i);
    UpperLim=Zones(i+1);
    t_idx=TrueYs>=LowerLim & TrueYs<UpperLim;
    SelectedResiduals=ResidualYs(t_idx);
    SelectedScaled=ScaledResidualYs(t_idx);
    %Assign to matrix
    AllResiduals(1:length(SelectedResiduals),i)=SelectedResiduals;
    ScaledAllResiduals(1:length(SelectedResiduals),i)=SelectedScaled;
    %Labels
    ColName{i}=sprintf('%d%s%d',LowerLim,' to ',UpperLim);
end
%Find Zeros and Replace with Nans
tidx=AllResiduals==0;
AllResiduals(tidx)=nan;
ScaledAllResiduals(tidx)=nan;
%Now plot boxplots
figure('Name','Spread of Residuals per zone')
subplot(2,3,1)
mnl_boxplot2(AllResiduals,ColName,'Residual','y','y')
title('Residuals (signed)')
subplot(2,3,2)
SqAllResiduals=sqrt(AllResiduals.^2);
mnl_boxplot2(SqAllResiduals,ColName,'Residual','y','y')
title('Residuals (unsigned)')
subplot(2,3,3)
SqScaledAllResiduals=sqrt(ScaledAllResiduals.^2);
mnl_boxplot2(SqScaledAllResiduals,ColName,'Residual','y','y')
title('Residuals (unsigned and scaled)')
subplot(2,3,4)
mnl_boxplot2(AllResiduals,ColName,'Residual','y','y')
title('Residuals (signed)')
ylim([-1000 1000])
subplot(2,3,5)
SqAllResiduals=sqrt(AllResiduals.^2);
mnl_boxplot2(SqAllResiduals,ColName,'Residual','y','y')
title('Residuals (unsigned)')
ylim([0 1000])           
subplot(2,3,6)
SqScaledAllResiduals=sqrt(ScaledAllResiduals.^2);
mnl_boxplot2(SqScaledAllResiduals,ColName,'Residual','y','y')
title('Residuals (unsigned and scaled)')
ylim([0 1])
end