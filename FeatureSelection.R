#Set seed for reproducibility
set.seed(7)
library(mlbench)
library(caret)
library(tidyverse)

#Read in the data
mydata=read.csv("LFQ_Data_Num.csv", stringsAsFactors = FALSE)
#Create a vector that contains "pos" or "neg" for each row to indicate presence or lack of cancer
cancerPosNeg=unlist(mydata[1,], use.names=FALSE)
#Transpose data frame so that the rows are patients and the columns are proteins
flipped=as.data.frame(t(mydata), stringsAsFactors = FALSE)
proteinNamesTable=read.csv("ProteinNames.csv")
proteinNames=as.vector(proteinNamesTable[,1])
colnames(flipped)=proteinNames
#Ensure that all the data is numeric in the data frame and convert to a matrix
flipped[]=lapply(flipped[,2:6795], as.numeric)
#Calculate correlation matrix
correlationMatrix <- cor(flipped[,2:6795])
#Remove rows and columns that are entirely made of zeroes
correlationMatrix=correlationMatrix[-c(943,980,984,1109,1479,1548,1927,2511,2707,3366,3510,4095,4152,4252,4341,4366,4533,5093,5148,5993,6007),]
correlationMatrix=correlationMatrix[,-c(943,980,984,1109,1479,1548,1927,2511,2707,3366,3510,4095,4152,4252,4341,4366,4533,5093,5148,5993,6007)]
#Indentify the highly correlated features. These will be removed
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75, verbose=FALSE)
#Remove highly correlated features
withoutHighlyCorrelated=flipped[,-highlyCorrelated]
#Train a control
control=trainControl(method="repeatedcv", number=10, repeats=3)
#Add the cancer presence data back to the data frame
withCancerData=cbind(cancerPosNeg, withoutHighlyCorrelated)
#identify start and end column
i=1
j=2580
#Train the model using ROC curve method
model=train(cancerPosNeg~., data=withCancerData[,i:j], method="rocc", trControl=control)
#Rank by importance
importance=varImp(model, scale=FALSE)
print(importance)
