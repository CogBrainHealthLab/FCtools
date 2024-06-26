% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CPM.R
\name{cpm.lesion}
\alias{cpm.lesion}
\title{cpm.lesion}
\usage{
cpm.lesion(train.data, test.data, train.outcome, test.outcome, p = 0.05)
}
\arguments{
\item{train.data}{An N x E matrix containing the vectorized edges in the training dataset; where N = number of subjects, E=number of edges}

\item{test.data}{An N x E matrix containing the vectorized edges in the training dataset}

\item{train.outcome}{The outcome variable to predict, within the training dataset}

\item{test.outcome}{The outcome variable to predict, within the testing dataset}

\item{p}{The p-value threshold of a Pearson's correlation test between the feature and outcome that determines if the feature is selected. Set to 0.05 by default.}
}
\value{
Returns a matrix with the following columns
\itemize{
\item \code{lesion.model} The CPM model in which the listed network is left out
\item \code{positive} The predicted-actual correlations for each of the \code{lesion.model}s' positive network model
\item \code{negative} The predicted-actual correlations for each of the \code{lesion.model}s' negative network model
\item \code{both} The predicted-actual correlations for each of the \code{lesion.model}s' combined positive + negative network model
}
}
\description{
CPM with the leave-one-network-out lesion approach
}
\details{
This function runs the \code{cpm.train()} and \code{cpm.predict()} while leaving out a network of edges each time.
The changes in predicted-actual correlation while a network being left out can be interpreted to be an
indication of the contribution of the network to predicting the outcome
}
\examples{
\dontrun{
model1.lesion=cpm.lesion(train.data=FC_data.train,test.data=FC_data.test,train.outcome=train_dat$age, test.outcome=test_dat$age,p=0.05)
}

}
