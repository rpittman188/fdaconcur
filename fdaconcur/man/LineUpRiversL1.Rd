\name{LineUpRiversL1}
\alias{LineUpRiversL1}
\title{LineUpRiversL1 Function}
\usage{
Uses L1 difference as the criteria for selecting the best starting and ending index of a new curve to best resemble the desired curve. Used to select the best time to start and stop a differing flood event to use in the Concurrent FDA model.
}
\description{
This function takes the Main, target, curve and a alternative event that should include a much longer period than one would expect to be necessary and gives the best indices from the alternative curve to start and stop the event so it best resembles the target curve. It starts at the first and last index, interpolates the curve to have the same length as the Main curve and then uses the trapezoid rule to calculate the area between the absolute difference (L1 difference) of these two curves. Then it trims the New Curve by 1 index from the front of the dataset and repeats the process to find a new L1 difference (call this "a"). It returns the 1 cut point from the front and cuts one off from the end of the event and finds the L1 difference again (call this "b"). Whichever is smaller, "a" or "b" is the edit that is made to the dataset. This process repeats until there are essentially no more values to cut off. Then, the starting and ending indices that give the lowest L1 difference are returned to show the user where to best trim the data. The entire process' results are also included to ensure there are not other local minuma that should be considered.
}
\examples{
NA 
}