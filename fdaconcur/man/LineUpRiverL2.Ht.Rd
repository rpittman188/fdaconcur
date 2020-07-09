\name{LineUpRiversL2.Ht}
\alias{LineUpRiversL2.Ht}
\title{LineUpRiversL2.Ht Function}
\usage{
Uses a modified L2 difference (L2*) as the criteria for selecting the best starting and ending index of a new curve to best resemble the desired curve. The only modification is instead of simply taking the squared difference between the two curves, it weights the difference based on the squared height of the main curve. It is used to select the best time to start and stop a differing flood event to use in the Concurrent FDA model, placing extra emphasis on predicting the higher values of the main curve.
}
\description{
This function takes the Main, target, curve and a alternative event that should include a much longer period than one would expect to be necessary and gives the best indices from the alternative curve to start and stop the event so it best resembles the target curve. It starts at the first and last index, interpolates the curve to have the same length as the Main curve and then uses the trapezoid rule to calculate the area between the modified squared difference of these two curves and then takes the square root (L2* difference). The modification is as mentioned above. Then it trims the New Curve by 1 index from the front of the dataset and repeats the process to find a new L2* difference (call this "a"). It returns the 1 cut point from the front and cuts one off from the end of the event and finds the L2* difference again (call this "b"). Whichever is smaller, "a" or "b" is the edit that is made to the dataset. This process repeats until there are essentially no more values to cut off. Then, the starting and ending indices that give the lowest L2* difference are returned to show the user where to best trim the data. The entire process' results are also included to ensure there are not other local minuma that should be considered.
}
\examples{
NA 
}