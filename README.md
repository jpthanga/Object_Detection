# __Assignment 1 -Report__

### __Part 1 - Watermarking__

#### Parameters

1. Alpha = 11
2. L (length of binary vector) = random with a modulo 37 applied to it.
3. Radius = <\row size> - 40
4. Threshold = 0.2

#### Add watermark  

_./watermark 1.3 </name of input image> </output image> add N_'''

* Using a circle of the aforementioned radius we modify "l" real value parts of the DFT of the image. Here "l" is the length of the binary vector we randomly generate using the random function in C. Essentially, in the spectrogram, if the bit is set to "1" it gets altered by a factor of alpha else it retains it's original value.

* Owing to the randomness of the bits that will be generated through the seed, we had to try a number of alpha values and radius values such that they wouldn't greatly distort the image and still be detectable through the checking process. After testing across upto 40 images, we settled at the previously mentioned values for the constants.

* An important observation we made was related to the alpha values - images that have high frequency areas work well with low values of alpha - sometimes even set to 1 - which is essentially doubling the value when the bit in the vector V is set to 1.

#### Check watermark

_./watermark 1.3 <\name of input image> <\any name> check N_

* Here we seed the random number generator with N in order to regenerate the binary vector of length "l".

* Tracing the "circle" using the same parameters as in the addition process, we collect all real values of the DFT of the image. Then, we run a Pearson's Correlation Coefficient (PCC) comparison between the two vectors.

* The threshold t has been arbitratily set at 0.2 to determine whether the watermark exists or not. Since the values of the PCC vary between -1 to +1 the only viable values to determine whether there's a match is between 0 to +1.


#### Qualitative Analysis
* As described earlier, there were images that suffered from noise/distortion in the watermarking process. On an overall basis, the implementation managed to identify 23 out of the 35 images correctly.

![](part1/markedImage3.png)
![](part1/markedImage1.png)

#### Quantitative Analysis
* Using the same paramters, we tested our test set of 35 images using N = 36. The results were as follows:

1. True positives = 23
2. False positives = 4
3. False negatives = 80
4. True negatives = 2

----------------------------------------------------------------------------------------------------------------------------------

### __Part 2 - Detecting Objects__

To run the application, run make command on part2. 
Then run _./detect \<filename of image\>._

This will create three files:
1. edges.png - Contains the output of the sobel edge detector.
2. detected.png - Contains Input Image with boundaries around the ICs
3. detected.txt - Description of boundaries in the required format.

#### Algorithm:

* The black and white image is take and sibel detector is applied to find the edge points.

![](part2/E_ex.png)

* Lines passing through these points are found using Hough Transform.

![](part2/lines.png)

* The corners on these lines are found using Harris corner detector.

![](part2/corners.png)

* All the rhombi formed with the corner points are considered.
* The image between the image is clipped and compared in the Fourier space with test images(present in repository) of IC's.
* The boundaries are decided based in the score from the comparison.

#### Assumptions:

* The algorithm assumes that the IC board with be parallel to the edges of the image.
* Only the vertical and horizontal lines are comsidered for corner detection.

#### Test Results:

* The algorithm struggle on IC's with a dimension of less that 32 px is one or more directions.
* Accuracy: The algorithm detects nearly all big square IC's and struggles a little on the rectangualr ones.

#### Future Improvements:

* Non Maximal Suppression of boundaries need to be implemented, edges and lines are implemented.
* More varied sample test images needs to be used for ideal accuracy.

