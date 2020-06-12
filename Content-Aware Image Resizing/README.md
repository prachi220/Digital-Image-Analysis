This involves the implementation of [“Seam carving for content-aware image resizing.” ACM Transactions on graphics (TOG). Vol. 26. No. 3. ACM, 2007](http://graphics.cs.cmu.edu/courses/15-463/2012_fall/hw/proj3-seamcarving/imret.pdf). This focuses on changing the size of an image by gracefully carving-out or inserting pixels in different parts of the image. Seam carving uses
an energy function defining the importance of pixels. A seam is a connected path of low energy pixels crossing the image from top to bottom, or from left to right.  By repeatedly carving out or inserting seams in one direction we can change the aspect ratio of an image.