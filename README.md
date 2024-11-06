# Project Overview

An extinct spieces of [Linx](https://github.com/kabasset/Linx).

# Purpose

This project is a work in progress to gather ND array related features
of [EleFits](https://cnes.github.io/EleFits/)' `Raster` -- an extensible ND data storage class with pixel-wise operations --
and PhiFun's operations -- e.g. Fourier transforms, linear filtering, and interpolation.

The target is a header-only library focused on ease of use, extensibility,
and which interfaces seamlessly with the standard C++ library.
At some point, the ND array implementation of EleFits should be droped in favor of Linx.

The image processing feature set is expected to grow steadily on demand.
Currently supported operations are:

* Pixel-wise operations (`+`, `*`, `exp()`, `apply()`...) on `Raster`s and on `Raster` regions;
* Noise generation (`GaussianNoise`, `PoissonNoise`, `ImpulseNoise`...);
* Discrete Fourier transforms (FFTW wrapper) through the `DftPlan` class;
* Linear, rank-order, and morphological filtering through the `SimpleFilter`, `FilterSeq` and `FilterAgg` classes;
* Interpolation and extrapolation (`Nearest`, `Linear`, `Constant`...);
* Affine transformations.

# License

The Linx library is licensed under [Apache-2.0](LICENSE.txt).
Some executables including some unit tests are licensed under [GPL-3.0-or-later](https://www.gnu.org/licenses/gpl-3.0-standalone.html).

# Installation

The simplest way to use Linx, which is a header-only library, is to copy the deepest `Linx` folder into your include directory:

```sh
git clone https://github.com/kabasset/Linx.git
cp -r Linx/Linx/Linx /usr/include/
```

# Alternatives

The following libraries offer features similar to Linx.
Most of the libraries which support ND arrays are focused on linear algebra,
while most image processing libraries limit the number of axes.
ITK (and SimpleITK) is the notable exception.
Linx aims at being simpler, less verbose and more extensible,
although with a more limited feature set.

* Armadillo
* Blitz++
* Boost.MultiArray
* CImg
* Eigen
* ITK, SimpleITK
* ndarray
* OpenCV
* STL's valarray
* XTensor

Here is a quick comparison of ITK, CImg, Linx and NumPy/SciKit for the following use case:
read an image, dilate it with an L2-ball structuring element, and write the output.

ITK:

```cpp
using T = unsigned char;
static constexpr unsigned int N = 2;
using Image = itk::Image<T, N>;

auto raw = itk::ReadImage<ImageType>(input);

using StructuringElement = itk::FlatStructuringElement<N>;
StructuringElement::RadiusType strelRadius;
strelRadius.Fill(radius);
StructuringElementType ball = StructuringElement::Ball(strelRadius);
using GrayscaleDilateImageFilter = itk::GrayscaleDilateImageFilter<Image, Image, StructuringElement>;
GrayscaleDilateImageFilter::Pointer dilateFilter = GrayscaleDilateImageFilter::New();
dilateFilter->SetInput(input);
dilateFilter->SetKernel(ball);

itk::WriteImage(dilateFilter->GetOutput(), output);
```

CImg (limited to N <= 3):

```cpp
using T = unsigned char;

auto raw = cimg::CImg<T>().load(input);

cimg::CImg<bool> ball(2 * radius + 1, 2 * radius + 1, 1, 1, false);
bool color[1] = {true};
ball.draw_circle(radius, radius, radius, color);
auto dilated = raw.get_dilate(ball, 0, true);

dilated.write(output);
```

Linx:


```cpp
using T = unsigned char;
static constexpr Linx::Index N = 2;

auto raw = Linx::read<T, N>(input);

auto ball = Linx::Mask<N>::ball<2>(radius);
auto dilated = dilation<T>(ball) * extrapolation(raw);

Linx::write(dilated, output);
```

NumPy/SciKit:

```python
raw = np.load(input)

ball = skimage.morphology.disk(radius)
dilated = skimage.morphology.dilation(raw, ball)

np.save(output, dilated)
```
