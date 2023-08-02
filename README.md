# parallel-rasterizer


## Table of Contents

- [Requirements](#Requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)


## Requirements

The goal of this project is to implement a parallel 3D software rendering pipeline with programmable fragment shader.

## Installation

In order to run this project a basic C++ compiler with a bash/shell or an IDE is required.
Also, the provided input images are required in a bitmap format (.bmp).

## Usage

To run the projects it's required an IDE suitable to run the Jupyter Notebook.

## Contributing

I implemented the code exploiting the problem-specific characteristics (e.g., semi-constant lightning, same shadow on the background fabric, ...), hence using image processing techniques specific wrt to it, at least for the custom segmentation, hence in order to readapt the custom segmentation algorithm many modifications might be required. On the other hand, the Watershed-based segmentation, followed a more generic approach, and it's easier to readapt to other problems.

```bash
git clone https://github.com/jgurakuqi/parallel-rasterizer
```

## License

MIT License

Copyright (c) 2021 Jurgen Gurakuqi

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS," WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
