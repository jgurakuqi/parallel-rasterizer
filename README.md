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

stb_image.h and stb_image_write.h are a freely available library I decided to import to load easily bitmap files. If a better alternative is available they could be replaced.
Also the multi-threading level could be improved, to achieve even better performance.

```bash
git clone https://github.com/jgurakuqi/parallel-rasterizer
```

## License

MIT License

Copyright (c) 2021 Jurgen Gurakuqi

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS," WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
