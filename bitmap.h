//
// Created by grego on 01/12/2020.
//

#ifndef ASSIGNMENT_2_AA2_BITMAP_H
#define ASSIGNMENT_2_AA2_BITMAP_H

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define BYTE_BOUND(value) value < 0 ? 0 : (value > 255 ? 255 : value)
#include "stb_image_write.h"
#include "stb_image.h"

#include <string>
#include <sstream>
#include <memory>


namespace bitmap {
    struct Pixel {
        uint8_t r;
        uint8_t g;
        uint8_t b;
    };

    struct BitmapImage {
        std::shared_ptr<uint8_t[]> data;
        int width;
        int height;
        int channels;

        BitmapImage(int w, int h) : width(w), height(h), channels(3) {}

        BitmapImage(const BitmapImage &img) : BitmapImage(img.width, img.height) {
            data = img.data;
        }

        explicit BitmapImage(const char *bitmap_filename) : data(
                stbi_load(bitmap_filename, &width, &height, &channels, 0), std::default_delete<uint8_t[]>()) {
            assert(channels == 3 && "The loaded image doesn't have 3 channels!");
        }

        void write(const char *filename) {
            stbi_write_bmp(filename, width, height, channels, data.get());
        }

        std::string str() const {
            std::stringstream ss;
            ss << "[image: channels(" << channels << "), width(" << width << "), height(" << height << ")]";
            return ss.str();
        }

        int positive_bound(int value, int max) {
            if (value < 0)
                return 0;
            if (value > max)
                return max;
            return value;
        }

        Pixel pixel_for_uv(float u, float v) {
            int x = positive_bound(int(u * (width - 1)), width);
            int y = positive_bound(int(v * (height - 1)), height);
            uint8_t red = data.get()[channels * x + channels * width * y + 0];
            uint8_t green = data.get()[channels * x + channels * width * y + 1];
            uint8_t blue = data.get()[channels * x + channels * width * y + 2];

            return Pixel({red, green, blue});
        }
    };

}


#endif //ASSIGNMENT_2_AA2_BITMAP_H
