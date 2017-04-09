#ifndef IO_INTEGER_H
#define IO_INTEGER_H

#include <cstring>
#include <cstdint>
#include <string>
#include <opencv2/core/version.hpp>
#if CV_MAJOR_VERSION == 2
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#elif CV_MAJOR_VERSION == 3
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#endif

class IntegerImage {
public:
    // default constructor
    IntegerImage() {
        data_   = 0;
        width_  = 0;
        height_ = 0;
    }

    // construct image from png file
    IntegerImage(const std::string file_name) {
        readIntegerMap(file_name);
    }

    // copy constructor
    IntegerImage(const IntegerImage &I) {
        width_  = I.width_;
        height_ = I.height_;
        data_   = (int32_t*)malloc(width_*height_ * sizeof(int32_t));
        memcpy(data_, I.data_, width_*height_ * sizeof(int32_t));
    }

    // construct image from data
    IntegerImage(const int32_t* data, const int32_t width, const int32_t height) : width_(width), height_(height) {
        data_ = (int32_t*)malloc(width*height * sizeof(int32_t));
        memcpy(data_, data, width*height * sizeof(int32_t));
    }

    // construct empty (= all pixels 0) integer map of given width / height
    IntegerImage(const int32_t width, const int32_t height) : width_(width), height_(height) {
        data_ = (int32_t*)malloc(width*height * sizeof(int32_t));
        for (int32_t i = 0; i < width*height; i++) {
            data_[i] = 0;
        }
    }

    // deconstructor
    virtual ~IntegerImage() {
        if (data_) { free(data_); data_ = 0; }
    }

    // assignment operator, copies contents of D
    IntegerImage& operator= (const IntegerImage &I) {
        if (this != &I) {
            if (I.width_ != width_ || I.height_ != height_) {
                free(data_);
                width_  = I.width_;
                height_ = I.height_;
                data_   = (int32_t*)malloc(width_*height_ * sizeof(int32_t));
            }
            memcpy(data_, I.data_, width_*height_ * sizeof(int32_t));
        }
        return *this;
    }

    // read disparity image from png file
    void read(const std::string file_name) {
        if (data_) { free(data_); data_ = 0; }
        readIntegerMap(file_name);
    }

    // write disparity image to png file
    void write(const std::string file_name) {
        writeIntegerMap(file_name);
    }

    // get value at given pixel
    inline int32_t getValue(const int32_t u, const int32_t v) {
        return data_[v*width_ + u];
    }

    // set value at given pixel
    inline void setValue(const int32_t u, const int32_t v, const int32_t val) {
        data_[v*width_ + u] = val;
    }

    // direct access to private variables
    int32_t* data() { return data_; }
    int32_t  width() { return width_; }
    int32_t  height() { return height_; }

private:
    void readIntegerMap(const std::string file_name) {
        cv::Mat image = cv::imread(file_name, cv::IMREAD_UNCHANGED);
        width_  = image.cols;
        height_ = image.rows;
        data_   = (int32_t*)malloc(width_*height_ * sizeof(int32_t));
        for (int32_t v = 0; v < height_; v++) {
            for (int32_t u = 0; u < width_; u++) {
                setValue(u, v, image.at<uchar>(v, u));
            }
        }
    }

    void writeIntegerMap(const std::string file_name) {
        cv::Mat image(cv::Size(width_, height_), CV_16UC1, cv::Scalar(0));
        for (int32_t v = 0; v < height_; v++) {
            for (int32_t u = 0; u < width_; u++) {
                image.at<uchar>(v, u) = getValue(u, v);
            }
        }
        cv::imwrite(file_name, image);
    }

public:
    int32_t *data_;
    int32_t width_;
    int32_t height_;
};

#endif // INTEGER_IMAGE_H
