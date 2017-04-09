#ifndef IO_FLOW_H
#define IO_FLOW_H

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstring>
#include <cstdint>
#include <algorithm>
#include <opencv2/core/version.hpp>
#if CV_MAJOR_VERSION == 2
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#elif CV_MAJOR_VERSION == 3
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#endif

#include "log_colormap.h"

class FlowImage {
public:
    // default constructor
    FlowImage() {
        data_   = 0;
        width_  = 0;
        height_ = 0;
    }

    // construct flow image from png file
    FlowImage(const std::string file_name) {
        readFlowField(file_name);
    }

    // copy constructor
    FlowImage(const FlowImage &F) {
        width_  = F.width_;
        height_ = F.height_;
        data_   = (float*)malloc(width_*height_ * 3 * sizeof(float));
        memcpy(data_, F.data_, width_*height_ * 3 * sizeof(float));
    }

    // construct flow field from data
    FlowImage(const float* data, const int32_t width, const int32_t height) : width_(width), height_(height) {
        data_ = (float*)malloc(width*height * 3 * sizeof(float));
        memcpy(data_, data, width*height * 3 * sizeof(float));
    }

    // construct empty (= all pixels invalid) flow field of given width / height
    FlowImage(const int32_t width, const int32_t height) : width_(width), height_(height) {
        data_ = (float*)malloc(width*height * 3 * sizeof(float));
        for (int32_t i = 0; i < width*height * 3; i++) {
            data_[i] = 0;
        }
    }

    // deconstructor
    virtual ~FlowImage() {
        if (data_) { free(data_); data_ = 0; }
    }

    // assignment operator, copies contents of F
    FlowImage& operator= (const FlowImage &F) {
        if (this != &F) {
            if (F.width_ != width_ || F.height_ != height_) {
                free(data_);
                width_  = F.width_;
                height_ = F.height_;
                data_   = (float*)malloc(width_*height_ * 3 * sizeof(float));
            }
            memcpy(data_, F.data_, width_*height_ * 3 * sizeof(float));
        }
        return *this;
    }

    // read flow field from png file
    void read(const std::string file_name) {
        if (data_) { free(data_); data_ = 0; }
        readFlowField(file_name);
    }

    // write flow field to png file
    void write(const std::string file_name) {
        writeFlowField(file_name);
    }

    // write flow field to false color map using the Middlebury colormap
    void writeColor(const std::string file_name, float max_flow = -1.0f) {
        if (max_flow <= 1.0f) {
            max_flow = std::max(maxFlow(), 1.0f);
        }
        writeFalseColors(file_name, max_flow);
    }

    // get optical flow u-component at given pixel
    inline float getFlowU(const int32_t u, const int32_t v) {
        return data_[3 * (v*width_ + u) + 0];
    }

    // get optical flow v-component at given pixel
    inline float getFlowV(const int32_t u, const int32_t v) {
        return data_[3 * (v*width_ + u) + 1];
    }

    // check if optical flow at given pixel is valid
    inline bool isValid(const int32_t u, const int32_t v) {
        return data_[3 * (v*width_ + u) + 2] > 0.5;
    }

    // get optical flow magnitude at given pixel 
    inline float getFlowMagnitude(const int32_t u, const int32_t v) {
        float fu = getFlowU(u, v);
        float fv = getFlowV(u, v);
        return sqrt(fu*fu + fv*fv);
    }

    // set optical flow u-component at given pixel
    inline void setFlowU(const int32_t u, const int32_t v, const float val) {
        data_[3 * (v*width_ + u) + 0] = val;
    }

    // set optical flow v-component at given pixel
    inline void setFlowV(const int32_t u, const int32_t v, const float val) {
        data_[3 * (v*width_ + u) + 1] = val;
    }

    // set optical flow at given pixel to valid / invalid
    inline void setValid(const int32_t u, const int32_t v, const bool valid) {
        data_[3 * (v*width_ + u) + 2] = valid ? 1.0f : 0.0f;
    }

    // get maximal optical flow magnitude
    float maxFlow() {
        float max_flow = 0;
        for (int32_t u = 0; u < width_; u++) {
            for (int32_t v = 0; v < height_; v++) {
                if (isValid(u, v) && getFlowMagnitude(u, v) > max_flow) {
                    max_flow = getFlowMagnitude(u, v);
                }
            }
        }
        return max_flow;
    }

    // interpolate all missing (=invalid) optical flow vectors
    void interpolateBackground() {
        // for each row do
        for (int32_t v = 0; v < height_; v++) {
            // init counter
            int32_t count = 0;

            // for each pixel do
            for (int32_t u = 0; u < width_; u++) {
                // if flow is valid
                if (isValid(u, v)) {
                    // at least one pixel requires interpolation
                    if (count >= 1) {

                        // first and last value for interpolation
                        int32_t u1 = u - count;
                        int32_t u2 = u - 1;

                        // set pixel to min flow
                        if (u1 > 0 && u2 < width_ - 1) {
                            float fu_ipol = std::min(getFlowU(u1 - 1, v), getFlowU(u2 + 1, v));
                            float fv_ipol = std::min(getFlowV(u1 - 1, v), getFlowV(u2 + 1, v));
                            for (int32_t u_curr = u1; u_curr <= u2; u_curr++) {
                                setFlowU(u_curr, v, fu_ipol);
                                setFlowV(u_curr, v, fv_ipol);
                                setValid(u_curr, v, true);
                            }
                        }
                    }
                    // reset counter
                    count = 0;

                    // otherwise increment counter
                }
                else {
                    count++;
                }
            }

            // extrapolate to the left
            for (int32_t u = 0; u < width_; u++) {
                if (isValid(u, v)) {
                    for (int32_t u2 = 0; u2 < u; u2++) {
                        setFlowU(u2, v, getFlowU(u, v));
                        setFlowV(u2, v, getFlowV(u, v));
                        setValid(u2, v, true);
                    }
                    break;
                }
            }

            // extrapolate to the right
            for (int32_t u = width_ - 1; u >= 0; u--) {
                if (isValid(u, v)) {
                    for (int32_t u2 = u + 1; u2 <= width_ - 1; u2++) {
                        setFlowU(u2, v, getFlowU(u, v));
                        setFlowV(u2, v, getFlowV(u, v));
                        setValid(u2, v, true);
                    }
                    break;
                }
            }
        }

        // for each column do
        for (int32_t u = 0; u < width_; u++) {
            // extrapolate to the top
            for (int32_t v = 0; v < height_; v++) {
                if (isValid(u, v)) {
                    for (int32_t v2 = 0; v2 < v; v2++) {
                        setFlowU(u, v2, getFlowU(u, v));
                        setFlowV(u, v2, getFlowV(u, v));
                        setValid(u, v2, true);
                    }
                    break;
                }
            }

            // extrapolate to the bottom
            for (int32_t v = height_ - 1; v >= 0; v--) {
                if (isValid(u, v)) {
                    for (int32_t v2 = v + 1; v2 <= height_ - 1; v2++) {
                        setFlowU(u, v2, getFlowU(u, v));
                        setFlowV(u, v2, getFlowV(u, v));
                        setValid(u, v2, true);
                    }
                    break;
                }
            }
        }
    }

    // compute error map of flow field, given the non-occluded and occluded
    // ground truth optical flow maps. stores result as color png image.
    cv::Mat errorImage(FlowImage &F_noc, FlowImage &F_occ, bool log_colors = false) {
        cv::Mat image(cv::Size(width_, height_), CV_8UC3, cv::Scalar(0, 0, 0));
        for (int32_t v = 1; v < height() - 1; v++) {
            for (int32_t u = 1; u < width() - 1; u++) {
                if (F_occ.isValid(u, v)) {
                    cv::Vec3b val;
                    if (log_colors) {
                        float dfu = getFlowU(u, v) - F_occ.getFlowU(u, v);
                        float dfv = getFlowV(u, v) - F_occ.getFlowV(u, v);
                        float f_err = sqrt(dfu*dfu + dfv*dfv);
                        float f_mag = F_occ.getFlowMagnitude(u, v);
                        float n_err = std::min(f_err / 3.0f, 20.0f*f_err / f_mag);
                        for (int32_t i = 0; i < 10; i++) {
                            if (n_err >= LC[i][0] && n_err < LC[i][1]) {
                                val[2] = (uint8_t)LC[i][2]; //red
                                val[1] = (uint8_t)LC[i][3]; //green
                                val[0] = (uint8_t)LC[i][4]; //blue
                            }
                        }
                        if (!F_noc.isValid(u, v)) {
                            val[2] *= 0.5; //red
                            val[1] *= 0.5; //green
                            val[0] *= 0.5; //blue
                        }
                    }
                    else {
                        float dfu = getFlowU(u, v) - F_occ.getFlowU(u, v);
                        float dfv = getFlowV(u, v) - F_occ.getFlowV(u, v);
                        float f_err = std::min(sqrtf(dfu*dfu + dfv*dfv), 5.0f) / 5.0f;

                        val[2] = (uint8_t)(f_err*255.0); //red
                        val[1] = (uint8_t)(f_err*255.0); //green
                        val[0] = (uint8_t)(f_err*255.0); //blue

                        if (!F_noc.isValid(u, v)) {
                            val[1] = 0; //green
                            val[0] = 0; //blue
                        }
                    }
                    for (int32_t v2 = v - 1; v2 <= v + 1; v2++) {
                        for (int32_t u2 = u - 1; u2 <= u + 1; u2++) {
                            image.at<cv::Vec3b>(v2, u2) = val;
                        }
                    }
                }
            }
        }
        return image;
    }

    // direct access to private variables
    float*  data() { return data_; }
    int32_t width() { return width_; }
    int32_t height() { return height_; }

private:
    void readFlowField(const std::string file_name) {
        cv::Mat image = cv::imread(file_name, cv::IMREAD_UNCHANGED);
        width_  = image.cols;
        height_ = image.rows;
        data_   = (float*)malloc(width_*height_ * 3 * sizeof(float));

        for (int32_t v = 0; v < height_; v++) {
            for (int32_t u = 0; u < width_; u++) {
                cv::Vec3w val = image.at<cv::Vec3w>(v, u);
                if (val[0] > 0) {
                    setFlowU(u, v, ((float)val[2] - 32768.0f) / 64.0f);
                    setFlowV(u, v, ((float)val[1] - 32768.0f) / 64.0f);
                    setValid(u, v, true);
                }
                else {
                    setFlowU(u, v, 0);
                    setFlowV(u, v, 0);
                    setValid(u, v, false);
                }
            }
        }
    }

    void writeFlowField(const std::string file_name) {
        cv::Mat image(cv::Size(width_, height_), CV_16UC3, cv::Scalar(0, 0, 0));
        for (int32_t v = 0; v < height_; v++) {
            for (int32_t u = 0; u < width_; u++) {
                cv::Vec3w val;
                val[0] = 0; //brue
                val[1] = 0; //green
                val[2] = 0; //red

                if (isValid(u, v)) {
                    val[2] = (uint16_t)std::max(std::min(getFlowU(u, v)*64.0f + 32768.0f, 65535.0f), 0.0f);
                    val[1] = (uint16_t)std::max(std::min(getFlowV(u, v)*64.0f + 32768.0f, 65535.0f), 0.0f);
                    val[0] = 1;
                }
                image.at<cv::Vec3w>(v, u) = val;
            }
        }
        cv::imwrite(file_name, image);
    }

    inline void hsvToRgb(float h, float s, float v, float &r, float &g, float &b) {
        float c = v*s;
        float h2 = 6.0f*h;
        float x = c*(1.0f - fabs(fmod(h2, 2.0f) - 1.0f));
        if (0 <= h2&&h2 < 1) { r = c; g = x; b = 0; }
        else if (1 <= h2&&h2 < 2) { r = x; g = c; b = 0; }
        else if (2 <= h2&&h2 < 3) { r = 0; g = c; b = x; }
        else if (3 <= h2&&h2 < 4) { r = 0; g = x; b = c; }
        else if (4 <= h2&&h2 < 5) { r = x; g = 0; b = c; }
        else if (5 <= h2&&h2 <= 6) { r = c; g = 0; b = x; }
        else if (h2 > 6) { r = 1; g = 0; b = 0; }
        else if (h2 < 0) { r = 0; g = 1; b = 0; }
    }

    void writeFalseColors(const std::string file_name, const float max_flow) {
        float n = 8; // multiplier
        cv::Mat image(cv::Size(width_, height_), CV_8UC3, cv::Scalar(0, 0, 0));
        for (int32_t v = 0; v < height_; v++) {
            for (int32_t u = 0; u < width_; u++) {
                float r = 0, g = 0, b = 0;
                if (isValid(u, v)) {
                    float mag = getFlowMagnitude(u, v);
                    float dir = atan2(getFlowV(u, v), getFlowU(u, v));
                    float h = fmod(dir / (2.0f*M_PI) + 1.0f, 1.0f);
                    float s = std::min(std::max(mag*n / max_flow, 0.0f), 1.0f);
                    float v = std::min(std::max(n - s, 0.0f), 1.0f);
                    hsvToRgb(h, s, v, r, g, b);
                }
                image.at<cv::Vec3b>(v, u) = cv::Vec3b(b*255.0f, g*255.0f, r*255.0f);
            }
        }
        cv::imwrite(file_name, image);
    }

public:
    float  *data_;
    int32_t width_;
    int32_t height_;
};

#endif // IO_FLOW_H
