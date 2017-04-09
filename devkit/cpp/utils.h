#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <opencv2/core/version.hpp>
#if CV_MAJOR_VERSION == 2
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#elif CV_MAJOR_VERSION == 3
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#endif

bool imageFormat(std::string file_name, int type, int32_t width, int32_t height) {
    cv::Mat image = cv::imread(file_name, cv::IMREAD_UNCHANGED);
    if (image.type() != type) return false;
    if (image.cols != width)  return false;
    if (image.rows != height) return false;
    return true;
}

float statMean(std::vector< std::vector<float> > &errors, int32_t idx) {
    float err_mean = 0;
    for (int32_t i = 0; i < errors.size(); i++) {
        err_mean += errors[i][idx];
    }
    return err_mean / (float)errors.size();
}

float statWeightedMean(std::vector< std::vector<float> > &errors, int32_t idx, int32_t idx_num) {
    float err = 0;
    float num = 0;
    for (int32_t i = 0; i < errors.size(); i++) {
        err += errors[i][idx];
        num += errors[i][idx_num];
    }
    return err / std::max(num, 1.0f);
}

float statMin(std::vector< std::vector<float> > &errors, int32_t idx) {
    float err_min = 1;
    for (int32_t i = 0; i < errors.size(); i++) {
        if (errors[i][idx] < err_min) err_min = errors[i][idx];
    }
    return err_min;
}

float statMax(std::vector< std::vector<float> > &errors, int32_t idx) {
    float err_max = 0;
    for (int32_t i = 0; i < errors.size(); i++) {
        if (errors[i][idx] > err_max) err_max = errors[i][idx];
    }
    return err_max;
}

#endif // UTILS_H
