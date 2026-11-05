// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _data_management_kernels_h
#define _data_management_kernels_h


#include <vector>

namespace bear{


template <typename T> void allocateOnDevice(T*& device_data, size_t nb_double_values);
template <typename T> void moveToDevice(T*& device_data, std::vector<T>& host_data, const bool alloc_memory);
template <typename T> void moveToDevice(T*& device_data, std::vector<T>& host_data);
void moveToDevice(float*& device_data, std::vector<double>& host_data);
template <typename T> void moveToHost(T*& device_data, std::vector<T>& host_data);
template <typename T> void moveToHostAndDelete(T*& device_data, std::vector<T>& host_data);
template <typename T> void deleteFromDevice(T*& device_data);
template <typename T> void initializeOnDevice(T*& device_data, const size_t nb_points);
template <typename T> void copyOnDevice(T* dst, const T* src, size_t nb_elements);



}


#endif
