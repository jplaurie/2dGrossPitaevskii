#include "hdf5_io.hpp"

#include <hdf5.h>

#include <chrono>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
template <class T> void check(T status) {
    if (status < 0)
        throw std::runtime_error("could not construct HDF5 metadata fixture");
}

void expectFailure(const std::filesystem::path &path, const std::string &message) {
    try {
        (void)readHdf5Field(path);
    } catch (const std::runtime_error &error) {
        if (std::string(error.what()).find(message) != std::string::npos)
            return;
        throw;
    }
    throw std::runtime_error("malformed HDF5 attribute was accepted");
}

void replaceAttribute(const std::filesystem::path &path, const char *name, hid_t type, hid_t space,
                      const void *value) {
    const hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    check(file);
    check(H5Adelete(file, name));
    const hid_t attribute = H5Acreate2(file, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
    check(attribute);
    check(H5Awrite(attribute, type, value));
    check(H5Aclose(attribute));
    check(H5Fclose(file));
}

void runChecks(const std::filesystem::path &path) {
    Parameters parameters;
    parameters.nx = parameters.ny = 8;
    const std::vector<Complex> physical(64, Complex(0.4, 0.2));
    const auto reset = [&] { writeHdf5Field(path, parameters, 0.0, 0, physical); };
    reset();
    if (readHdf5Field(path).wavefunction != physical)
        throw std::runtime_error("valid HDF5 field did not round trip");

    const hsize_t count = 256;
    const hid_t arraySpace = H5Screate_simple(1, &count, nullptr);
    check(arraySpace);
    const std::vector<unsigned long long> dimensions(count, 8);
    replaceAttribute(path, "nx", H5T_NATIVE_ULLONG, arraySpace, dimensions.data());
    expectFailure(path, "attribute must be scalar: nx");
    reset();

    const hid_t stringType = H5Tcopy(H5T_C_S1);
    check(stringType);
    check(H5Tset_size(stringType, 32));
    std::vector<char> strings(count * 32, '\0');
    replaceAttribute(path, "format", stringType, arraySpace, strings.data());
    expectFailure(path, "attribute must be scalar: format");
    reset();

    const hid_t scalarSpace = H5Screate(H5S_SCALAR);
    check(scalarSpace);
    const char text[32] = "8";
    replaceAttribute(path, "nx", stringType, scalarSpace, text);
    expectFailure(path, "attribute must be numeric: nx");
    reset();

    const unsigned long long number = 8;
    replaceAttribute(path, "format", H5T_NATIVE_ULLONG, scalarSpace, &number);
    expectFailure(path, "attribute must be a fixed-length string: format");
    reset();

    check(H5Tset_size(stringType, H5T_VARIABLE));
    const char *variableText = "gp2d_wavefunction_v1";
    replaceAttribute(path, "format", stringType, scalarSpace, &variableText);
    expectFailure(path, "attribute must be a fixed-length string: format");

    check(H5Tclose(stringType));
    check(H5Sclose(scalarSpace));
    check(H5Sclose(arraySpace));
}
} // namespace

int main() {
    const auto suffix = std::chrono::steady_clock::now().time_since_epoch().count();
    const auto path = std::filesystem::temp_directory_path() /
                      ("gp2d_hdf5_validation_" + std::to_string(suffix) + ".h5");
    int status = 0;
    try {
        runChecks(path);
        std::cout << "HDF5 attribute shapes and types are validated before reading\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        status = 1;
    }
    std::error_code ignored;
    std::filesystem::remove(path, ignored);
    return status;
}
