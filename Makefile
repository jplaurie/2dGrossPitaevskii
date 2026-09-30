BUILD_DIR ?= build/release
BUILD_TYPE ?= Release
JOBS ?= 4
CUDA_ARCHITECTURES ?= 120

.PHONY: all cpu mpi cuda initial test clean configure

all: configure
	cmake --build $(BUILD_DIR) -j$(JOBS)

configure:
	cmake -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) \
		-DGP2D_CUDA_ARCHITECTURES=$(CUDA_ARCHITECTURES)

cpu: configure
	cmake --build $(BUILD_DIR) --target gross_pitaevskii_cpu -j$(JOBS)

mpi: configure
	cmake --build $(BUILD_DIR) --target gross_pitaevskii_mpi -j$(JOBS)

cuda: configure
	cmake --build $(BUILD_DIR) --target gross_pitaevskii_cuda -j$(JOBS)

initial: configure
	cmake --build $(BUILD_DIR) --target gp2d_vortex_imprint gp2d_relax -j$(JOBS)

test: all
	ctest --test-dir $(BUILD_DIR) --output-on-failure

clean:
	cmake -E remove_directory build
