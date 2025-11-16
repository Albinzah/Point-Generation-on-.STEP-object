# ==== Compiler setup ====
CXX = g++
CXXWARN = -Wall -Wextra
STD = -std=c++20
INCLUDE_DIRS = -Iinclude

# ==== OpenCascade paths (local) ====
OCC_DIR = $(HOME)/projects/OCCT/install
OCC_INC = $(OCC_DIR)/include/opencascade
OCC_LIB = $(OCC_DIR)/lib

# ==== OpenCascade paths (conda) ====
OCC_DIR = $(CONDA_PREFIX)
OCC_INC = $(OCC_DIR)/include/opencascade
OCC_LIB = $(OCC_DIR)/lib

# ==== OpenCascade libs (local) ====
OCC_LIBS = TKernel TKMath TKBRep TKGeomBase TKGeomAlgo TKPrim TKTopAlgo \
           TKG2d TKG3d TKService TKMesh TKHLR TKV3d TKOpenGl \
           TKCDF TKCAF TKXCAF TKDESTEP TKXSBase \
           TKBO TKShHealing TKLCAF TKVCAF TKDE

# ==== OpenCascade libs (conda) ====
OCC_LIBS = TKernel TKMath TKBRep TKGeomBase TKGeomAlgo TKPrim TKTopAlgo \
           TKG2d TKG3d TKService TKMesh TKHLR TKV3d TKOpenGl \
           TKCDF TKCAF TKXCAF TKSTEP TKSTEPAttr TKSTEPBase TKIGES TKXSBase \
           TKBO TKShHealing TKLCAF TKVCAF TKXDESTEP TKXDEIGES TKXDE

LDFLAGS = $(foreach lib,$(OCC_LIBS),-l$(lib)) -L$(OCC_LIB) -lX11 -fopenmp

# ==== Leak suppression options ====
LEAKSUPP = ASAN_OPTIONS=detect_leaks=1:fast_unwind_on_malloc=0 \
           LSAN_OPTIONS=suppressions=suppress.txt:report_objects=1

# ==== PyBind11 ====
PYBIND11_INCLUDES := $(shell python3 -m pybind11 --includes)

# ==== Project setup ====
TARGET = main
SRCS = main.cpp src/util.cpp src/types.cpp src/geometryUtils.cpp src/debugVars.cpp
OBJS = $(SRCS:.cpp=.o)

# ==== Default flags ====
DEBUG_FLAGS = $(CXXWARN) $(STD) $(INCLUDE_DIRS) -I$(OCC_INC) -fsanitize=address -g
RELEASE_FLAGS = $(CXXWARN) $(STD) $(INCLUDE_DIRS) -I$(OCC_INC) -O3 -DNDEBUG

# ==== Default target ====
all: debug

# ==== Build types ====
debug: CXXFLAGS = $(DEBUG_FLAGS) -fopenmp
debug: $(TARGET)

release: CXXFLAGS = $(RELEASE_FLAGS)
release: $(TARGET)

# ==== Target rule ====
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ==== Bindings ====
bindings: clean_bindings
	$(CXX) $(RELEASE_FLAGS) $(PYBIND11_INCLUDES) -shared -fPIC \
	src/geometryUtils.cpp src/types.cpp src/util.cpp tests/bindings.cpp \
	-o geomlib/geom_cpp$(shell python3-config --extension-suffix) \
	$(LDFLAGS)

bindings_debug:
	$(CXX) $(DEBUG_FLAGS) $(PYBIND11_INCLUDES) -shared -fPIC \
	src/geometryUtils.cpp src/types.cpp src/util.cpp tests/bindings.cpp \
	-o geomlib/geom_cpp_debug$(shell python3-config --extension-suffix) \
	$(LDFLAGS)

# ==== Run ====
run: debug
	LD_LIBRARY_PATH=$(OCC_LIB):$$LD_LIBRARY_PATH $(LEAKSUPP) ./$(TARGET)

# ==== Cleanup ====
clean:
	rm -f $(OBJS) $(TARGET) geom_cpp*.so

clean_bindings:
	rm -f geom_cpp*.so
