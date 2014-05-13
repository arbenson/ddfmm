# Distributed Directional Fast Multipole Method
#   Copyright (C) 2014 Austin Benson, Lexing Ying, and Jack Poulson
#
# This file is part of DDFMM.
#
#    DDFMM is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    DDFMM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with DDFMM.  If not, see <http://www.gnu.org/licenses/>.
HOST=stampede
ifndef HOST
  $(error "$${HOST} is not defined")
endif
include makeinc/${HOST}
CXXFLAGS += -std=c++11

# default rule
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -Iinclude -c $< -o $@

LIB_SRC = src/wave3d.cpp \
          src/kernel3d.cpp \
          src/mlib3d.cpp \
          src/wave3d_setup.cpp \
          src/wave3d_eval.cpp \
          src/wave3d_check.cpp \
          src/vecmatop.cpp \
          src/parallel.cpp \
          src/global.cpp \
          src/utility.cpp \
          src/translations.cpp \
          src/file_io.cpp \
          src/communication.cpp \
          src/DataCollection.cpp \
          src/data_distrib.cpp \

ACOU_SRC = examples/acoustic/trmesh.cpp \
           examples/acoustic/acou3d.cpp

ACOU_OBJ = $(ACOU_SRC:.cpp=.o)
LIB_OBJ = $(LIB_SRC:.cpp=.o)

libwave.a: ${LIB_OBJ}
	$(AR) $(ARFLAGS) $@ $(LIB_OBJ)
	$(RANLIB) $@

libacoustic.a: ${ACOU_OBJ}
	$(AR) $(ARFLAGS) $@ $(ACOU_OBJ)
	$(RANLIB) $@

tt: src/tt.o libwave.a
	${CXX} -o $@ $^ ${LDFLAGS}

acoustic: examples/acoustic/main.o libwave.a libacoustic.a
	${CXX} -o $@ $^ ${LDFLAGS}

file_io_test: src/file_io_test.o libwave.a
	${CXX} -o $@ $^ ${LDFLAGS}

#------------------------------------------------------
clean:
	rm -rf *~ src/*.d src/*.o *.a tt

tags:
	etags include/*.hpp src/*.cpp

