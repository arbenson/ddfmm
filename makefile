ifndef HOST
  $(error "$${HOST} is not defined")
endif
include makeinc/${HOST}

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
          src/parallel.cpp

LIB_OBJ = $(LIB_SRC:.cpp=.o)

libwave.a: ${LIB_OBJ}
	$(AR) $(ARFLAGS) $@ $(LIB_OBJ)
	$(RANLIB) $@

tt: src/tt.o libwave.a 
	${CXX} -o $@ $^ ${LDFLAGS}

#------------------------------------------------------
clean:
	rm -rf *~ src/*.d src/*.o *.a tt 

tags:
	etags include/*.hpp src/*.cpp

