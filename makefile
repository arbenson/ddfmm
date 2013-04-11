ifndef HOST
  $(error "$${HOST} is not defined")
endif
include makeinc/${HOST}

# default rule
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -Iinclude -c $< -o $@

LIB_SRC = src/wave3d.cpp \
          src/kernel3d.cpp \
          src/mlib3d.cpp \
          src/wave3d_setup.cpp \
          src/wave3d_eval.cpp \
          src/wave3d_check.cpp \
          src/vecmatop.cpp \
          src/parallel.cpp

LIB_OBJ = $(LIB_SRC:.cpp=.o)

TST_SRC = src/tt.cpp src/ss.cpp

DEP = $(LIB_SRC:.cpp=.d) $(TST_SRC:.cpp=.d)

libwave.a: ${LIB_OBJ}
	$(AR) $(ARFLAGS) $@ $(LIB_OBJ)
	$(RANLIB) $@

tt: src/tt.o libwave.a 
	${CXX} -o $@ $^ ${LDFLAGS}

ss: src/ss.o
	${CXX} -o $@ $^ ${LDFLAGS}

-include $(DEP)

#------------------------------------------------------
tilde:
	rm -f *~

clean:
	rm -rf *~ src/*.d src/*.o *.a tt 

# TODO (Austin): not sure if these should be included in the default clean or not
cleandata:
	rm -rf data/*.bin data/*.wrl_*_*_*

tags:
	etags include/*.hpp src/*.cpp

